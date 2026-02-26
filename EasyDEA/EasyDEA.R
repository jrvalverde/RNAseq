#!/usr/bin/env Rscript
#
#
# We start from data that has already been cleaned by performing a quality
# check with FastQC and subsequent edge trimming.
#
# NOTE: for debugging look into Rbase_tools.R

# we need this function here so we can include additional files
# until we make this into a package
getScriptPath <- function()
{
     
    # this works if we were called with 'source()'
    src.path <- utils::getSrcFilename(function() {}, full.names=T)
    if ((! is.null(src.path)) && 
        (length(src.path) != 0)) {
        return(src.path)
	}
    # this works for Rscript, it may match more than one path
    # if Rscript was used with several --file= arguments, in
    # which case we will return *only* the first one
    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
    script.path <- regmatches(cmd.args, m)
    if (length(script.path) >= 1) return(script.path[1])
    if (length(script.path) == 1) return(script.path)	# may return multiple matches

    # this works for 'R -f', it will return only the first -f argument
    for (i in 1:length(cmd.args) ) {
        print(i); print(cmd.args[i])
        if (cmd.args[i] == '-f') return(cmd.args[i+1])
    }
    
    # if we arrive here, we didn't match anything, turn to last resort
    return (sys.frame(1)$ofile)
}

# get my location
my.source <- getScriptPath()
#my.dir <- '.'
my.name <- basename(my.source)
my.dir <- normalizePath(dirname(my.source))
my.lib <- normalizePath(file.path(my.dir, "lib"))
# inside my location there should be a 'lib' directory with the needed
# auxiliary scripts: we'll source them all (after loading the sourceDir
# function from lib/Rbase_tools.R)
source(file.path(my.lib, "Rbase_tools.R"))
sourceDir(my.lib, VERBOSE = F)

# ensure we have all the packages we need
# we do it here to ensure all is installed before we are run and that we do
# not get an obscure comment somewhere while being run that might escape
# out attention among the flow of messages printed
use.package("optparse")
use.package("ensembldb")		# needs to be first to avoid S4 inconsistencies
use.package("Rsubread")		# for read mapping
use.package("AnnotationHub")		# to seek annotation packages
use.package("AnnotationForge")	# to build our annotation package
use.package("GO.db")
use.package("PFAM.db")

use.package("biomaRt")		# an alternate approach to retrieve annotation
use.package("STRINGdb")
use.package("igraph")
#use.package("netSmooth")

use.package("tibble")			# general tools
use.package("tidyr")
use.package("stringr")
use.package("dplyr")
use.package("readr")
use.package("keypress")

use.package("edgeR")			# RNAseq with edgeR
use.package("limma")
use.package("RColorBrewer")
use.package("gplots")

use.package("DESeq2")			# RNAseq with DESeq2
use.package("IHW")			# for p-value adjustment with IHW
use.package("ggplot2")

use.package("cluster") 
use.package("factoextra")
use.package("fpc")
use.package("NbClust")
use.package("clusterProfiler")
use.package("enrichplot")
use.package("pathview")

use.package("tcltk")
use.package("gWidgets2")


# Parser arguments
options <- get.options()

# For convenience, we will assign options to specific names
#	we could get a similar efect if we were to use attach()
#	but this makes it evident which variables correspond to 
# 	configurable options, whereas with attach() they would
#   seemingly pop up out of nowhere as if by magic and it
#   wouldn't be evident were did they get their value from
ALIGN <-                  options$ALIGN         # we need to align the reads
PAIRED <-                 options$PAIRED        # we have paired reads
BOTH <-					  options$BOTH          # require BOTH reads in a pair to match
USE.ONLINE.ANNOTATION <-  options$USE.ONLINE.ANNOTATION
USE.EDGER <-              options$USE.EDGER
USE.DESEQ2 <-             options$USE.DESEQ2
reference <-              options$reference		# reference genome sequence
annotation <-             options$annotation	# reference annotation file
release <-                options$release		# release name
target.organism <-        options$target.organism
ens.version <-            options$ens.version	# version of ENSEMBL to use
mart.name <-              options$mart.name     # mart from BiomaRt to use
org.package <-            options$org.package   # name of the Org package for this organism
ncbi.taxid <-             options$ncbi.taxid
kegg.organism <-		  options$kegg.organism # name of the organism in KEGG
n.genes <-                options$n.genes		# number of top genes to revise
fastq.dir <-              options$fastq.dir     # directory with fastq files
alignment.dir <-          options$alignment.dir # guess what
feature.count.dir <-	  options$feature.count.dir
rnaseq.out <-             options$rnaseq.out    # directory for the results
my.name <-                options$my.name		# used to identify maintainer of
my.email <-               options$my.email		# any created annotation package
my.user <-                options$my.user       # used to access MySQL
my.password <-            options$my.password   # used to access MySQL
metadata <-               options$metadata      # name of the metadata file
cpm.threshold <-          options$cpm.threshold # you may need to do a test run
significance.threshold <- options$significance.threshold # to determine these
design.column <-          options$design.column # column in the metadata used to guide the comparisons
config.file <-            options$config.file   # not really needed or used later
INTERACTIVE <-            options$INTERACTIVE   # whether we want to control the run
VERBOSE <-                options$VERBOSE       # produce additional output
# or we could simply attach(options)
# only I (JR) am not too happy with non-explicit variable names

# convenience variables
by.rows=1
by.columns=2

##############################################################################
#
# DO THE WORK
#
##############################################################################

short.title('EasyDEA: RNA-Seq Analysis')		# Print a visible title

# Generate Output Directory Hierarchy
# folder is the working folder, rnaseq.out/{both_ends|any_end}
folder <- create.output.hierarchy(rnaseq.out, use.both.reads=BOTH)

# To ensure reproducibility, make a copy of the metadata
# into the output rnaseq folder
system(paste("cp", metadata, rnaseq.out))
system(paste("cp", metadata, folder))
# also save our own code
if (! file.exists(file.path(folder, my.name))) {
    system(paste("cp -R", my.source, my.lib, rnaseq.out))
}

# save the options so we keep a record of how the analysis was generated
opt.file <- file.path(folder, 'RNAseq_options.conf')
if (!save.options(options, file=opt.file)) {
    cat.warn("Could not save options file", opt.file)
}
# keep a log file for tracking and reporting
logfile <- file.path(folder, "log", "RNAseq.log")
log <- openLogFile(logfile)	# opens with sink(), may be closed
						    # specifically with closeLogFile(log)
                            
# when a script ends or when stop() is called, there are a number
# of housekeeping tasks to do. on.exit() allows us to add additional
# tasks so that they, too, are executed at the end of the script;
# adding this here we do not need to keep track of all the sink()s
on.exit(sink.titanic, add=T)



##############################################################################
#
#  PREPARE ANNOTATION SO IT IS READY WHEN WE GET THE RESULTS
#
##############################################################################

annotation.dir <- file.path(rnaseq.out, "annotation")

#---------------------------------------------------------------------------

# Now is time to prepare to connect all the results we get with the existing
# information  we know from the literature, so we will retrieve infos from
# ENSEMBL and connect with the genes  we have kept. We are interested only in
# genes and transcriptomes (non-characterized genes) from the organism
# Gallus galllus. there is no database available for this organism so we
# have to create it by ourselves

# we create two different databases one from the gtf file (small archive)
# and one from the data we extracted from ENSEMBL and stored it in a sqlite 
# file

# we need to load the package that enables us to connect a file as an
# external database ensembldb. After that we create a variable for the DB
# storing it in the memory simple process: we extract the files from ENSEMBL
# and store the data to MySQL account and then with the script
# 'generate-EnsDBs.R' we create a sqlite file that has all the information
# needed to continue

cat('\n\n')     # leave some room for clarity
short.title('Ensembl Annotation')

# THIS SHOULD MATCH THE LOCATION IN CREATE.OUTPUT.HIERARCHY
#       XXX JR XXX We need to think a better way of maintaining consistence
#       maybe collect all output directories and files in a single object
#       'out'
# NOTE: to be used to store/search annotation from now on
#   used to be 'folder'/'annotation', which now is a relative symlink to
ens.db <- NULL

if (USE.ONLINE.ANNOTATION == TRUE) {
    # get access to annotation hub
    ah <- AnnotationHub()
    # get ENSEMBL annotations for our target
    qr <- query(ah, c("EnsDb", target.organism))
    if (length(qr) > 0) {
        # there is at least one, choose the most recent (last) one
        #	NOTE: this might not be desired sometimes
        last.ref <- last(names(qr))
        ens.db <- qr[[last.ref]]
    }
} else {
    # build annotation from GTF file
    DBFile <- build.offline.annotation( 
                              annotation.dir,
                              db.dir=paste('EnsDb', release, ens.version, sep = "_"),		##ADRIAN
                              target.organism=target.organism,
                              reference.gtf=reference.gtf,
                              release=release,
                              ens.version=ens.version,
                              user=my.user,
                              pass=my.password,
                              author=my.name,
                              maintainer=paste(my.name, my.email), 
                              license="Artistic-2.0")

    ens.db <- EnsDb(DBFile)
}

if (VERBOSE) {
    # check it
    cat.info("EnsDB:")
    print(ens.db)
    columns(ens.db)
}


###################################################
###   A L T E R N A T E   A N N O T A T I O N   ###
###################################################

# ---------------------------------------------------------------
# M A K E     O R G     P A C K A G E
# ---------------------------------------------------------------

cat('\n\n')
short.title('Org DB Annotation')

org.db <- NULL

if (USE.ONLINE.ANNOTATION == TRUE) {
    # use AnnotationHub to seek a suitable package
    #ah <- AnnotationHub()	# we have already open the connection
    # check if an Org.db is available for our organism)
    qo <- query(ah, c("Orgdb", org.package))
    last.ref <- last(names(qo))
    org.db <- qo[[last.ref]]
    
	#if ( last.ref %in% names(qo) ) {
    #    org.db <- query(ah, "Orgdb")[[ last.ref ]]
    #} 
	
	# If AnnotationHub does not work, try to retrieve NCBI annotation
	# directly from org.package if it exists.

	if ( is.null(org.db)) {			
		if ( ! require(org.package, character.only = TRUE, quietly = TRUE)) {
			BiocManager::install(org.package) }

		# Character.only allows to require package form variable name
		require(org.package, character.only = TRUE)
		org.db <- get(org.package)
	}
	
} else {
    # Build offline annotation if none is available locally
    if ( ! require(org.package, character.only=T)) {
        #-----------------------------------------------
        # Create Org Package
        #-----------------------------------------------
        #
        # The datasets were first downloaded by hand from NCBI
        # and SwissProt into org.Xx.eg.db.
        # Then the following command had to be used:
        #
        makeOrgPackageFromNCBI(			## ADRIAN
                author  = "J. R. Valverde <jrvalverde@cnb.csic.es>", 
                maintainer = "J. R. Valverde <jrvalverde@cnb.csic.es>", 
                tax_id  = ncbi.taxid, # from NCBI Taxonomy browser
                genus   = split(target.organism, ' ')[1], 
                species = split(target.organism, ' ')[2], 
                version = "0.1",        # or better, the genome's version 
                outputDir = "./org", 
                NCBIFilesDir = "./ncbi"#, 
                #rebuildCache=FALSE
                )
        # We specify a directory to save locally the files used (and 
        # retrieved from NCBI) to create the Org database.
        # Normally, if the files are older than one day, they will be
        # downloaded again. Since they are very large, the download may
        # take too long and be interrupted frequently. This implies that
        # if we are unable to download everything in a single day, we are
        # doomed. The 'rebuildCache' option should avoid the re-downloading,
        # but it is described as an option for "internal use only and for
        # testing", and it may imply that no files are downloaded at all
        # and only local files in the cache directory are used, so we will
        # not use it unless strictly necessary.
        # GET
        # ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/*
        # ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/ARCHIVE/gene2unigene
        # ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
        # if that fails try using https://... instead  <<< PREFERRED!!!
        # if that fails try curl

        # We can download the files by hand first and then run this command
        # so it uses the already downloaded files.
        
        # install the package
        install.packages(org.package, repos=NULL)
        if (require(org.package, character.only=T)) {
            org.db <- eval(parse(text=org.package))
            columns(org.db)
        } else {
            org.db <- NULL
            cat.warn('\nCould not generate org.db\n')
        }
        ### NOTE
        # This fails in Gallus gallus due to download failures from NCBI
        # because of excesses in download times.
        # We need to use the latest GitHub version installed with
        # library(devtools)
        # install_github("Bioconductor/AnnotationHub")
    }
}

# At this point ens.db should contain the ENSEMBL data and 
# org.db the Org data.


# ---------------------------------------------------------------
# O B T A I N   B I O M A R T   A N N O T A T I O N
# ---------------------------------------------------------------
cat('\n\n')
short.title('BiomaRt Annotation')

bm.1.file <- file.path(annotation.dir, 'biomaRt.annotation.1st.txt')
print(bm.1.file)
print(annotation.dir)

if ( file.exists(bm.1.file) ) {
    cat.nl("Loading existing biomaRt annotation")
   # we have aready retrieved and saved the annotation, use it
    bm.annot.1 <- read.table(
    		file=file.path(annotation.dir, 'biomaRt.annotation.1st.txt'), 
	        sep='\t', header=T)
} else {
    cat.nl("Building biomaRt annotation")
    bm.annot.1 <- get.biomart.annotation(annotation.dir)
}


# ---------------------------------------------------------------
# O B T A I N   S T R I N G   P P I   A N N O T A T I O N
# ---------------------------------------------------------------
cat('\n\n')
short.title('Getting STRINGdb PPI')
ppi <- get.stringsdb.ppi(ncbi.taxid=ncbi.taxid, 
                         mart.name=mart.name, 
                         save.dir=annotation.dir)

# at this point we have all the different annotation subsets we would
# like to use.

# And now we are ready with our reference info at hand...


##########################################
#										 #
# Get the alignments and feature counts	 #
#										 #
##########################################
#
# The next step is to align the reads and calculate the counts of reads
# that map to each gene after obtaining alignments in BAM format: 
# the .bam files (binary files containing the reads aligned to the reference 
# genome sequence) and the .bam.bai files containing the indexes for each 
# .bam one. These two steps take the most of the time.
# 
# We need to calculate the counts of the reads that map to gene regions,
# summing up the reads that match each gene. For this, we use each .bam file
# and process it with the function featureCounts from R package Rsubread, and
# the reference genome data, which is in a GTF file. This last file
# contains the information about the features annotated in the genome,
# including the coordinates of each feature in the reference genome. 
# The function featureCounts will use these coordinates in the genome and
# the coordinates where each read maps to the genome in the BAM file to know 
# to which feature each read maps.
#
if ( ALIGN == TRUE ) {
    cat('\n\n')
	short.title("Aligning Reads")
	dir.create(alignment.dir, showWarnings = FALSE)
    bam.files <- align.fastq(fastq.dir = fastq.dir,
							reference = reference,
							alignment.dir = alignment.dir,
                            paired=PAIRED) 
}

# get the TABLE with the feature counts
cat('\n\n')
short.title("Getting feature counts")
fc <- get.feature.counts(file.path(folder, feature.count.dir), alignment.dir, PAIRED, BOTH)

#########################################################################################

# We are ready; we have
#	annotation (ensembldb, org.db, biomart)
#	feature counts table

########################################################################################

# PREFILTERING STEP
# Discard genes with no counts in any of the samples
counts <- fc[rowSums(fc)>1, ]
head(counts)
#dim(counts)

## LOAD SAMPLE METADATA AND SET DESIGN COLUMN
sampleInfo <- read.table(metadata, header = TRUE)
target <- as.factor(sampleInfo[, design.column])
sampleInfo[, design.column] <- target

# ---------------------------------------------------------------
# A N A L Y S I S     W I T H     E D G E R
# ---------------------------------------------------------------
if (USE.EDGER) {

	cat("\n\n")
    short.title('EdgeR Analysis')
	Sys.sleep(1)
    
	dge <- eR.differential.gene.expression(counts, 
    			metadata = sampleInfo,
				design.column = design.column,
                cpm.threshold = cpm.threshold,
				n.genes = n.genes, 
                ens.db = ens.db, 
                folder = folder)
    
    fit <- eR.voom.variation.analysis(dge, folder)

    fit <- eR.fit.annotate.save(fit,
    					   		ens.db, 
                                biomart.ann, 
                                folder, 
                                n.genes)

    # Now create Volcano plot for each coefficient (comparison) from our fitted
	# model, for only the top 1/2 genes using gene annotation
	for (i in 1:ncol(fit)){
		coefficient <- colnames(fit)[i]
		out.png <- file.path(folder, 'edgeR', 'img', 
                   paste('edgeR_volcanoplot_', coefficient,'.png', sep=''))
		as.png(volcanoplot(fit, highlight = n.genes/2, coef = i, names=fit$genes$SYMBOL, cex=0.8, pch=1),
        	out.png)
	}
    
    # Testing Differentially expressed genes relative to a threshold
    fit.thres <- eR.fit.treat(fit, threshold = significance.threshold, folder)
    
    # Save all the contents of 'fit' in an RDS file
    saveRDS(fit, file=file.path(folder, 'edgeR', 'annotatedVOOMfit+.rds'))
    saveRDS(fit.thres, 
            file=file.path(folder, 'edgeR', 
                       paste('annotatedVOOMfit+.gt.',
                       significance.threshold, '.rds', sep='')))

    # and save as well as Rdata file
    save(fit, file=file.path(folder, 'edgeRannotatedVOOMfit+.RData'))
    save(fit.thres, 
         file = file.path(folder, 'edgeR', 
              paste('annotatedVOOMfit+.gt.',
                    significance.threshold, '.RData', sep='')))


    # -----------------------------------------------------------------
    # Getting beyond here is easier with DESeq2
    # -----------------------------------------------------------------
    # 
    # # redefine the design so we can make any kind of comparison by
    # # using a null reference ("0 + ")
    # 
    # this carries out all comparisons and annotates each fit obtained
    # saving all of them in a list
    eR.data <- eR.dge.all.comparisons(	dge = dge,
										design.column,
										ens.db = ens.db,
										biomart.db = bm.annot.1,
										folder = folder)

    # now eR.data is a list where each element is a comparison A-B (A minus B),
    # i.e. each A-B is a list of tables, one of them named "table" and containing
    # logFC, logCPM, F and PValue
    #
    # We would like to have FDR-corrected p-values as well, which can be got by
    # defining a threshold. But that implies we know of a meaningful one, which
    # we don't yet.

	short.title('EdgeR Analysis Finished !')
	Sys.sleep(1)

}    # end if (USE.EDGER)



#################################################################
#
#################################################################
                                                               
#################################################################

# ---------------------------------------------------------------
# A N A L Y S I S     W I T H     D E S E Q 2
# ---------------------------------------------------------------

if (USE.DESEQ2) {

    cat("\n\n")
	short.title('DESeq2 Analysis')
	Sys.sleep(1)

	countData <- counts     # use prefiltered (non-zero) counts
	colData <- sampleInfo
    if (VERBOSE) {
        cat.info("Obtained gene counts (first two rows)")
        print(head(countData, 2))
	}

    # save normalized counts in addition to raw feature counts
    # so we can use them for AI analysis
    # SAVE NORMALIZED COUNTS FOR AI ANALYSIS
    #
    # We do not need normalized counts for differential gene expression
    # because DESeq2 calculates them automatically.
    #
    # But we do need the normalized counts for the AI analysis, so we 
    # will get and save them.
    #
    # edgeR-type TMM normalization depends on the comparison design
    # try both and see if there are differences (tried, no relevant diffs)
    # sort metadata
    #
    # Note: DESeq2 doesn't actually use normalized counts, rather it uses
    # the raw counts and models the normalization inside the Generalized Linear
    # Model (GLM). These normalized counts will be useful for downstream
    # visualization of results, but cannot be used as input to DESeq2 or any other
    # tools that peform differential expression analysis which use the negative
    # binomial model.
    #
    # To use these counts directly for variable reduction, we should first 
    # select significant genes to reduce the number of explanatory variables
    # and we should make sure that in the design, the reference value (e.g.
    # "wt" for comparisons) is the first element in countData.
    #
    # IMPORTANT NOTE: !!!
    # We need to have metadata and feature counts in the appropriate
    # name order for DESeq2 to work correctly.
    #
    norm_counts_file <- paste(folder, "/feature_counts/normalizedCounts_by_", column, ".csv", sep='')
    dds.col <- DESeqDataSetFromMatrix(countData=countData, 
                colData=colData,  
    dds.col <- estimateSizeFactors(dds.col)	# add norm factors to dds.col
    sizeFactors(dds.col)
    normalized_counts <- counts(dds.col, normalized=TRUE)	# get normalized counts
                design=as.formula(paste("~", design.column)))
    write.table(normalized_counts, 
                file=norm_counts_file,
	            sep='\t', row.names=T, col.names=T) # better if row.names=F
    
    
	# Create DESeq objects to compare samples by design.column
    out.base <- file.path(folder, "DESeq2",
        paste("dds.DESeq2.", design.column, sep=''))
    
	if ( file.exists(paste(out.base, "rds", sep='.')) ) {
        dds <- readRDS(file=paste(out.base, "rds", sep='.'))	
		
	} else {
	
		cat("\n\tGenerating DESeq Object\n\n")
        dds <- DESeqDataSetFromMatrix(countData = countData, 
                                      colData = colData,  
                                      design = as.formula(paste("~", design.column)), 
                                      tidy=F)
        dds <- DESeq(dds)		
		
        # SAVE DESEQ ANALYSIS
        # -------------------        
		save(dds, file = paste(out.base, "RData", sep='.'))
        saveRDS(dds, file = paste(out.base, "rds", sep='.'))
    }
	
	# Do some diagnostic plots - PCA
	vst <- varianceStabilizingTransformation(dds)
	out.png <- file.path(folder, "DESeq2", "img", 
               paste("dds.DESeq2.", design.column, "_PCA_plot.png", sep=''))
	as.png(plotPCA(vst, intgroup = design.column),
        	out.png)
	
	# Plot Dispersion Estimates
	out.png <- file.path(folder, "DESeq2", "img",
               paste("dds.DESeq2.", design.column, "_DispEstimates.png", sep=''))
	as.png(plotDispEsts(dds, main = "DESeq2 Per-gene Dispersion Estimates"),
        	out.png)
	
	
    cat("\n\tDESeq2 analysis produced the following models\n\n")
    print(resultsNames(dds))
	
	
    ### XXX JR XXX ### THIS SHOULD BE A SINGLE FUNCTION LIKE WITH EDGER!!!

	## Perform all pairwise comparisons available in DESeq2 model
	## DDS objects only include comparisons that take the first level
	## of the design column as reference.
	## If we wan to do all possible comparisons, we need to relevel.
    ##  relevel re-orders the levels of a factor so that the level
    ##  specified by "ref" becomes the first and all others are moved 
    ##  down
	
	dds.results <- list()
	levels <- levels(dds[[design.column]])	
	
	for (ref in levels) {

		# Reorder levels of dds$sample (not dds$sample)
		if (length(levels(dds$sample)) > 2) {
            dds$sample <- relevel(dds$sample, ref)
		}
        # and repeat the DE analysis using 'ref' as reference level
        dds <- DESeq(dds)
		
		# All comparisons for that reference
		comparisons <- resultsNames(dds)[2:length(resultsNames(dds))]

		for (comparison in comparisons){
			if (comparison == "Intercept") next

        	cat("\n                                    ")						
			cat("\nComputing DGE:",	  comparison,	"\n")
			cat(  "==================================\n")

            contrast <- str_split(comparison, pattern = "_")[[1]][c(1, 2, 4)]	
            # Eg: c("Sample", "DF1.PC", "DF1") or c(col, val1, val2)
            # This fails if values have '_' in them. 
            # Assuming that the column name doesn't contain any '_', 
            # (otherwise there is in principle no safe way to separate them)
            # we could use instead
            c  <- str_split(comparison, pattern="_")[[1]][1]
            v1 <- gsub("^.*?_","", str_split(comparison, pattern="_vs_")[[1]][1])
            # the '?' makes matching lazy so it only goes up to the first '_'
            # or we could use
            #v1 <- gsub(paste("^", c, "_" sep=''), '', str_split(comparison, pattern="_vs_")[[1]][1])
            v2 <- str_split(comparison, pattern="_vs_")[[1]][2]
            contrast <- c(c, v1, v2) 
            if (VERBOSE) cat.info(">>> comparison =", comparison)
            if (VERBOSE) cat.info(">>> contrast =", contrast)
			out.file.base <- paste("dds", comparison, sep = ".")				#Ex: "dds.Sample.DF1_PC.vs.DF1"

			# Perform comparison, shrinkage, annotate and identify significant changes.
			dds.results[[comparison]] <- DESeqCompareAndPlot(
								dds = dds,
								contrast = contrast,
								filterFun = ihw,
								alpha = 0.01,
								shrnk.type = 'apeglm',
								annotate = TRUE, 
								ensembl.db = ens.db,
								biomart.db = bm.annot.1,
								org.db = org.db,
								annotation.dir = file.path(folder, "annotation"),
								save = TRUE,
								out.file.base = out.file.base,
								out.dir = folder)

			if (VERBOSE) dds.results[[comparison]]$result.ann %>% head(10)
		}
	}

    #annot <- ds2.get.annotation.ens.org(dds, ens.db, org.db)
    #ens.ann <- annot$ens.ann		### NOTE we could use attach here
    #ens.ann.1 <- annot$ens.ann.1	# but this is more explicit
    #geneSymbols <- annot$geneSymbols
    #go.ann <- annot$go.ann
    #GOdescription <- annot$GOdescription

    # Add annotation to dds to keep everything in one place
    #dds$ens.annot.1 <- ens.ann.1
    #dds$bm.annot.1 <- bm.annot.1


    # Now result list contains annotation of EnsemblDB, BiomaRt and ORG.Db (if provided)
	# Org.DB packages are not available for all organisms, so by default, is set to NULL.

	if (FALSE){
	if ( ! is.null(org.db) ) {

		for (res in names(dds.results)){

			res.ann <- dds.results[[res]]$result.ann

			# And now we should be able to use the Org package if we successfully built it
			# at the beginning.		
			gene.id <- rownames(res.ann)
			entrez <- as.character(res.ann$ENTREZID)

			# Check the information we can retrieve from the database
			keytypes(org.db)

			# Retrieve Gene Ontology IDs
			ncbi.ann <- AnnotationDbi::select(org.db, 
				    		keys = gene.id, 
				    		columns = c("ENTREZID", "GENENAME", "GENETYPE", "SYMBOL",
										#"ALIAS", "ENZYME", "ACCNUM", "PMID",
										#"UNIPROT", "PROSITE", "PFAM",
										"GO", "GOALL","ONTOLOGYALL", "PATH"), 
				    		keytype = "SYMBOL", 
				    		multiVals = "CharacterList")

			# Retrieve Gene Ontology descriptions
			GO.ann <- AnnotationDbi::select(GO.db,
							keys = ncbi.ann$GO, 
							columns = c(	"GOID", "TERM", "DEFINITION", "ONTOLOGY"),
							keytype = "GOID")

			for (i in colnames(GO.id)){
				go.ann[i] <- Go.id[match(res.ann$gene.id, ensembl.ann[,by]), i]
			}	
		}
	}
	}
    # SAVE SELECTED COMPARISONS
    # -------------------------
    # We will get all comaprisons
    #for (a in levels(as.factor(target[ , design.column])) ) {
    #    for (b in levels(as.factor(target[ , design.column])) ) {
    #        print(paste(a, b))
    #        if (a == b) next
    #        ds2.dds.compare.annotate.save( 
    #                        dds,
    #                        column=design.column,
    #                        x=a, y=b,
    #                        filterFun=ihw, alpha=0.01,
    #                        ensembl.ann=ens.ann.1,
    #                        biomart.ann=bm.annot.1,
    #                        outDir=folder
    #                        )
	#
    #    }
    #}
	#
	#
	#
    # SAVE GENES WITH SIGNIFICANT CHANGES
    # -----------------------------------
	#
    #ds.data <- list()
    #x <- 0
    #contr <- design.column
    #grps <- levels(colData[ , contr ])
    #for (a in grps) {
    #    for (b in grps) {
    #        print(paste(a, b))
    #        if (a == b) next
    #        res <- ds2.dds.plot.and.save( 
    #                        dds, contr,
    #                       x=a, y=b,
    #                        filterFun=ihw, alpha=0.01,
    #                        ensembl.ann=ens.ann.1,
    #                        biomart.ann=bm.annot.1,
    #                        outDir=folder
    #                        )
    #        print(names(res))
	#        name <- paste(contr, "_", a, "_", b, sep='')
    #        #x <- x + 1
    #        #ds.data[[x]] <- res
    #        ds.data[[name]] <- res
    #        #names(ds.data)[x] <- name
    #        #stop()
    #    }
    #}

    if (INTERACTIVE) {
        ds2.interactively.print.n.most.significant.de.genes(dds.results, n=10)

        ds2.interactively.plot.top.up.down.regulated.gene( dds, 
                                                           dds.results, 
                                                           design.column)
    }


    # -------------------------------
    # Do Gene Set Enrichment Analysis
    # -------------------------------

	short.title('Gene Set Enrichment Analysis')

	ds.data <- dds.results

	# For GSEA, we need to define a maximum number of genes we want to include per group
	# when performing clustering. For large groups, the ontology will be hierarchically 
	# superior and the functions retrieved will be generic. In contrastm, for small groups,
	# the ontology will be more specific.
	# We do not know what is the best upper limit, so we'll try several

	########################
	#   CLUSTER PROFILER   #
	########################

    # minGSSize (minimal size of each geneSet for analyzing) is hardcoded (3-10)
    #   in the subroutines below, default is 10
    # maxGSSize (maximal size of genes annotated for testing) is switched here
    #   default is 500
    # these values will affect the results by detecting more or less bigger or
    # smaller groups: too small or too large groups have more chances of
    # being statistically significative; this may show as differences when
    # showing (e.g. emapplot) the highest p-valued groups only
    for (ms in c(100, 500)) {
    #for (ms in c(50, 100, 250, 500)) {
	#for (ms in c(5)) {
		for (cmp.name in names(ds.data)) {
			
			#if (str_detect(cmp.name, "_vs_DF1$")) next
            
			# for each comparison name
            cmp.data <- ds.data[[ cmp.name ]]
        	cat('\n\nDoing Enrichment with clusterProfiler on', cmp.name, '\n') 
            ann.shrunk.lfc <- cmp.data[[ "shrunk.ann" ]] %>% as.data.frame()

            out.dir <- sprintf("%s/DESeq2/GO+KEGG_cProf/max.size=%03d/%s",
                    folder, ms, cmp.name)
            dir.create(out.dir, showWarning = FALSE, recursive = TRUE)

			GO_KEGG_clusterProfiler(ann.data = ann.shrunk.lfc,
                                    ranking.column = 'log2FoldChange',
									max.size = ms,
									out.dir = out.dir, 
									out.name = 'cProf',
						  			use.description = TRUE,
						  			OrgDb = org.db,
									kegg_organism = kegg.organism,
						  			top.n = 10,
						  			top.biblio = 5,
						  			verbose = VERBOSE)
		}
    }


	#####################
	#   FGSEA PACKAGE   #
	#####################

	for (ms in c(100, 500)) {
#	for (ms in c(50, 100, 250, 500)) {
    	for (cmp.name in names(ds.data)) {
        	# for each comparison name
        	cmp.data <- ds.data[[ cmp.name ]]
        	cat('\n\nDoing Enrichment with FGSEA on', cmp.name, '( max size =', ms, ')\n') 
        	ann.shrunk.lfc <- cmp.data[[ "shrunk.ann" ]] %>% as.data.frame()
        	#out.dir <- file.path(folder, "DESeq2", "GO_fgsea", cmp.name)
            # XXX JR XXX CHANGE TO USE file.path()
        	out.dir <- sprintf("%s/DESeq2/GO_fgsea/max.size=%03d/%s", folder, ms, cmp.name)
        	dir.create(out.dir, showWarning = FALSE, recursive = TRUE)
        	gogsea <- GO_fgsea(ann.data = ann.shrunk.lfc, 
                       ranking.column = 'log2FoldChange',
                	   max.size = ms,
                	   out.dir = out.dir, 
                	   out.name = 'fgsea',
                	   use.description = TRUE,
                       top.n = 10,
                       top.biblio = 5,
                	   verbose = VERBOSE)
     	}
	}


    # -------------------------
    # Analyze GO representation
    # -------------------------
    ds2.analyze.go.representation(ds.data = ds.data, bm.go.annot = bm.annot.1, folder)

    ds2.analyze.go.ancestry(ds.data = ds.data, l2fc.threshold = 0, folder)

    ds2.analyze.pfam.representation(ds.data = ds.data, bm.fam.annot = bm.annot.1, folder)


    # ------------------------------
    # DO PCA AND CLUSTERING ANALYSES
    # ------------------------------

    ## CREATE THE COLUMNS FOR THE PCA ANALYSIS 
    # we need to retrieve from the dataframes only the log2Fold change and the
    # gene names
    # we will take the data from the unsorted tables signif_* 
    # we save the rownames as a distinct column and we pass it as a column 
    # we delete column1 (gene-names) before clustering because it is not needed

    # now that we have the data we need to keep only the common rows to all of them
    # we use the function intersect by pairs and then all together
    # the common genes are rows from common2
    # finally we create a dataframe with everything we have


    # This is here in case we decide to loop over several columns later 
    contrasts.column <- design.column
    references <- levels(as.factor(colData[ , contrasts.column ]))

    for (ref in references) {
        # CLUSTER GENES

        # 1. Find all the genes that are common to all samples
		
		# Create a convenience text variable to simplify/unify filenames below
        ccol_ref <- paste(contrasts.column, ref, sep='_')

        # First define all genes detected in results
        common <- rownames(ds.data[[1]]$result.ann)
        
		# Now find common genes for all comparisons 
        name=''
		for (i in references) {
            if (i == ref) next
            name <- paste(ccol_ref, "vs", i, sep='_')
            cat(name)
            common <- intersect(common, rownames(ds.data[[name]]$signif))
			#print(length(common))
        }
        cat(paste("\n', name, '\t", length(common), "common genes detected\n\n"))
        if (length(common) < 2) {
            cat.warn("=======================================================")
            cat.warn(">>>>>>>>>>>>>>>>        ", name)
            cat.warn("NOT ENOUGH SIGNIFICANTLY DIFFERENT GENES FOR CLUSTERING")
            cat.warn("=======================================================")
            next
        }

		# Create data table with common genes
        data.table <- data.frame(genes = common)
        rownames(data.table) <- common
        
		# 2. Extract LFC values for each comparison and
		# add it to the common genes table
        
		for (i in references) {
            if (i == ref) next
            name <- paste(ccol_ref, "vs", i, sep='_')
            data.table[name] <- ds.data[[name]]$signif[common, "log2FoldChange"]
        }
        
        # 3. Add annotation to the data
        
        data.annot <- ds.data[[name]]$signif.annot[common, ]
        
        # 4. Prepare for clustering
        par(mfrow=c(1,1))
        set.seed(1963)

        # 5. Have a general look at the methods to get a feeling for the
        # best number of clusters

        dif <- data.table[ , -1]
        # if there is only one reference, then dif only has one column
        #   this may happen if we only have two levels (e.g. Y/N)
        #if (dim(dif)[1] == 0 || dim(dif)[2] == 0) {
        # This check should no longer be needed since we check length(common)
        if (isEmpty(dif)) {
            cat('\n')
            cat.warn("=========================================================")
            cat.warn(">>>>>>>>>>>>>>>>        ", name)
            cat.warn("INSUFFICIENT SIGNIFICANTLY DIFFERENT GENES FOR CLUSTERING")
            cat.warn("=========================================================")
            next
        }
        by.row <- 1
        by.col <- 2

		# 5.1 Normalize LFC by calculating Z-score (x-mean)
        if (is.vector(dif)) {
            #means <- mean(dif)
            #sds <- sd(dif)
            #nor <- scale(dif,center = means,scale = sds)
            dif <- as.matrix(dif)
        } 
        if (dim(dif)[2] == 1) {     # n.common.genes
            cat.warn("========================================================")
            cat.warn("ONLY ONE SIGNIFICANTLY DIFFERENT GENE FOR CLUSTERING")
            cat.warn("========================================================")
            next
        }
        means <- apply(dif, by.col, mean)
        sds <- apply(dif, by.col, sd)
        nor <- scale(dif,center = means,scale = sds)
		
        if (INTERACTIVE) {
            # Do a scatterplot matrix
            car::scatterplotMatrix(dif)
        }
        out.png <- file.path(folder, "DESeq2", "img", 
                   paste('scatterplot_matrix_', ccol_ref, ".png", 
                   sep=''))
        as.png( {
                print(car::scatterplotMatrix(dif))
                }, out.png)


        ## Try to guess the optimum number of K-means clusters with NBClust
        #out.log <- paste(folder, "/DESeq2/cluster/NBClust_", ccol_ref, ".log", sep='')
        #openlog(out.log)
        #
        ## NbClust helps us predict best number of clusters using
		## Hubert index and D index
        #tryCatch(
        #    nbc <- NbClust(dif, diss=NULL, 
        #            distance="euclidean", method="complete", 
        #            min.nc=3, max.nc=10, 
        #            index="all", alphaBeale=0.1),
        #    error=function(err) { cat.nl("NBClust failed (",err$message,")") }   
        #)
        #cat("NbClust recommends", nbc, "clusters\n")
        #sink()  # close last log file

        # 5.2. cluster genes by k-means
        # Let the user see the various clusters and make a decision
        #
        out.log <- file.path(folder, "DESeq2", "cluster", "kmeans",
                         paste("DESeq2_kmeans_all_", ccol_ref, ".log", sep=''))
        openlog(out.log)
        for (i in 2:10) {
            cat("\n\tClustering with K-means (", i, " clusters)\n\n", sep = "")
            cl <- kmeans(nor, i)				# NOTE: nor
            cat("Cluster membership summary:\n", table(cl$cluster), "\n")
            tryCatch(
                print(fviz_cluster(cl$cluster, geom = "point", data=nor)),
                error=function(e) {
                        cat.err("Cannot plot cluster", i, '\n', error=e)
                        #str(e)
                      }
            )
            out.png <- sprintf(
                    "%s/DESeq2/cluster/kmeans/DESeq2_kmeans_%s_nc=%03d.png", 
        	    folder, ccol_ref, i)
            as.png(fviz_cluster(cl, geom = "point", data=nor,
                        main=paste("K-means clsutering nc =", i) ),
                    out.png)
            if (INTERACTIVE)
                continue.on.enter("Press [ENTER] to continue ")
        }
        sink()

        cat("The log and plots for K-means clustering have been saved to\n", 
            "\t", folder, "/DESeq2/cluster/kmeans/\n", 
            "please, inspect them and select the best number of clusters\n",
            sep='')

        #kmeans.nc <- readline("How many clusters should I use for k-means? ")
        #kmeans.nc <- as.numeric(kmeans.nc)
        # for gg
        kmeans.nc <- 5

        # 5.3 cluster genes by PAM
        out.log <- file.path(folder, "DESeq2", "cluster", "pam",
                   paste("DESeq2_pam_all_", ccol_ref, ".log", sep=''))
        openlog(out.log)
        for (i in 2:10) {
            cat("Clustering with PAM (", i, " clusters)\n")
            cl <- pam(nor, i)					# NOTE:nor
            print(table(cl$cluster))
            print(fviz_cluster(cl, geom = "point"))
            out.png <- sprintf("%s/DESeq2/cluster/pam/DESeq2_pam_%s_nc=%03d.png", 
                folder, ccol_ref, i)
            as.png(fviz_cluster(cl, geom = "point",
                        main=paste("Partitioning Around Medoids clsutering nc =", i) ), 
                    out.png)
            #continue.on.enter("Press [ENTER] to continue ")
        }
        sink()

        cat("The log and plots for PAM clustering have been saved to\n", 
            "\t", folder, "/DESeq2/cluster/pam/\n", 
            "please, inspect them and select the best number of clusters\n",
            sep='')

        #pam.nc <- readline("How many clusters should I use for PAM? ")
        #pam.nc <- as.numeric(pam.nc)
        # for gg
        pam.nc <- 7

        # 5.4 cluster genes with DBscan
        out.log <- file.path(folder, "DESeq2", "cluster", "dbscan",
                   paste("DESeq2_dbscan_all_", ccol_ref, ".log", sep=''))
        openlog(out.log)
        for (e in c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.0)) {
            cat("\nClustering with DBScan ( eps =", e, " )\n")
            cl <- dbscan(dif, eps=e, MinPts=10, showplot=1)	# NOTE: dif
            print(table(cl$cluster))
            print(fviz_cluster(cl, data=dif, 
                         show.clust.cent=FALSE, labelsize=4,
                         ellipse=TRUE, ellipse.type="convex"))
            out.png <- sprintf("%s/DESeq2/cluster/dbscan/DESeq2_dbscan_%s_eps=%03.2f.png", 
                         folder, ccol_ref, e)
            as.png(fviz_cluster(cl, data=dif, 
                          show.clust.cent=FALSE, 
                          geom="point", #labelsize=4,
                          ellipse=TRUE, ellipse.type="convex",
                          main=paste("Density based clustering eps =", e) ),
                    out.png)
            #continue.on.enter("Press [ENTER] to continue ")
        }
        sink()

        cat("The log and plots for DBscan clustering have been saved to\n", 
            "\t", folder, "/DESeq2/cluster/dbscan/\n", 
            "please, inspect them and select the best epsilon value\n",
            sep='')

        #dbscan.eps <- readline("Which epsilon should I use for DBscan? ")
        #dbscan.eps <- as.numeric(dbscan.eps)
        # for gg
        dbscan.eps <- 1.5

        # 5.5 cluster genes with hierarchical clustering
        hclust.nc <- 5	# we'll set it by hand for now

        # 6 Repeat in detail with each method
        # Now, proceed in detail, method by method, with a deeper and more 
        # detailed analysis

        # 6.1 hierarchical clustering (hclust)
        # possible distances are c("euclidean", "maximum", "manhattan", "canberra",
        #	"binary", "minkowski", "pearson", "spearman", "kendall")
        # possible methods are c("ward.D", "ward.D2", "single", "complete",
        #	"average" (UPGMA), "mcquitty" (WPGMA), "median" (WPGMC), "centroid" (UPGMC))
        out.log <- file.path(folder, "DESeq2", "cluster", "hcluster", 
                   paste("DESeq2_hcluster_", contrasts.column, "_", ref, ".txt", sep=''))
        sink(out.log, split=T)
        h.cluster.changes(
            #data.table=pca.table[ , -1],
            #annotated.data=ds.data[[name]]$signif.annot[common, ],
            clusters=hclust.nc,
            data.table=dif,
            annotated.data=data.annot
        )
        sink()

        hclust_folder <- file.path(folder, "ESeq2", "cluster",
                         paste("hclust_", ccol_ref, sep=''))
        dir.create(hclust_folder, showWarnings=FALSE)

        # 6.2 hierarchical clustering (hcut)
        hcut_folder <- file.path(folder, "DESeq2", "cluster", 
                       paste("hcut_", ccol_ref, sep=''))
        dir.create(hcut_folder, showWarnings=FALSE)

        out.log <- file.path(hcut_folder, 
                   paste("/DESeq2_hcut_", ccol_ref, ".txt", sep=''))
        sink(out.log, split=T)
        
        boot <- 100

        hc.cl <- cluster.changes(
            #data.table=pca.table[ , -1],
            #annotated.data=sa.data[[name]]$signif.annot[common, ],
            data.table=dif,
            annotated.data=data.annot,
            FUN=hcut,
            clusters=hclust.nc,
            estimate=TRUE,
            gap_bootstrap=boot,
            output.folder=hcut_folder
            )
        sink()

        # 6.3 k-means
        kmeans_folder <- file.path(folder, "DESeq2", "cluster",
                 paste("kmeans_", ccol_ref, sep=''))
        dir.create(kmeans_folder, showWarnings=FALSE)

        out.log <- file.path(kmeans_folder, 
                paste("DESeq2_kmeans_", ccol_ref, ".txt", sep=''))
        sink(out.log, split=T)
        km.cl <- cluster.changes(
            #data.table=pca.table[ , -1],
            #annotated.data=ds.data[[name]]$signif.annot[common, ],
            data.table=dif,
            annotated.data=data.annot,
            FUN=kmeans,
            algorithm="Hartigan-Wong",
            clusters=kmeans.nc,
            nstart=kmeans.nc*10,
            estimate=T,
            gap_bootstrap=boot,
            output.folder=kmeans_folder
            )
        sink()

        # 6.4 PAM
        pam_folder <- file.path(folder, "DESeq2", "cluster",
            paste("pam_", ccol_ref, sep=''))
        dir.create(pam_folder, showWarnings=FALSE)

        out.log <- file.path(pam_folder, paste("DESeq2_pam_", 
                         ccol_ref, ".txt", sep=''))
        sink(out.log, split=T)
        pam.cl <- cluster.changes(
            #data.table=pca.table[ , -1],
            #annotated.data=ds.data[[name]]$signif.annot[common, ],
            data.table=dif,
            annotated.data=data.annot,
            FUN=pam,
            distance="euclidean",	# see dist()
            clusters=pam.nc,	# n. of clusters
            nstart=pam.nc*10,	# n. of random start sets to choose
            estimate=T,
            gap_bootstrap=boot,
            output.folder=pam_folder
            )
        sink()

        # 6.5 DBscan
        dbscan_folder <- file.path(folder, "DESeq2", "cluster",
               paste("dbscan_", ccol_ref, sep=''))
        dir.create(dbscan_folder, showWarnings=FALSE)

        out.log <- file.path(dbscan_folder, paste("DESeq2_dbscan_", 
                         ccol_ref, ".txt", sep=''))
        sink(out.log, split=T)
        dbs.cl <- cluster.changes(
            #data.table=pca.table[ , -1],
            #annotated.data=ds.data[[name]]$signif.annot[common, ],
            data.table=dif,
            annotated.data=data.annot,
            FUN=dbscan,
            distance="euclidean", # see dist()
            eps=dbscan.eps,
            normalize=F,
            estimate=T,
            gap_bootstrap=boot,
            output.folder=dbscan_folder
            )
        sink()
    }


    for (ref in references) {
        # Now cluster by experiment
        # -------------------------

        # 1. Find all the genes that are common to all samples
		# Create a convenience text variable to simplify/unify filenames below
        ccol_ref <- paste(contrasts.column, ref, sep='_')
        # First define all genes detected in results
        ## common <- rownames(ds.data[[name]]$signif)
        common <- rownames(ds.data[[1]]$result.ann)
        
		# Now find common genes for all comparisons 
        name=''
		for (i in references) {
            if (i == ref) next
            name <- paste(ccol_ref, "vs", i, sep='_')
            cat(name)
            common <- intersect(common, rownames(ds.data[[name]]$signif))
			#print(length(common))
        }
        cat(paste("\n', name, '\t", length(common), "common genes detected\n\n"))
        if (length(common) < 2) {
            cat.warn("=========================================================")
            cat.warn(">>>>>>>>>>>>>>>>        ", name)
            cat.warn("NOT ENOUGH SIGNIFICANTLY DIFFERENT GENES FOR CLUSTERING")
            cat.warn("=========================================================")
            next
            # XXX JR XXX Note that this implies that clustering
            # by sample will not be done either !!!
        }


		# Create data table with common genes
        data.table <- data.frame(genes = common)
        rownames(data.table) <- common
        
		# 2. Extract LFC values for each comparison and
		# add it to the common genes table
        for (i in references) {
            if (i == ref) next
            name <- paste(ccol_ref, "vs", i, sep='_')
            data.table[name] <- ds.data[[name]]$signif[common, "log2FoldChange"]
        }
        
        # 3. Add annotation to the data
        data.annot <- ds.data[[name]]$signif.annot[common, ]
        
        # 4. Prepare for clustering
        par(mfrow=c(1,1))
        set.seed(1963)
        by.row <- 1
        by.col <- 2
        
        dif <- data.table[ , -1]
        # we use similar but different messages to be able to track 
        # a problem
        if (isEmpty(dif)) {
            cat('\n')
            cat.warn("===================================================")
            cat.warn(">>>>>>>>>>>>>>>>        ", name)
            cat.warn("THERE IS NO DATA FOR EXPERIMENT CLUSTERING")
            cat.warn("===================================================")
            next
        }
        if (is.vector(dif)) dif <- as.matrix(dif)
        if (dim(dif)[2] == 1) {     # n.common.genes
            cat.warn("====================================================")
            cat.warn(">>>>>>>>>>>>>>>>        ", name)
            cat.warn("THERE IS NOT ENOUGH DATA FOR EXPERIMENT CLUSTERING")
            cat.warn("====================================================")
            next
        }
        # transpose the table so we work by experiment instead 
        fid <- t(dif)
        if (dim(fid)[2] == 1) {     # n.common.genes
            cat.warn("=============================================================")
            cat.warn(">>>>>>>>>>>>>>>>        ", name)
            cat.warn("THERE IS TOO FEW DATA FOR EXPERIMENT CLUSTERING")
            cat.warn("=============================================================")
            next
        }
        means <- apply(fid, by.col, mean)
        sds <- apply(fid, by.col, sd)
        ron <- scale(fid,center=means,scale=sds)

        # here we have a small number of rows and can set a maximum number of clusters
        maxclust <- nrow(fid) - 1
        
        # 5 CLUSTER BY EXPERIMENT
        # 5.1 hclust
        distan = dist(ron, method="euclidean")
        hcl <- hclust(distan)
        plot(hcl,labels=rownames(fid),main='Default from hclust')
        out.png <- sprintf(
                "%s/DESeq2_hcluster_grps_%s.png",
	        hclust_folder, ccol_ref)
        as.png(plot(hcl,labels=rownames(fid),main='Default from hclust'), out.png)

        # 5.2 kmeans
        fviz_nbclust(ron, kmeans, method="silhouette", k.max=maxclust)
        out.png <- sprintf(
                "%s/DESeq2_kmeans_grps_silhouette_%s.png",
	        kmeans_folder, ccol_ref, ccol_ref)
        as.png(fviz_nbclust(ron, kmeans, method="silhouette", k.max=maxclust), out.png)
        kcl <- kmeans(ron, centers=3, nstart=100)
        fviz_cluster(kcl, data=ron)
        out.png <- sprintf(
                "%s/DESeq2_kmeans_grps_%s_nc=%03d.png",
	        kmwans_folder, ccol_ref, ccol_ref, 3)
        as.png(fviz_cluster(kcl, data=ron), cout.png)

        # 5.3 pam
        fviz_nbclust(ron, pam, method="silhouette", k.max=maxclust)
        out.png <- sprintf(
                "%s/DESeq2_pam_grps_silhouette_%s.png",
	        pam_folder, ccol_ref)
        as.png(fviz_nbclust(ron, pam, method="silhouette", k.max=maxclust), out.png)
        pcl <- pam(ron, k=3, diss=F)
        fviz_cluster(pcl, data=ron)
        out.png <- sprintf(
                "%s/DESeq2_pam_grps_%s_nc=%03d.png",
	        folder, ccol_ref, 3)
        as.png(fviz_cluster(pcl, data=ron), out.png)

        # 5.4 dbscan
        # this results in the same two PCs but one cluster
        fviz_nbclust(ron, dbscan, method="silhouette", k.max=maxclust)
        out.png <- sprintf(
                "%s/DESeq2_dbscan_grps_silhouette_%s_%s.png",
	        dbscan_folder, ccol_ref)
        as.png(fviz_nbclust(ron, dbscan, method="silhouette", k.max=maxclust), out.png)
        
        ### JR ### eps should be tuned for each experiment
        ### we should likely use for loop
        dcl <- dbscan(ron, eps=0.2, MinPts=2, showplot=1)
        fviz_cluster(dcl, data=ron)
        out.png <- sprintf(
                "%s/DESeq2_dbscan_grps_%s_eps=%03.2f.png",
	        dbscan_folder, ccol_ref, 0.2)
        as.png(fviz_cluster(dcl, data=ron), out.png)

        # 5.5 this fails
        tryCatch( {
            NbClust(ron, diss=NULL, 
                    distance="euclidean", method="complete", 
                    min.nc=3, max.nc=10, 
                    index="all", alphaBeale=0.1)
        } )
    }

}



# ---------------------------------------- #
# GENERATE HTML REPORT FOR clusterProfiler #
# ---------------------------------------- #

## Run the Bash script to generate HTML index files
short.title("Generating report")

report.scr.dir <- file.path(my.dir, 'lib', "report")
wd <- getwd()

cat("\n\tGenerating HTML Report for clusterProfiler\n\n")

cProf.dir <- file.path(folder, "DESeq2", "GO+KEGG_cProf")

for (size.dir in list.dirs(cProf.dir, recursive = F, full.names = T)) {
	
    grp.size <- basename(size.dir)
	cat("\n\tGene Set", grp.size, "\n")
	cat(  "\t=====================")
    for (sample.dir in list.dirs(size.dir, recursive = F, full.names = T)) {
    
    	sample <- basename(sample.dir)
        cat("\n",sample, "\n\n", sep = "")

        cur.dir <- setwd(sample.dir)
        if (VERBOSE) cat.info(">>> working in", cur.dir)
		system(paste("bash", file.path(report.scr.dir, "/make-clusterprof-index.sh")))
		file.copy(file.path(report.scr.dir, "style.css"), ".")
		cur.dir <- setwd(wd)
        if (VERBOSE) cat.info(">>> returning to", wd)
    }
}


    # ------------------------------ #
    # GENERATE HTML REPORT FOR FGSEA #
    # ------------------------------ #

    ## Run the Bash script to generate HTML index files

report.scr.dir <- file.path(my.lib, "report")
wd <- getwd()       # so we can get back here when using setwd() below

cat("\n\tGenerating HTML Report for FGSEA\n\n")

fgsea.dir <- file.path(folder, "DESeq2", "GO_fgsea")

for (size.dir in list.dirs(fgsea.dir, recursive = F, full.names = T)) {
	
    grp.size <- basename(size.dir)
	cat("\n\tGene Set", grp.size, "\n")
	cat(  "\t=====================")
    for (sample.dir in list.dirs(size.dir, recursive = F, full.names = T)) {
    
    	sample <- basename(sample.dir)
        cat("\n",sample, "\n\n", sep = "")

        cur.dir <- setwd(sample.dir)
        if (VERBOSE) cat.info(">>> working in", cur.dir, '\n')
		system(paste("bash", file.path(report.scr.dir, "make-fgsea-index.sh")))
        file.copy(file.path(report.scr.dir, "style.css"), ".")
		cur.dir <- setwd(wd)    
        if (VERBOSE)  cat.info(">>> returning to", wd, '\n')
    }
	
    cur.dir <- setwd(size.dir)
    if (VERBOSE) cat.info(">>> working in", cur.dir, '\n')
	system(paste("bash", file.path(report.scr.dir, "rename-to-export.sh")))
	file.copy(file.path(report.scr.dir, "style.css"), ".")
	cur.dir <- setwd(wd)
    if (VERBOSE) cat.info(">>> returning to", wd, '\n')
}


# close all log files
sink.titanic()


