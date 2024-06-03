library(ensembldb)		# needs to be first to avoid S4 inconsistencies
library(Rsubread)		# for read mapping
library(AnnotationHub)		# to seek annotation packages
library(AnnotationForge)	# to build our annotation package
library(GO.db)
library(PFAM.db)

library("biomaRt")		# an alternate approach to retrieve annotation

library(tibble)			# general tools
library(tidyr)
library(stringr)
library(dplyr)
library(readr)
library(keypress)

library(edgeR)			# RNAseq with edgeR
library(limma)
library(RColorBrewer)
library(gplots)

library(DESeq2)			# RNAseq with DESeq2
library(IHW)			# for p-value adjustment with IHW
library(ggplot2)
library(ashr)

library(GSEABase)
library(fgsea)
library(reactome.db)
library(pathview)

library(gprofiler2)
library(clusterProfiler)
require(DOSE)
library(enrichplot)
library(europepmc)

library(cluster) 
library(factoextra)
library(fpc)
library("NbClust")

library(tcltk)
library(gWidgets2)


ALIGN <- FALSE		# do we start from the aligned data?
BOTH <- TRUE
use.online.annotation <- TRUE
VERBOSE <- TRUE
INTERACTIVE <- FALSE

# we will move inside './data/.' to work, so all paths are relative to it.
reference <- 'Cjaponica'               # coturnix
release <- "Coturnix_japonica_2.0"
target.organism <- 'Coturnix_japonica'
ncbi.taxid <- '93934'
ncbi.genus <- 'Coturnix'
ncbi.species <- 'japonica'
ncbi.version <- '2.0'
alignment.dir <- 'Alignments_Coturnix'
ens.db.pkg <- "EnsDb.Cjaponica"
ens.version <- "105"
mart.name <- 'cjaponica_gene_ensembl'
rnaseq.out <- 'rnaseq-test'
org.package <- "org.Cjaponica.eg.db"
KEGG_org <- 'cjo'

#reference <- "GRCg6a"			# gallus
reference <- 'Gg_GRGc7b'
release <- "GRCg6a"
target.organism <- 'Gallus_gallus'
ncbi.taxid <- '9031'
genus <- 'Gallus'
species <- 'gallus'
version <- '6a'
alignment.dir <- 'gg-aln-Rsubread-6a'
ens.db.pkg <- "EnsDb.Ggallus"
ens.version <- '106'
mart.name <- 'ggallus_gene_ensembl'
rnaseq.out <- 'rnaseq-test'
org.package <- "org.Ggallus.eg.db"
org.package <- "org.Gg.eg.db"
KEGG_org <- 'gga'

n.genes <- 1000		# number of top genes to revise

fastq.data <- 'fastq-qc-paired'
#data <- fastq.data
db.passwd <- '1_4_All+All_4_1'

out <-  alignment.dir

if (BOTH == T) {
    requireBothEnds <- T
    folder <- paste(rnaseq.out, "both_ends", sep='/')	# whether both ends should match or just any end
    logfile <- paste(folder, "log", "RNAseq.log", sep='/')
} else {
    requireBothEnds <- F
    folder <- paste(rnaseq.out, "any_end", sep='/')
    logfile <- paste(folder, "log", "RNAseq.log", sep='/')
}

# create needed directory hierarchy
dir.create(rnaseq.out, showWarnings=FALSE)
dir.create(folder, showWarnings=FALSE)
dir.create(file.path(folder, "log"), showWarning=FALSE)
dir.create(file.path(folder, "img"), showWarning=FALSE)
dir.create(file.path(folder, "edgeR"), showWarning=FALSE)
dir.create(file.path(folder, "edgeR/img"), showWarning=FALSE)
dir.create(file.path(folder, "edgeR/go"), showWarning=FALSE)
dir.create(file.path(folder, "edgeR/pfam"), showWarning=FALSE)
dir.create(file.path(folder, "edgeR/cluster"), showWarning=FALSE)
dir.create(file.path(folder, "edgeR/cluster/kmeans"), showWarning=FALSE)
dir.create(file.path(folder, "edgeR/cluster/pam"), showWarning=FALSE)
dir.create(file.path(folder, "edgeR/cluster/dbscan"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/raw"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/shrunk"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/signif"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/go"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/pfam"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/img"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/cluster"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/cluster/kmeans"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/cluster/pam"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/cluster/dbscan"), showWarning=FALSE)

as.png <- function(PLOT=NULL, 
                   file='out.png', width=1024, height=1024, 
                   overwrite=FALSE) {
    
    if (is.null(file)) {
        tryCatch(print(PLOT))
    } else if ( overwrite | ! file.exists(file)  ) {
        if (VERBOSE)
            cat("as.png(): creating", file, '\n')
        tryCatch( {
                png(file, width=width, height=height)
                print(PLOT)
            },
            finally=dev.off()
        )
    }
    return ()
}


align.fastq <- function(path, reference, aln.out, save.dir) {
    # get the fastq file names
    R1.fastq.files <- list.files(path=fastq.data, pattern='R1', full.names=TRUE)
    R2.fastq.files <- list.files(path=fastq.data, pattern='R2', full.names=TRUE)
    print(R1.fastq.files)
    print(R2.fastq.files)
    
    # we'll check if the output files exist to avoid repeating
    # work already done
    setwd(aln.out)
    dbname <- reference
    reference.fasta <- paste("./ref/", reference, ".fna", sep='')
    
    if (! file.exists(paste(reference, '.0.b.tab', sep=''))) {
        cat('\nBUILDING INDEX\n')
        # build the reference index inside the 'aln.out' directory
        buildindex(basename=dbname,reference=reference.fasta)
        dir()
    }
    
    r11 <- basename(R1.fastq.files[1])
    if (! file.exists(paste(aln.out, '/', r11, '.subread.BAM', sep=''))) {
        cat('\nALIGNING\n')
        # align the reads
        # IMPORTANT NOTE: WE NEED TWO FILES LISTING ALL THE FASTQ FILES TO ALIGN
        #	R1.fastq.files and R2.fastq.files
        align(index=dbname, readfile1=R1.fastq.files, readfile2=R2.fastq.files)
        
        # align will generate the output in the data directory, we
        # do not want to mix datasets, so we move the alignment results
        # to the corresponding output directory
        system(paste("mv ", fastq.data, "/*.BAM .", sep=""))
        system(paste("mv ", fastq.data, "/*.vcf .", sep=""))
        system(paste("mv ", fastq.data, "/*.summary .", sep=""))
    }
    setwd('..')
    
    # get the bam file names and inspect them to count the number
    # of reads that map to each genome position
    if ( ! file.exists(paste(save.dir, 'bam_files_stats.txt', sep='/'))) {
        bam.files <- list.files(path=aln.out, pattern='.BAM$', full.names = TRUE)
        bam.files
        props <- propmapped(files=bam.files)
        props
        write.table(props, file=paste(save.dir, 'bam_files_stats.txt', sep='/'), 
    	    row.names=T, col.names=T, sep='\t')
    }

}


compute.feature.counts <- function(files, reference, requireBothEnds, save.dir) {

    reference.gtf   <- paste("./ref/", reference, ".gtf", sep='')
    reference.gff   <- paste("./ref/", reference, ".gff", sep='')

    fc <- featureCounts(files=bam.files, 
	        annot.ext=reference.gtf, 
            isGTFAnnotationFile=T,
	        isPairedEnd=T, 
            requireBothEndsMapped=requireBothEnds, 
            primaryOnly=T, 
            ignoreDup=T, 
            useMetaFeatures=T)
    #  this means: process all BAM files
    #	Use as annotation the external file reference.gtf which is GTF
    #	BAM files contain paired end reads and we will only consider those
    #	where both ends match against the genome
    #    -------------------------------------------------------------
    #          ----->        <-----
    #	We will filter matches so that if a read may match more than one
    #	place in the genome, we will only consider the primary match and
    #	ignore other, duplicate matches
    #	We use meta-features
    #	We match against any feature in the GTF file, not only genes
    #	(this implies we will need to check the annotation carefully later)
    #	We do not remove chimeric fragments, but do actually count them too

    # SAVE FEATURE COUNTS
    # -------------------
    # now the variable fc contains 4 columns: annotation, target, counts and stats
    # but they exist only in the RAM memory, they are not stored somewhere safe,
    # so, we save them and create new variables so as to make it easier for us to 
    # manipulate the data
    write.csv(fc$counts, file=paste(save.dir, 'featureCounts.csv', sep='/'),
	      row.names=T, col.names=T, quote=F)
    write_delim(fc$stat, file=paste(save.dir, 'featureCounts_stat.txt', sep='/'), 
	        delim='      ')
    # save all the contents of 'fc' in an RDS file
    saveRDS(fc, file=paste(save.dir, '/featureCounts.rds', sep=''))
    #	'fc' can later be recovered with: fc <- readRDS(file='featureCounts.rds'
    # and save as well as Rdata file
    save(fc, file=paste(save.dir, '/featureCounts.RData', sep=''))
    #	'fc' can later be recovered with: fc <- load(file='featureCounts.RData')
}



##############################################################################
#
#  PREPARE ANNOTATION SO IT IS READY WHEN WE GET THE RESULTS
#
##############################################################################
#  THIS NEEDS GENERALIZING ### JR ###
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

cat("\n================================================\n")
cat(  "     B U I L D I N G   A N N O T A T I O N\n")
cat(  "================================================\n")

ah <- AnnotationHub()
# get ENSEMBL annotations for our target
#qr <- query(ah, c("EnsDb", target.organism))
# choose the most recent (last) one
#last.ref <- names(qr)[length(names(qr))]
#net.edb <- qr[[last.ref]]
#net.edb <- qr[["AH97993"]]	# Cjaponica v105

# with this we do not need to build it from GTF/GFF/mysql,

if (use.online.annotation == FALSE) {
    ###
    ### JR ### THIS NEEDS TO BE HAND TUNED !!!
    ###
    # We'll try to build an EnsDb Package nevertheless
    # so we can keep it locally for future use
    #EnsDbPackageDir <- paste(folder, 'EnsDb.Ggallus.v106', sep='/')
    #EnsDbPackageDir <- paste(folder, 'EnsDb.Cjaponica.v105', sep='/')
    EnsDbPackageDir <- paste(folder, '/', 
        ens.db.pkg, '.v', ens.version, sep='')
    
    scr <- system.file("scripts/generate-EnsDbs.R", package = "ensembldb")
    source(scr)

    # use our local GTF file
    if ( ! dir.exists(EnsDbPackageDir) ) {
        # generate SQLite database in place from GTF. Produces
        #    ./gallus_gallus.GRCg6a.106.sqlite
        gtf.edb <- ensDbFromGtf(gtf=reference.gtf, 
	        organism=target.organism,
                genomeVersion=release,
                version=ens.version,
                destDir=folder)			# lacks entrezid
        dir()
        # file name of the SQLite database created
        #DBFile <- paste(folder, 'gallus_gallus.GRCg6a.106.sqlite', sep='/')
        sqlite <- paste(target.organism, ".", release, ".", ens.version, ".sqlite", sep='/')
        DBFile <- paste(folder, '/', sqlite)

        # we'll select the SQLite database file generated from the GTF to
        # create an EnsDb R package 
        system(paste("mv", sqlite, DBFile))

        makeEnsembldbPackage(DBFile, version=ens.version, 
    		             maintainer="J. R. Valverde <jrvalverde@cnb.csic.es>", 
    		             author="J. R. Valverde",
                             destDir=folder, license="Artistic-2.0")
    } else {
        sqlite <- paste(target.organism, ".", release, ".", ens.version, ".sqlite", sep='/')
        DBFile <- paste(folder, '/', sqlite, sep='')
    }

    # retry with GFF if that didn't work
    if ( ! dir.exists(EnsDbPackageDir) ) {
        # As an alternative, we may use the GFF3 file (which should contain
        # the same information
        gff.edb <- ensDbFromGff(gff=reference.gff, 
           organism=target.organism,
                genomeVersion=release,
                version=ens.version)                     # lacks entrezid and transcript_id
        dir()
        # file name of the SQLite database created
        sqlite <- paste(target.organism, ".", release, ".", ens.version, ".sqlite", sep='/')
        DBFile <- paste(folder, '/', sqlite)

        # we'll select the SQLite database file generated from the GFF to
        # create an EnsDb R package 
        system(paste("mv", sqlite, DBFile))

        makeEnsembldbPackage(DBFile, version=ens.version, 
    		             maintainer="J. R. Valverde <jrvalverde@cnb.csic.es>", 
    		             author="J. R. Valverde",
                             destDir=folder, license="Artistic-2.0")
    } else {
        sqlite <- paste(target.organism, ".", release, ".", ens.version, ".sqlite", sep='')
        DBFile <- paste(folder, '/', sqlite, sep='')
    }

    # we still have a third option, build it from a MySQL file
    if ( ! file.exists(EnsDbPackageDir) ) {
        # This shouln't be needed because we have already done it above
        # using the GTF/GFF3 file.
        # This is an alternate way to do generate the EnsDb package from
        # MySQL data downloaded from ENSEMBL
        local.mysql.db <- paste("mysql/", 
        	tolower(target.organism), 
                "_core_", ens.version, "_6", sep='')
        #local.mysql.db <- "mysql/coturnix_japonica_core_104_2"
        #local.mysql.db <- "mysql/gallus_gallus_core_105_6"
        createEnsDbForSpecies(ens_version = ens.version,
                user = "sci", pass=db.passwd,
                host = "localhost", 
                local_tmp=local.mysql.db, 
                #species="gallus_gallus", 
                species=target.organism,
                dropDb=FALSE)
        # This should create a .sqlite file like ensDbFromGtf above,
        # named
        #	Gallus_gallus.GRCg6a.106.sqlite
        #       ^
        # with similar contents to the one produced from ensDBFromGtf/Gff;
        # if none of these does work, then tweak the process by hand in
        # source("scripts/generate-EnsDBs.R")
        sqlite <- paste(target.organism, ".", release, ".", ens.version, ".sqlite", sep='/')
        DBFile <- paste(folder, '/', sqlite)
        system(paste("mv", sqlite, DBFile))
        makeEnsembldbPackage(DBFile, version=ens.version, 
    		             maintainer="J. R. Valverde <jrvalverde@cnb.csic.es>", 
    		             author="J. R. Valverde",
                             destDir=folder, license="Artistic-2.0")

    } 
}


# SELECT ENSEMBL ANNOTATION DATABASE TO USE
if (use.online.annotation == TRUE) {
    dir.create('net.EnsDb', showWarnings=FALSE)
    if ( ! file.exists(paste('net.EnsDb', '/net.ens.db.RData', sep='')) ) {
        # get ENSEMBL annotations for our target
        qr <- query(ah, c("EnsDb", target.organism))
        # choose the most recent (last) one
        last.ref <- names(qr)[length(names(qr))]
        net.edb <- qr[[last.ref]]
        ens.db <- net.edb

        # save for future reference (downloading is too slow)
        saveRDS(ens.db, file=paste('net.EnsDb', '/net.ens.db.rds', sep=''))
        #	'ens.db' can later be recovered with: fc <- readRDS(file='net.ens.db.rds'
        # and save as well as Rdata file
        save(ens.db, file=paste('net.EnsDb', '/net.ens.db.RData', sep=''))
        #	'ens.db' can later be recovered with: fc <- load(file='net.ens.db.RData')
    } else {
        ens.db <- readRDS(file=paste('net.EnsDb', '/net.ens.db.rds', sep=''))
    }
} else {
    ens.db <- EnsDb(DBFile)
}

# check it
print(ens.db)
head(keys(ens.db, 'GENEID'))
columns(ens.db)


get_ens_ann <- function(gene_names, ens.db, folder='.') {
    out.file <- paste(folder, '/ensembl.annotation.txt', sep='')
    out.file.1 <- paste(folder, '/ensembl.annotation.1st.txt', sep='')
    if ( file.exists(out.file)) {
        ens.ann <- read.table(file=out.file, 
	                       sep='\t', header=T)
        ens.ann.1 <- read.table(file=out.file.1, 
	                       sep='\t', header=T)

    } else {
        ens.ann <- ensembldb::select(ens.db, 
                      column='GENID', keytype='GENEID', keys=gene_names, 
                      columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION', # no longer available
                                 'GENENAME', 'GENEID', 'ENTREZID', # empty
                                 'TXID', 'TXBIOTYPE', # these make the call fail
                                 'PROTEINID', 'UNIPROTID' # no longer available
                                )) 
        # SAVE ANNOTATION
        # ---------------
        write.table(ens.ann , file=out.file, 
	        sep='\t', row.names=T, col.names=T)

        # check if the amount of genes we have is the same as the number of 
        # the annotations that we have extracted
        if ( ! table(ens.ann$GENEID==rownames(fit)) ) {
            cat("Annotation does not match genes\n")
            cat("Using only one entry per gene\n")
            ens.ann.1 <- ens.ann[ ! duplicated(ens.ann$GENEID), ]
        } else {
            ens.ann.1 <- ens.ann
        }
        write.table(ens.ann.1, file=paste(folder, '/ensembl.annotation.1st.txt', sep=''), 
	        sep='\t', row.names=T, col.names=T)
    }
    return( list(ens.ann=ens.ann, ens.ann.1=ens.ann.1) )

}

# these should both give similar results
get_ensembl_ann <- function(gene_names) {
    return( ensembldb::select(ens.db, 
              keytype= 'GENEID', keys=rownames(fit), 
              columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION',
                         'GENENAME', 'GENEID', 'ENTREZID', 
                          'TXNAME', 'TXBIOTYPE', 
                          'PROTEINID', 'UNIPROTID')) )

} 



###################################################
###   A L T E R N A T E   A N N O T A T I O N   ###
###################################################

# ---------------------------------------------------------------
# M A K E     O R G     P A C K A G E
# ---------------------------------------------------------------

org.db <- NULL

if (use.online.annotation == TRUE) {
    # use AnnotationHub to seek a suitable package
    qo <- query(ah, "Orgdb")
    if ( last.ref %in% names(qo) ) {
        org.db <- query(ah, "Orgdb")[[ last.ref ]]
    } 
} 

if (is.null(org.db)) {
     # then try to build off-line annotation
     if ( ! require(org.package, character.only=T)) {
        #-----------------------------------------------
        # Create Org Package
        #-----------------------------------------------
        #
        # The datasets were first downloaded by hand from NCBI
        # and SwissProt into org.Gg.eg.db.
        # Then the following command had to be used:
        #
        makeOrgPackageFromNCBI(
                author  = "J. R. Valverde <jrvalverde@cnb.csic.es>", 
                maintainer = "J. R. Valverde <jrvalverde@cnb.csic.es>", 
                tax_id  = ncbi.taxid, # from NCBI Taxonomy browser (ncbi:txid9031)
                genus   = ncbi.genus, 
                species = ncbi.species, 
                version = ncbi.version, 
                outputDir = "./org", 
                #NCBIFilesDir=".", rebuildCache=FALSE)
                NCBIFilesDir = "./ncbi"#, 
                #rebuildCache=FALSE
                )
        # We specify a firectory to save locally the files used (and 
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
        }
        # This fails in Gallus gallus due to download failures from NCBI.
        # Should use the latest GitHub version installed with
        # library(devtools)
        # install_github("Bioconductor/AnnotationHub")
    }
}

# org.package has already been loaded
# convert the packahe name to variable name and save the variable
org.db <- get(org.package)

# at this point ens.db should contain the ENSEMBL data and 
# *.org.db the Org type data.


# ---------------------------------------------------------------
# O B T A I N   B I O M A R T   A N N O T A T I O N
# ---------------------------------------------------------------
if ( ! file.exists(paste(folder, 'biomaRt.annotation.1st.txt', sep='/')) ) {
    if ( ! file.exists(paste(folder, 'biomaRt.annotation.txt', sep='/')) ) {

        marts <- listMarts()
        head(marts)
        datasets <- listDatasets(useMart("ensembl"))
        head(datasets)

        ## set up connection to ensembl database
        ensembl=useMart("ENSEMBL_MART_ENSEMBL")

        # list the available datasets (species)
        listDatasets(ensembl) %>%  filter(str_detect(description, release))

        # already defined above
        #mart.name <- 'cjaponica_gene_ensembl'
        #mart.name <- 'ggallus_gene_ensembl'

        mart.db <- useMart("ensembl", mart.name)
        attributes <- listAttributes(mart.db)
        head(attributes)
        filters <- listFilters(mart.db)
        head(filters)

        # check the available "filters" - things you can filter for
        listFilters(mart.db) %>% 
            filter(str_detect(name, "ensembl"))

        # check the available "attributes" - things you can retreive
        attr <- listAttributes(mart.db)[,1:2]	# the first 200 are general
        listAttributes(mart.db)[,1:2] %>% 
            head(20)

        # we cannot get all the annotation at once because it times out
        #full.annot <- getBM(attributes=
        #                       c("ensembl_gene_id", "ensembl_transcript_id", 
        #			   "start_position", "end_position", 
        #                          "chromosome_name", "gene_biotype", 
        #                          "description", 
        #                          "entrezgene_id", "entrezgene_accession", "entrezgene_description", 
        #                          "go_id", "name_1006", "definition_1006", "go_linkage_type", "namespace_1003", 
        #                          "goslim_goa_accession", "goslim_goa_description", 
        #                          "pdb", 
        #                          "reactome", "uniprotswissprot"), 
        #                       mart=mart.db)

        # so we will retrieve the data in pieces, including ensembl_gene_id in
        # each piece so we can use it as key for merging the annotation
        if ( ! file.exists( paste(folder, '/biomart.ensembl.tab', sep='')) ) {
            bm.ensembl.annot <- getBM(attributes=c(
                                   "ensembl_gene_id", 
                                   "ensembl_transcript_id", 
                                   "start_position", "end_position", 
                                   "chromosome_name", "gene_biotype", 
                                   "description"),
                                mart=mart.db)
            write.table(bm.ensembl.annot, paste(folder, '/biomart.ensembl.tab', sep=''), 
	            row.names=T, col.names=T, sep='\t')
        } else {
            bm.ensembl.annot <- read.table(paste(folder, '/biomart.ensembl.tab', sep=''), 
	            header=T, sep='\t')
        }

        if ( ! file.exists(paste(folder, '/biomart.entrez.tab', sep='')) ) {
            bm.entrez.annot <- getBM(attributes=c(
                                   "ensembl_gene_id", 
                                    "entrezgene_id", "entrezgene_accession", "entrezgene_description"),
                                mart=mart.db)
            write.table(bm.entrez.annot, paste(folder, '/biomart.entrez.tab', sep=''), 
	            row.names=T, col.names=T, sep='\t')
        } else {
            bm.entrez.annot <- read.table(paste(folder, '/biomart.entrez.tab', sep=''), 
	            header=T, sep='\t')
        }

        if ( ! file.exists(paste(folder, '/biomart.go.tab', sep='')) ) {
            bm.go.annot <- getBM(attributes=c(
                                   "ensembl_gene_id", 
                                    "go_id", "name_1006", "definition_1006", "go_linkage_type", "namespace_1003"), 
                               mart=mart.db)
            write.table(bm.go.annot, paste(folder, '/biomart.go.tab', sep=''), 
	            row.names=T, col.names=T, sep='\t')
        } else {
            bm.go.annot <- read.table(paste(folder, '/biomart.go.tab', sep=''), 
	            header=T, sep='\t')
        }

        if ( ! file.exists(paste(folder, '/biomart.goslim.tab', sep='')) ) {
            bm.goslim.annot <- getBM(attributes=c(
                                   "ensembl_gene_id", 
                                    "goslim_goa_accession", "goslim_goa_description"),
                               mart=mart.db)
            write.table(bm.goslim.annot, paste(folder, '/biomart.goslim.tab', sep=''), 
	            row.names=T, col.names=T, sep='\t')
        } else {
            bm.goslim.annot <- read.table(paste(folder, '/biomart.goslim.tab', sep=''), 
	            header=T, sep='\t')
        }

        if ( ! file.exists(paste(folder, '/biomart.fam.tab', sep='')) ) {
            bm.fam.annot <- getBM(attributes=c(
                                   "ensembl_gene_id",
                                   "pfam",
                                   "pirsf",                                   "prints",
                                   "tigrfam"
                                   ), 
                               mart=mart.db)
            write.table(bm.fam.annot, paste(folder, '/biomart.fam.tab', sep=''), 
	            row.names=T, col.names=T, sep='\t')
        } else {
            bm.fam.annot <- read.table(paste(folder, '/biomart.fam.tab', sep=''), 
	            header=T, sep='\t')
        }
        if ( ! file.exists(paste(folder, '/biomart.prosite.tab', sep='')) ) {
            bm.prosite.annot <- getBM(attributes=c(
                                   "ensembl_gene_id",
                                   "scanprosite",
                                   "pfscan" 
                                   ), 
                               mart=mart.db)
            write.table(bm.prosite.annot, paste(folder, '/biomart.prosite.tab', sep=''), 
	            row.names=T, col.names=T, sep='\t')
        } else {
            bm.prosite.annot <- read.table(paste(folder, '/biomart.prosite.tab', sep=''), 
	            header=T, sep='\t')
        }
        
        if ( ! file.exists(paste(folder, '/biomart.sfam.tab', sep='')) ) {
            bm.sfam.annot <- getBM(attributes=c(
                                   "ensembl_gene_id",
                                   "superfamily"
                                   ), 
                               mart=mart.db)
            write.table(bm.sfam.annot, paste(folder, '/biomart.sfam.tab', sep=''), 
	            row.names=T, col.names=T, sep='\t')
        } else {
            bm.sfam.annot <- read.table(paste(folder, '/biomart.sfam.tab', sep=''), 
	            header=T, sep='\t')
        }

        if ( ! file.exists(paste(folder, '/biomart.extra.tab', sep='')) ) {
            bm.extra.annot <- getBM(attributes=c(
                                   "ensembl_gene_id", 
                                   "pdb",
                                   #"reactome", 
                                   "uniprotswissprot"), 
                               mart=mart.db)
            write.table(bm.extra.annot, paste(folder, '/biomart.extra.tab', sep=''), 
	            row.names=T, col.names=T, sep='\t')
        } else {
            bm.extra.annot <- read.table(paste(folder, '/biomart.extra.tab', sep=''), 
	            header=T, sep='\t')
        }

        # Now that we have all the pieces, merge them all together
        # into a single annotation variable
        #bm.annot <- merge(bm.ensembl.annot, bm.entrez.annot, by="ensembl_gene_id")
        #bm.annot <- merge(bm.annot, bm.go.annot, by="ensembl_gene_id")
        #bm.annot <- merge(bm.annot, bm.goslim.annot, by="ensembl_gene_id")
        #bm.annot <- merge(bm.annot, bm.extra.annot,  by="ensembl_gene_id")

        # SAVE ANNOTATION
        # ---------------
        write.table(bm.annot, file=paste(folder, '/biomaRt.annotation.txt', sep=''), 
	        sep='\t', row.names=T, col.names=T)
    } else {
        bm.ensembl.annot <- read.table(paste(folder, '/biomart.ensembl.tab', sep=''), 
	        header=T, sep='\t')
        bm.entrez.annot <- read.table(paste(folder, '/biomart.entrez.tab', sep=''), 
	        header=T, sep='\t')
        bm.go.annot <- read.table(paste(folder, '/biomart.go.tab', sep=''), 
	        header=T, sep='\t')
        bm.goslim.annot <- read.table(paste(folder, '/biomart.goslim.tab', sep=''), 
	        header=T, sep='\t')
        bm.fam.annot <- read.table(paste(folder, '/biomart.fam.tab', sep=''), 
	        header=T, sep='\t')
        bm.prosite.annot <- read.table(paste(folder, '/biomart.prosite.tab', sep=''), 
	        header=T, sep='\t')
        bm.sfam.annot <- read.table(paste(folder, '/biomart.sfam.tab', sep=''), 
	        header=T, sep='\t')
        bm.extra.annot <- read.table(paste(folder, '/biomart.extra.tab', sep=''), 
	        header=T, sep='\t')

        # this is too big to use
        #bm.annot <- read.table(file=paste(folder, '/biomaRt.annotation.txt', sep=''), 
	#    sep='\t', header=T)
    }

    # One possible way to do it would be to filter the queries above
    #	to retrieve the annotation matching ensembl_ids
    #
    # We can set a field to use to filter the output data
    # Set the filter type and values
    #ourFilterType <- "ensembl_gene_id"
    # and the values to select from that field
    #filterValues <- rownames(fit)
    #
    # and then obtain the specified annotation from records that match the values
    # specified in the filter field
    #fit.bm.extra.annot <- getBM(attributes=c(
    #                       "ensembl_gene_id", 
    #                       "pdb",
    #                       "reactome", 
    #                       "uniprotswissprot"), 
    #                   mart=mart.db,
    #                   filters=ourFilterType,
    #                   values=filterValues)
    #                   
    # deduplicate selecting the first annotation
    #fit.bm.extra.annot.1 <- fit.bm.extra.annot[ ! duplicated(fit.bm.extra.annot$ensembl_gene_id), ]
    # then we would repeat this for each annotation subset and merge all of them
    # at the end...

    # or we could deduplicate everything first and match aftwerards
    bm.ensembl.annot.1 <- bm.ensembl.annot[ ! duplicated(bm.ensembl.annot$ensembl_gene_id), ]
    bm.entrez.annot.1 <- bm.entrez.annot[ ! duplicated(bm.entrez.annot$ensembl_gene_id), ]
    bm.go.annot.1 <- bm.go.annot[ ! duplicated(bm.go.annot$ensembl_gene_id), ]
    bm.goslim.annot.1 <- bm.goslim.annot[ ! duplicated(bm.goslim.annot$ensembl_gene_id), ]
    bm.extra.annot.1 <- bm.extra.annot[ ! duplicated(bm.extra.annot$ensembl_gene_id), ]

    bm.annot.1 <- merge(bm.ensembl.annot.1, bm.entrez.annot.1, by="ensembl_gene_id")
    bm.annot.1 <- merge(bm.annot.1, bm.go.annot.1, by="ensembl_gene_id")
    bm.annot.1 <- merge(bm.annot.1, bm.goslim.annot.1, by="ensembl_gene_id")
    bm.annot.1 <- merge(bm.annot.1, bm.extra.annot.1,  by="ensembl_gene_id")

    write.table(bm.annot.1, 
                file=paste(folder, '/biomaRt.annotation.1st.txt', sep=''), 
	        sep='\t', row.names=T, col.names=T)
    # we save it to avoid repeating this in the future

    # or even do it all at once?
    ##bm.annot.1 <- bm.annot[ ! duplicated(bm.annot$ensembl_gene_id), ]
    ##write.table(bm.annot.1, 
    ##            file=paste(folder, '/biomart.annotation.1st.txt', sep=''), 
    ##            sep='\t', row.names=T, col.names=T)
} else {
    bm.annot.1 <- read.table(
    		file=paste(folder, '/biomaRt.annotation.1st.txt', sep=''), 
	        sep='\t', 
                header=T)
}
# as a byproduct we know at this point annotation files have been saved



# And now we are ready with our reference info at hand...


#################################################################################
#
# Get the alignments and feature counts
#
#################################################################################

# Note: since we will be looking at infected cells with aberrant RNA viral genomes
# which might present recombinations of virus-host reads, we will restrict ourselves
# to matches in both ends.

if ( ALIGN == TRUE) {
    ### 
    ### JR ### This needs to be generalized
    ###
    out <- 'gg-aln-Rsubread-6a'
   # out <- 'Alignments_Coturnix'
    align.fastq(path=path, reference=reference, aln.out=out, save.dir=folder) 
} else {
    out <- 'gg-aln-Rsubread-6a'
    #out <- 'Alignments_Coturnix'
}

#bam.files <- list.files(path = out, pattern = '.bam$', full.names = TRUE)
bam.files <- list.files(path = out, pattern = '.BAM$', full.names = TRUE)

# get the feature counts
#	we use EnsEMBL annotation from release 6a
#	this will tae the position counts and match them against the 
#	genes annotated in the GFF3 file for the reference, obtaining
#	a list of genes and the number of reads that mapped to them
#	as an indicator of their expression level
cat('\nCOMPUTING FEATURE COUNTS\n')

if (! file.exists(paste(folder, '/featureCounts.rds', sep=''))) {
    compute.feature.counts(files=bam.files, 
                           reference=reference, 
                           requireBothEnds=requireBothEnds, 
                           save.dir=folder) 
} else {
    fc <- readRDS(file=paste(folder, '/featureCounts.rds', sep=''))
}

#################################################################
#
#################################################################
                                                                #
#################################################################
#
#################################################################
                                                                #
#################################################################

# ---------------------------------------------------------------
# A N A L Y S I S     W I T H     D E S E Q 2
# ---------------------------------------------------------------


countDataFile <- paste(folder, "/featureCounts.csv", sep='')

countData <- read.csv(countDataFile,
		row.names=1) %>%
		as.matrix()
countData <- countData[rowSums(countData)>1, ]
head(countData)

colData <- read.delim(paste(rnaseq.out, "SampleInfo.tab", sep='/'))

colData$PFU <- as.factor(colData$PFU)
colData$Src <- as.factor(colData$Src)
colData$viral.dose <- as.factor(colData$viral.dose) #coturnix

# before doing the analysis, we will calculate normalized counts and
# save them for our record
# edgeR-type TMM normalization depends on the comparison design
# try both and see if there are differences
# sort metadata
metadata <- rbind(colData[13:15,], colData[1:3,], colData[10:12,], colData[7:9,], colData[4:6,])
rownames(metadata) <- 1:15
# sort countdata
counts <- cbind(countData[,13:15], countData[,1:3], countData[,10:12], countData[,7:9], countData[,4:6])

dds.pfu <- DESeqDataSetFromMatrix(counts, metadata,  design=~PFU)
dds.pfu <- estimateSizeFactors(dds.pfu)	# add norm factors to dds.pfu
sizeFactors(dds.pfu)
normalized_counts <- counts(dds.pfu, normalized=TRUE)	# get normalized counts
write.table(normalized_counts, file="rnaseq-test/normalized_counts.tab",
	        sep='\t', row.names=T, col.names=T) # better if row.names=F
#
tc <- t(normalized_counts)
tc <- cbind(c(rep(0, 3), rep(0.1, 3), rep(1.0,3), rep(10.0,3), rep(100.0,3)), tc)
# Note: DESeq2 doesn't actually use normalized counts, rather it uses
# the raw counts and models the normalization inside the Generalized Linear
# Model (GLM). These normalized counts will be useful for downstream
# visualization of results, but cannot be used as input to DESeq2 or any other
# tools that peform differential expression analysis which use the negative
# binomial model.
#
# To use these counts directly for variable reduction, we should first 
# select significant genes to reduce the number of explanatory variables
# and we should make sure that in the design, wt is the first element in
# countData.


name <- paste(folder, "DESeq2/dds.DESeq2.PFU", sep='/')

if ( ! file.exists(paste(name, "rds", sep='.')) ) {

    dds.pfu <- DESeqDataSetFromMatrix(countData, colData,  design=~PFU, tidy=F)

    dds.pfu <- DESeq(dds.pfu)

    # SAVE DESEQ ANALYSIS
    # -------------------
    #name <- paste(folder, "/DESeq2/dds.DESeq2", sep='')
    save(dds.pfu, file=paste(name, "RData", sep='.'))
    #	'dds' can later be recovered with: 
    #		dds <- load(file=paste(folder, "/dds.DESeq2.RData", sep=''))
    saveRDS(dds.pfu,  file=paste(name, "rds", sep='.'))
    # read as follows:    
} else {
    dds.pfu <- readRDS(file=paste(name, "rds", sep='.'))
}

cat("DESeq2 analysis produced the following fit")
print(resultsNames(dds.pfu))


name <- paste(folder, "DESeq2/dds.DESeq2.viral.dose", sep='/')
if ( ! file.exists(paste(name, "rds", sep='.')) ) {

    dds.vd <- DESeqDataSetFromMatrix(countData, colData,  design=~viral.dose, tidy=F)

    dds.vd <- DESeq(dds.vd)

    # SAVE DESEQ ANALYSIS
    # -------------------
    #name <- paste(folder, "/DESeq2/dds.DESeq2", sep='')
    save(dds.vd, file=paste(name, "RData", sep='.'))
    #	'dds' can later be recovered with: 
    #		dds <- load(file=paste(folder, "/dds.DESeq2.RData", sep=''))
    saveRDS(dds.vd,  file=paste(name, "rds", sep='.'))
    # read as follows:    
} else {
    dds.vd <- readRDS(file=paste(name, "rds", sep='.'))
}

cat("DESeq2 analysis produced the following fit")
print(resultsNames(dds.vd))

# get ens.ann from above

ens.ann <- get_ens_ann(rownames(countData), ens.db, folder)$ens.ann

ens.ann.1 <- ens.ann[ ! duplicated(ens.ann$GENEID), ]

# Add annotation to dds to keep everything in one place
#dds$ens.annot.1 <- ens.ann.1
#dds$bm.annot.1 <- b .annot.1


# SAVE SELECTED COMPARISONS
# -------------------------
# We will concentrate on the results between adjacent experiments
# where wt - PC - P - 1.0 - 0.1

dds_compare_annotate_save <- function(dds, 
	column, 
        x, 
        y, 
        filterFun=ihw, 
        alpha=0.01, 
        ensembl.ann, 
        biomart.ann,
        outDir='.' ) {

    out.file.base <- paste("raw", column, x, "vs", y, sep='_')
    # obtain comparison results
    cmp <- results(dds, 
                   contrast=c(column, x, y),
                   filterFun=filterFun,
                   alpha=alpha)
    # convert to data frame and add row names (ensembl.gene.id) 
    # as an additional column
    cmp.df.a <- data.frame(ensembl.gene.id=rownames(cmp), cmp)
    # ensembl.ann matches all ENSMEBL-ID in the cmp results 
    # (but lacks description)
    cmp.df.a <- cbind(cmp.df.a, 
                      ensembl.ann[ match(cmp.df.a$ensembl.gene.id, ensembl.ann$GENEID), ])
    # bm.annot fails to annotate some entries (why?)
    cmp.df.a <- cbind(cmp.df.a, 
                      biomart.ann[ match(cmp.df.a$ensembl.gene.id, biomart.ann$ensembl_gene_id), ])

    # and now save
    # save summary
    sink(paste(outDir, '/DESeq2/', out.file.base, '_summary.txt', sep=''), split=T)
    summary(cmp)
    sink()
    # unnanotated results object
    saveRDS(cmp, paste(outDir, "/DESeq2/", out.file.base, ".rds", sep=""))
    # annotated results as data frame (table)
    write.table(cmp.df.a,
	    paste(outDir, "/DESeq2/", out.file.base, "_annotated.tab", sep=""),
	    row.names=T, col.names=T, sep='\t')
    # histogram plot
    out.png <- paste(folder, '/DESeq2/img/DESeq2_', out.file.base, '_hist.png', sep='')
    as.png( {
        margins <- par("mar")
        par(mar=c(5, 5, 5, 5))
        hist(cmp.df.a$pvalue, main=out.file.base, breaks=1/alpha, xlab="p-value")
        par(mar=margins)
    }, out.png )

    return (cmp.df.a)
}


dds <- dds.pfu		# comparisons by column "PFU"

for (a in levels(as.factor(target$PFU)) ) {
    for (b in levels(as.factor(target$PFU)) ) {
        print(paste(a, b))
        if (a == b) next
        dds_compare_annotate_save( 
                        dds,
                        column="PFU",
                        x=a, y=b,
                        filterFun=ihw, alpha=0.01,
                        ensembl.ann=ens.ann.1,
                        biomart.ann=bm.annot.1,
                        outDir=folder
                        )

    }
}



# SAVE GENES WITH SIGNIFICANT CHANGES
# -----------------------------------

plot_and_save <- function(dds, 
                          column, x, y, 
                          filterFun=ihw, alpha=0.01,
        		  ensembl.ann, biomart.ann,
                          outDir='.',
                          save=TRUE ) {
    # base output file name
    out.base <- paste(column, ':_', x, '_x_', y, sep='')
    
    # Do the comparison and get the raw results (with baseMean, log2FC, p, padjusted)
    raw <- results(dds, 
                   contrast=c(column, x, y),
                   filterFun=filterFun,
                   alpha=alpha)
    # raw contains the data after comparing x and y as a DESEq2 result object
    # convert to a data frame and add row names (ensembl.gene.id) 
    # as an additional column
    raw.df <- data.frame(ensembl.gene.id=rownames(raw), raw)
    # ensembl.ann matches all ENSMEBL-ID in the raw results 
    # (but lacks description)
    raw.df <- cbind(
                raw.df, 
                ensembl.ann[ match(raw.df$ensembl.gene.id, ensembl.ann$GENEID), ]
                )
    # bm.annot fails to annotate some entries (why?)
    raw.df <- cbind(
                raw.df, 
                biomart.ann[ 
                  match(raw.df$ensembl.gene.id, biomart.ann$ensembl_gene_id),
                   ]
                )

    # and now save
    # save summary
    sink(paste(outDir, '/DESeq2/raw/raw_', out.base, '_summary.txt', sep=''), split=T)
    summary(raw)
    sink()
    # unnanotated results object
    saveRDS(raw, paste(outDir, "/DESeq2/raw/raw_", out.base, ".rds", sep=""))
    # annotated results as data frame (table)
    write.table(raw.df,
	    paste(outDir, "/DESeq2/raw/raw_", out.base, "_annotated.tab", sep=""),
	    row.names=T, col.names=T, sep='\t')
    # histogram plot
    out.png <- paste(folder, '/DESeq2/img/DESeq2_raw2', out.base, '_hist.png', sep='')
    as.png( {
        margins <- par("mar")
        par(mar=c(5, 5, 5, 5))
        hist(raw.df$pvalue, main=out.base, breaks=1/alpha, xlab="p-value")
        par(mar=margins)
    }, out.png )
    
    
    # now we'll shrink the data to improve visualization and ranking
    shrunk.lfc <- lfcShrink(dds, contrast=c(column, x, y), type="ashr")
    if (save) {
        ofile <- paste(outDir, "/DESeq2/img/DESeq2_", out.base, "_raw+shrunk_MA.png", sep='')
        cat("    plotting", ofile, '\n')
    } else ofile <- NULL
    as.png( {
        dim <- par("mfrow")
        par(mfrow=c(1,2))
        # plotMA shows the log2 fold changes attributable to a given variable
        # over the mean of normalized counts for all the samples in the dataset
        # Points above alpha are colored, outliers are shown as directed triangles    
        plotMA(raw, alpha=alpha)
        # it is useful to look at the shrunk l2fc values, which removes the noise
        # associated with l2fc changes from low-count genes
        plotMA(shrunk.lfc, alpha=alpha)
        # after plotMA, one may identify interesting points interactively using
        # identify() and clicking on them:
        # idx <- identify(res$baseMean, res$log2FoldChange)
        # rownames(res[idx, ])
        #
        # Alternatively, looking at the plot and deciding which coordinates are
        # of interest, one can use, e.g. 
        # res_wt_vs_PC[ (res_wt_vs_PC$log2FoldChange) > 4) 
        #               & (res_wt_vs_PC$baseMean > 1000), ]
        # or
        # res_wt_vs_PC[  (abs(res_wt_vs_PC$log2FoldChange) > 4) 
        #              & (res_wt_vs_PC$baseMean > 1000), ]
        par(mfrow=dim)
    }, ofile )
    
    # prepare for ggplot:
    #	convert to data frame
    #	add rownames as two new columns GENEID and ensembl_gene_id
    #   add annotation from EnsDb using GENEID
    #   add biomaRt annotation using ensembl_gene_id
    #   rename some columns for plotting
    shrunk.lfc$ensembl_gene_id <- rownames(shrunk.lfc)
    ann.shrunk <- as.data.frame(shrunk.lfc) %>%
    rownames_to_column("GENEID") %>% 
    left_join(ensembl.ann, "GENEID") %>% 
    left_join(biomart.ann, "ensembl_gene_id") %>% 
    rename(logFC=log2FoldChange, FDR=padj)   
    
    # ggplot does not work inside a function, so this code seems useless
    if (FALSE) {
        if (save) {
            ofile <- paste(outDir, "/DESeq2/img/DESeq2_", out.base, "_l2FCshrunk.png", sep='')
            cat("    plotting", ofile, '\n')
        } else ofile <- NULL
        as.png( {
            ggplot(ann.shrunk, 
              aes(x = log2(baseMean), y=logFC),
              environment=environment()) + # this is supposed to make it work in a local env
                geom_point(aes(colour=FDR < alpha), shape=20, size=0.5) +
                geom_text(data=~top_n(.x, 10, wt=-FDR), aes(label=SYMBOL)) +
                labs(x="mean of normalised counts", y="log fold change")
        }, ofile )
    }
    
    if (save) {
        ofile <- paste(outDir, "/DESeq2/shrunk/shrunk_", out.base, ".rds", sep='')
        cat("    saving", ofile, '\n')
        saveRDS(shrunk.lfc, file=ofile)
        ofile <- paste(outDir, "/DESeq2/shrunk/shrunk_", out.base, "_annotated.tab", sep='')
        cat("    saving", ofile, '\n')
        write.table(ann.shrunk, 
                    file=ofile,
                    row.names=T, col.names=T, sep='\t')
    }
    
    # find statistically significant changes
    signif <- raw[ raw$padj < alpha, ]
    signif$abs_lfc <- abs(signif$log2FoldChange)
    if (save) {
        ofile <- paste(folder, '/DESeq2/signif/signif_', out.base, "_α<", alpha, ".tab", sep='')
        cat("    saving", ofile, '\n')
        write.table(signif, 
                file=ofile, 
                sep='\t', row.names=T, col.names=T)
    }
    
    # order the data, firstly by the decreasing absolute 
    # value of log2FoldChange and secondly by the increasing pvalue....
    srt <- signif[ order(signif$abs_lfc,
                        signif$padj,
                        decreasing=c(T, F)), ]
    if (save) {
        ofile <- paste(folder, "/DESeq2/signif/signif_sorted_", out.base, ".tab", sep='')
        cat("    saving", ofile,'\n')
        write.table(srt, 
                file=ofile, 
	        sep='\t', row.names=T, col.names=T)
    }
    # annotate the sorted significant data
    srt.df <- data.frame(ensembl.gene.id=rownames(srt), srt)
    srt.df <- cbind(srt.df, 
                ensembl.ann[ match(srt.df$ensembl.gene.id, ensembl.ann$GENEID), ])
    srt.df <- cbind(srt.df, 
                biomart.ann[ match(srt.df$ensembl.gene.id, biomart.ann$ensembl_gene_id), ])
    if (save) {
        ofile <- paste(folder, "/DESeq2/signif/signif_sorted_", out.base, "_annotated.tab", sep='')
        cat("    saving", ofile,'\n')
        write.table(srt.df, 
                file=ofile, 
	        sep='\t', row.names=T, col.names=T)
    }

    
    return(list(result=raw, 
                shrunk=shrunk.lfc, 
                shrunk.annot=ann.shrunk, 
                signif=srt, 
                signif.annot=srt.df))
}

# Generating this takes long. We should save it and load from file
# when we run the script
ds.data <- list()
#x <- 0
contr <- "PFU"
grps <- levels(colData[ , contr ])
for (a in grps) {
    for (b in grps) {
        print(paste(a, b))
        if (a == b) next
        res <- plot_and_save( 
                        dds, contr,
                        x=a, y=b,
                        filterFun=ihw, alpha=0.01,
                        ensembl.ann=ens.ann.1,
                        biomart.ann=bm.annot.1,
                        outDir=folder,
						#save=TRUE
                        )
        print(names(res))
	name <- paste(contr, "_", a, "_", b, sep='')
        #x <- x + 1
        #ds.data[[x]] <- res
        ds.data[[name]] <- res
        #names(ds.data)[x] <- name
        #stop()
    }
}
print(names(ds.data))


# -------------------------------
# Do Gene Set Enrichment Analysis
# -------------------------------

GO_fgsea <- function (ann.shrunk.lfc, 
		      max.size=250,
                      out.dir=paste(folder, 'go_fgsea', sep='/'), 
                      out.name='GO_fgsea',
                      use.description=TRUE,
                      top.n=20,
                      top.biblio=5,
                      verbose=FALSE) {
    # Do GSEA on GO terms using fgsea

    # Rank all genes on their fold change.
    #	Here we exclude genes for which we have no EntrezID and
    #	use shrunk LFC values
#    gseaDat <- filter(ann.shrunk.lfc, !is.na(ENTREZID))
    gseaDat <- filter(ann.shrunk.lfc, !is.na(GENEID))

    ranks <- gseaDat$lfc
    #names(ranks) <- gseaDat$ENTREZID
    names(ranks) <- gseaDat$GENEID
    head(ranks)
    #uranks <- ranks[!duplicated(sort(names(ranks)))]

    # plot all the ranked fold changes
    out.png <- paste(out.dir, '/', out.name, '.barplot.png', sep='')
    as.png(barplot(sort(ranks, decreasing=T)), out.png)
    
    # load pathways
    #pathways <- ann.shrunk.lfc$ENTREZID
    #pathways <- ann.shrunk.lfc$go_id
    #names(pathways) <- gseaDat$ENTREZID
    #names(pathways) <- gseaDat$GENEID
    #upathways <- pathways[!duplicated(sort(names(pathway)))]
    
    # create a list of go_id terms where each term contains a vector
    # of emsembl_gene_id in that term
    #   first recover the annotation (in case we are re-run and do not
    #   have it
    if ( ! exists(substitute(bm.go.annot)) ) {
        bm.go.annot <- read.table(paste(folder, '/biomart.go.tab', sep=''), 
	            header=T, sep='\t')
    }
    #if ( ! exists(substitute(bm.goslim.annot)) ) {
    #    bm.goslim.annot <- read.table(paste(folder, '/biomart.goslim.tab', sep=''), 
    #	            header=T, sep='\t')
    #}

    if (use.description == FALSE) {
        # split() will divide the gene-ids by their go-id
        pathways <- split(bm.go.annot$ensembl_gene_id, bm.go.annot$go_id)
    } else {
        # Do the same but with large gene ontology names
        # split() will divide the gene-ids by their go-name
        pathways.go <- split(bm.go.annot$ensembl_gene_id, bm.go.annot$name_1006)
    }
    # do analysis
    # the resulting table contains enrichment scores and p-values
    out.file <- sprintf("%s/go_fgsea_10-%d.RData", out.dir, max.size)
    if (file.exists(out.file)) {
        load(file=out.file)
    } else {
        fgseaRes <- fgsea(pathways=pathways.go, 
                          stats=ranks, 
		          minSize=10, 
		          maxSize=max.size, 
                          nPermSimple=100000 
                          )
        save(fgseaRes, file=out.file)
    }
    if (verbose == TRUE)
        head(fgseaRes[order(padj, -abs(NES)), ], n=10)

    if (verbose == TRUE) {
        # plot enrichment score
        sorted.fgsea.res <- fgseaRes[order(padj, -abs(NES)), ]
        sfr.names <- sorted.fgsea.res$pathway
        for (i in 1:top.n) {
            if (use.description == FALSE)
                descr <- bm.go.annot[bm.go.annot$go_id == sfr.names[i], "name_1006"]
            else
                descr <- sfr.names[i]
            print(
                plotEnrichment(pathways.go[[ sfr.names[i] ]], ranks) +
                    labs(title=descr)
                )
            ans <- readline("Press RETURN to continue: ")
            if (ans == "q") break
        }
    }
    # gsea table plot of top.n gene families
    #   top_n() is now deprecated
    topUp <- fgseaRes %>%
        filter(ES > 0) %>%
        top_n(top.n, wt=-padj)

    topDown <- fgseaRes %>%
        filter(ES < 0) %>%
        top_n(top.n, wt=-padj)

    topPathways <- bind_rows(topUp, topDown) %>%
        arrange(-ES)
    # last resort when "pos" is used    
    #topPathways <- sorted.fgsea.res[1:top.n, ]

    # do the plots and save descriptions
    out.file <- paste(out.dir, "/topUp.", top.n, '.txt', sep='')
    for (i in 1:top.n) {
        name <- topUp[i]$pathway
        if (use.description == FALSE)
            descr <- unique(bm.go.annot[bm.go.annot$go_id == name, "name_1006"])
        else
            descr <- name
        cat(i, descr, file=out.file, '\n', sep='\t', append=TRUE)
        #print(name)
        #out.png <- paste(out.dir, '/', out.name, '.topUp', i, '.enrichment.png', sep='')
        out.png <- sprintf("%s/%s.topUp.%03d.enrichment.png", 
        	out.dir, out.name, i)
        out.plt <-
            plotEnrichment(pathways.go[[ name ]] , ranks) +
            labs(title=descr)
        ggsave(out.png, out.plt)
    }
    out.file <- paste(out.dir, "/topDn.", top.n, '.txt', sep='')
    for (i in 1:top.n) {
        name <- topDown[i]$pathway
        if (use.description == FALSE)
            descr <- unique(bm.go.annot[bm.go.annot$go_id == name, "name_1006"])
        else
            descr <- name
        cat(i, descr, file=out.file, '\n', sep='\t', append=TRUE)
        #print(name)
        #out.png <- paste(out.dir, '/', out.name, '.topDown', i, '.enrichment.png', sep='')
        out.png <- sprintf("%s/%s.topDn.%03d.enrichment.png", 
        	out.dir, out.name, i)
        out.plt <-
            plotEnrichment(pathways.go[[ name ]] , ranks) +
            labs(title=descr)
        ggsave(out.png, out.plt)
    }

    out.png <- paste(out.dir, '/', out.name, '.GSEAtable.png', sep='')
    as.png(
        plotGseaTable(pathways.go[topPathways$pathway], 
                      ranks, 
                      fgseaRes, 
                      gseaParam = 0.5)
    , out.png, width=1024, height=100*top.n, overwrite=TRUE)
    
    
    # and now do a plot of the interest in citations during
    # the last ten years for the top 5 sets
    out.png <- paste(out.dir, '/', out.name, '.EUPMC.png', sep='')
    cur.year <- as.integer(format(Sys.Date(), "%Y"))
    terms <- topDown[1:top.biblio]$pathway
    as.png(
        pmcplot(terms, (cur.year-10):cur.year, proportion=FALSE),
        out.png)

    
    # finally, return fgseaRes
    return(fgseaRes)
}



GO_KEGG_clusterProfiler <- function(ann.shrunk.lfc, 
		      max.size=250,
                      out.dir=paste(folder, 'go_fgsea', sep='/'), 
                      out.name='GO_fgsea',
                      use.description=TRUE,
                      OrgDb = org.Cjaponica.eg.db,
                      kegg_organism = "cjo",	# (cjo = coturnix japonica)
                                             # (gga = gallus gallus)
                      top.n=10,
                      top.biblio=5,
                      verbose=FALSE) {
        
    gseaDat <- filter(ann.shrunk.lfc, !is.na(ENTREZID))
    ranks <- gseaDat$lfc
    names(ranks) <- gseaDat$ENTREZID
    ranks<-na.omit(ranks)

    s.ranks <- sort(ranks, decreasing=T)
    
    ### G O annotation
    #
    gse <- gseGO(geneList=s.ranks, 
                 ont ="ALL", 
                 keyType = "ENTREZID", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = OrgDb, 
                 pAdjustMethod = "fdr")
                 #pAdjustMethod = "none")
    if (dim(gse)[1] == 0) {
        # Try without correction issuing a warning
        gse <- gseGO(geneList=s.ranks, 
                 ont ="ALL", 
                 keyType = "ENTREZID", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = OrgDb, 
                 #pAdjustMethod = "fdr")
                 pAdjustMethod = "none")
        out.name <- paste(out.name, '.raw_p', sep='')
    }
    if (dim(gse)[1] == 0) {
        # give up
        cat('', file=paste(out.dir, '/',
                out.name, 'NO_ENRICHED_GO', sep=''))
#        return()
    } else {
        out.png <- paste(out.dir, '/', out.name, '.GOdotplot.png', sep='')
        as.png(
            dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)
            , out.png)

        out.png <- paste(out.dir, '/', out.name, '.GOemmaplot.png', sep='')
        as.png(
            emapplot(pairwise_termsim(gse), showCategory = top.n)
            , out.png)

        out.png <- paste(out.dir, '/', out.name, '.GOridgeplot.png', sep='')
        as.png(
            ridgeplot(gse) + labs(x = "enrichment distribution")
            , out.png)

        out.png <- paste(out.dir, '/', out.name, '.EUPMC.png', sep='')
        cur.year <- as.integer(format(Sys.Date(), "%Y"))
        terms <- gse$Description[1:top.biblio]
        pmcplot(terms, (cur.year-10):cur.year, proportion=FALSE)

        for (i in 1:dim(gse)[1] ) {
        out.png <- paste(out.dir, '/', 
                out.name, '.GOcnetplot.', i, '.', gse$ID[i], '.png', sep='')
        # categorySize can be either 'pvalue' or 'geneNum'
        as.png(
            cnetplot(gse, categorySize="pvalue", foldChange=s.ranks, 
                     showCategory=i)
            , out.png)

            # Use the `Gene Set` param for the index in the title, and as the value for geneSetId
            out.png <- paste(out.dir, '/', 
                    out.name, '.GOgseaplot.', i, '.', gse$ID[i], '.png', sep='')
            gseaplot(gse, by = "all", title = gse$Description[i], geneSetID = i)
        }
        out.file <- paste(out.dir, '/', out.name, '.topGO.tab', sep='')
        write.table(gse, file=out.file, row.names=T, col.names=T, sep='\t')
    }
    
    
    ### K E G G annotation
    # let's try with KEGG (from the ENTREZID which is the same a ncbi-genid)
    kse <- gseKEGG(geneList     = s.ranks,
               organism     = kegg_organism,
               #nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "fdr",
               keyType       = "ncbi-geneid",
               nPermSimple = 100000)
               
    if (dim(kse)[1] == 0) {
        # Try without correction issuing a warning
        kse <- gseKEGG(geneList     = s.ranks,
                   organism     = kegg_organism,
                   #nPerm        = 10000,
                   minGSSize    = 3,
                   maxGSSize    = 800,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "none",
                   keyType       = "ncbi-geneid",
                   nPermSimple = 100000)
        out.name <- paste(out.name, 'raw_p', sep='')
    }
    if (dim(kse)[1] == 0) {
        # give up
        cat('', file=paste(out.dir, '/',
                out.name, 'NO_ENRICHED_KEGG', sep=''))
#        return()
    } else {
        out.png <- paste(out.dir, '/', out.name, '.KEGGdotplot.png', sep='')
        as.png( 
            dotplot(kse, showCategory = 10, title = "Enriched Pathways" , 
                    split=".sign") + facet_grid(.~.sign)
            , out.png)

        out.png <- paste(out.dir, '/', out.name, '.KEGGemmaplot.png', sep='')
        as.png(
            emapplot(pairwise_termsim(kse))
            , out.png)


        # categorySize can be either 'pvalue' or 'geneNum'
        out.png <- paste(out.dir, '/', out.name, '.KEGGcnetplot.png', sep='')
        as.png(
            cnetplot(kse, categorySize="pvalue", foldChange=s.ranks)
            , out.png)

        out.png <- paste(out.dir, '/', out.name, '.KEGGridgeplot.png', sep='')
        as.png(
            ridgeplot(kse) + labs(x = "enrichment distribution")
            , out.png)

        cur.dir <- getwd()
        setwd(out.dir)
        for (i in 1:dim(kse)[1]) {
            # for each of the pathways in kse

            # Use the `Gene Set` param for the index in the title, and as the value for geneSetId
            out.png <- paste(out.dir, '/', 
                out.name, '.GSEAplot.', i, '.', kse$ID[i], '.png', sep='')
            gseaplot(kse, by = "all", title = kse$Description[i], geneSetID = i)

            # Produce the native KEGG plot (PNG)
            dme <- pathview(gene.data=s.ranks, pathway.id=kse$ID[i], species = kegg_organism)

            # Produce a different plot (PDF) (different from the previous one)
            dme <- pathview(gene.data=s.ranks, pathway.id=kse$ID[i], species = kegg_organism, kegg.native = F)
        }
        setwd(cur.dir)

        out.file <- paste(out.dir, '/', out.name, '.topKEGG.tab', sep='')
        write.table(kse, file=out.file, row.names=T, col.names=T, sep='\t')
    }
}


# we do not know what is the best upper limit, so we'll try several
ms <- 500 # default value in GSEA, should work as well as the others
for (ms in c(50, 100, 250, 500)) {
    for (cmp.name in names(ds.data)) {
        # for each comparison name
        cmp.data <- ds.data[[ cmp.name ]]
        cat('Doing fGSEA on', cmp.name, '\n') 
        ann.shrunk.lfc <- cmp.data[[ "shrunk.annot" ]]
        #out.dir <- paste(folder, "DESeq2/GO_fgsea", cmp.name, sep='/')
        out.dir <- sprintf("%s/DESeq2/GO_fgsea/max.size=%03d/%s",
                folder, ms, cmp.name)
        dir.create(out.dir, showWarning=FALSE, recursive=TRUE)
        gogsea <- GO_fgsea(ann.shrunk.lfc, 
                   max.size=ms,
                   out.dir=out.dir, 
                   out.name='fgsea',
                   use.description=TRUE)
        # try dotplot(gogsea), emapplot(gogsea)... etc from clusterProfiler
        # we should add fgsea results to the cmp.data list
        ds.data[[cmp.name]][['go.fgsea']] <- gogsea

        out.dir <- sprintf("%s/DESeq2/GO+KEGG_cProf/max.size=%03d/%s",
                folder, ms, cmp.name)
        dir.create(out.dir, showWarning=FALSE, recursive=TRUE)
        GO_KEGG_clusterProfiler(
                      ann.shrunk.lfc, 
		      max.size=ms,
                      out.dir=out.dir, 
                      out.name='cProf',
                      use.description=TRUE,
                      verbose=FALSE)
    }
}



continue.on.enter("You may want to maximize your terminal before continuing ")
options(width=200)
n <- 10
for (i in names(ds.data)) {
    cat("\nMost significant', n, 'genes for", i, '\n')
    print(head(ds.data[[i]]$signif.annot[ , c("log2FoldChange", "entrezgene_description")] ), n)
    continue.on.enter("Press [ENTER] to continue ")
}
options(width=80)
continue.on.enter("Done, you can restore your terminal now ")


# Plot counts of the gene with maximal l2FC (up or down)
threshold <- 0
for (n in names(ds.data)) {
    #res <- data[["PFU_wt_0.1"]]$signif
    res <- ds.data[[n]]$signif.annot
    up <- res[ res$log2FoldChange > threshold, ]
    down <- res[ res$log2FoldChange < -threshold, ]

    most.up <- rownames(up)[which.max(up$log2FoldChange)]
    cat("\ncounts of most overexpressed gene in", n, ":\n",
        most.up, up$log2FoldChange[which.max(up$log2FoldChange)], up$entrezgene_description[which.max(up$log2FoldChange)], 
        "\n")
    #plotCounts(dds, gene=rownames(res)[which.max(res$log2FoldChange)], intgroup="PFU")
    print(plotCounts(dds, gene=most.up, intgroup="PFU"))
    k <- continue.on.key()
    if (k == "q") break
    
    most.down <- rownames(down)[which.min(down$log2FoldChange)]
    cat("\ncounts of most underexpressed gene in", n, ":\n",
        most.down, down$log2FoldChange[which.min(down$log2FoldChange)], down$entrezgene_description[which.min(down$log2FoldChange)],
        "\n")
    #plotCounts(dds, gene=rownames(res)[which.min(res$log2FoldChange)], intgroup="PFU")
    print(plotCounts(dds, gene=most.down, intgroup="PFU"))
    continue.on.key()
    if (k == "q") break
}
cat("\n")

# Analyze GO representation
# -------------------------
#
# we use the biomaRt database from above
#	mart.db might have been already assigned above, but not
#	necessarily
#ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
mart.db <- useMart("ensembl", mart.name)
#listAttributes(mart.db)

for (n in names(ds.data)) {
    cat("Processing GOs for", n, '\n')
    sig <- as.data.frame(ds.data[[n]]$signif)
    sig.a <- as.data.frame(ds.data[[n]]$signif.annot)
    sig$ensembl_gene_id <- rownames(sig)
    ourFilterType <- "ensembl_gene_id"
    filterValues <- sig$ensembl_gene_id
    # this will give us ALL GO annotations
    # we should use the saved searches from above to save bandwidth.
#    gos <- getBM(attributes=c("ensembl_gene_id", 
#                              "go_id", "name_1006", 
#                              "definition_1006" ), 
#                 mart=mart.db,
#                 filters=ourFilterType,
#                 values=filterValues)
    gos <- bm.go.annot[ bm.go.annot$ensembl_gene_id %in% filterValues,
                        c("ensembl_gene_id", 
                          "go_id", "name_1006", 
                          "definition_1006" ) ]
    sig.a <- merge(sig, gos, by="ensembl_gene_id")
    table(sig.a$go_id)	# this gives counts, but we'd like to multiply those
                            # counts by log2FC
    #sig.a[ , c(1, 3, 10)]
    # aggregate log2FC by go_id and sum the values
    print(
    aggregate(sig.a$log2FoldChange, 
              by=list(Category=sig.a$go_id),
              FUN=sum)
    )
    # or with formula interface
    gosums <- aggregate(log2FoldChange ~ go_id, sig.a, sum)
    gosums <- gosums[ order(gosums$log2FoldChange), ]
    # annotate gosums
    gosums.a <- cbind(gosums, bm.annot.1[ match(gosums$go_id, bm.annot.1$go_id), ])
    
    out.dir <- paste(folder, '/DESeq2/go', sep='')
    dir.create(out.dir, showWarning=FALSE)
    out.file <- paste(out.dir, '/GO_', n, "_sum.tab", sep='')
    write.table(gosums.a, out.file, sep='\t')
    
    goavgs <- aggregate(log2FoldChange ~ go_id, sig.a, mean)
    goavgs <- gosums[ order(gosums$log2FoldChange), ]
    # annotate gosums
    goavgs.a <- cbind(goavgs, bm.annot.1[ match(goavgs$go_id, bm.annot.1$go_id), ])
    out.file <- paste(out.dir, '/GO_', n, "_average.tab", sep='')
    write.table(goavgs.a, out.file, sep='\t')
}


# Reannotate GO with ancestry
# library(GO.db)

annotate.go.ancestry <- function (annot.data, l2fc.threshold, outDir) {
    res <- annot.data
    
    goBPanc <- as.list(GOBPANCESTOR)
    # remove GO terms that do not have any ancestor
    goBPanc <- goBPanc[ ! is.na(goBPanc) ]

    goCCanc <- as.list(GOCCANCESTOR)
    # remove GO terms that do not have any ancestor
    goCCanc <- goCCanc[ ! is.na(goCCanc) ]

    goMFanc <- as.list(GOMFANCESTOR)
    # remove GO terms that do not have any ancestor
    goMFanc <- goMFanc[ ! is.na(goMFanc) ]
    
    res <- res[ abs(res$log2FoldChange) > l2fc.threshold, ]
    go.l2fc <- res[ , c("ensembl_gene_id", "log2FoldChange", "GENENAME", "go_id") ]
    
    for (i in 1:nrow(res)) {
        if ((i %% 100) == 0) cat(".")
        if (is.null(res[i, "go_id"])) next
        if (is.na(res[i, "go_id"])) next
        if (res[i, "go_id"] == '') next
        anc <- goBPanc[[ res[i, "go_id"] ]]
        if (! is.null(anc) & ! length(anc) == 0) {
            if (! is.na(anc) & anc != '') {
            # add  these GoIds to the table of gene/l2FC/GoIds 
                for (anc.goid in anc) {
                    go.l2fc <- rbind(go.l2fc, data.frame(
                          ensembl_gene_id=res[i, "ensembl.gene.id"],
                          log2FoldChange=res[i, "log2FoldChange"],
                          GENENAME=res[i, "GENENAME"],
                          go_id=anc.goid))
                }
            }
        }
        anc <- goCCanc[[ res[i, "go_id"] ]]
        if (! is.null(anc) ) {
            if (! is.na(anc) & anc != '') {
            # add  these GoIds to the table of gene/l2FC/GoIds 
                for (anc.goid in anc) {
                    go.l2fc <- rbind(go.l2fc, data.frame(
                          ensembl_gene_id=res[i, "ensembl.gene.id"],
                          log2FoldChange=res[i, "log2FoldChange"],
                          GENENAME=res[i, "GENENAME"],
                          go_id=anc.goid))
                }
            }
        }
        anc <- goMFanc[[ res[i, "go_id"] ]]
        if (! is.null(anc) ) {
            if (! is.na(anc) & anc != '') {
            # add  these GoIds to the table of gene/l2FC/GoIds 
                for (anc.goid in anc) {
                    go.l2fc <- rbind(go.l2fc, data.frame(
                          ensembl_gene_id=res[i, "ensembl.gene.id"],
                          log2FoldChange=res[i, "log2FoldChange"],
                          GENENAME=res[i, "GENENAME"],
                          go_id=anc.goid))
                }
            }
        }
    }
    # aggregate data by GOID (will result in c(Category, x) columns
    go.l2fc.sum <- aggregate(go.l2fc$log2FoldChange, 
               by=list(Category=go.l2fc$go_id), FUN=sum)
    go.l2fc.avg <- aggregate(go.l2fc$log2FoldChange, 
               by=list(Category=go.l2fc$go_id), FUN=mean)
    
    # rename Category, x to go_id, l2fc
    colnames(go.l2fc.sum) <- c('go_id', 'sum.l2fc')
    colnames(go.l2fc.avg) <- c('go_id', 'avg.l2fc')
  
    
    # annotate
    godesc <- AnnotationDbi::select(GO.db, keys=go.l2fc$go_id, 
                  columns= c("GOID", "TERM", "DEFINITION", "ONTOLOGY"), 
                  keytype= "GOID")
    go.l2fc <- cbind(go.l2fc, godesc)
    godesc <- AnnotationDbi::select(GO.db, keys=go.l2fc.sum$go_id, 
                  columns= c("GOID", "TERM", "DEFINITION", "ONTOLOGY"), 
                  keytype= "GOID")
    go.l2fc.sum <- cbind(go.l2fc.sum, godesc)
    godesc <- AnnotationDbi::select(GO.db, keys=go.l2fc.avg$go_id, 
                  columns= c("GOID", "TERM", "DEFINITION", "ONTOLOGY"), 
                  keytype= "GOID")
    go.l2fc.avg <- cbind(go.l2fc.avg, godesc)


    # sort by l2FC and save
    go.l2fc <- go.l2fc[ order(go.l2fc$log2FoldChange, decreasing=T), ]
    go.l2fc.sum <- go.l2fc.sum[ order(go.l2fc.sum$sum.l2fc, decreasing=T), ]
    go.l2fc.avg <- go.l2fc.avg[ order(go.l2fc.avg$avg.l2fc, decreasing=T), ]
    # sort by abs(log2FC)
    go.l2fc.sum.abs <- go.l2fc.sum[ order(abs(go.l2fc.sum$sum.l2fc), decreasing=T), ]
    go.l2fc.avg.abs <- go.l2fc.avg[ order(abs(go.l2fc.avg$avg.l2fc), decreasing=T), ]
    
    # save
    out.file <- paste(outDir, '/GOANC_', n, ".tab", sep='')
    write.table(go.l2fc, out.file, sep='\t')
    out.file <- paste(outDir, '/GOANC_', n, "_sum.tab", sep='')
    write.table(go.l2fc.sum, out.file, sep='\t')
    out.file <- paste(outDir, '/GOANC_', n, "_average.tab", sep='')
    write.table(go.l2fc.avg, out.file, sep='\t')
    out.file <- paste(outDir, '/GOANC_', n, "_sum_abs.tab", sep='')
    write.table(go.l2fc.sum.abs, out.file, sep='\t')
    out.file <- paste(outDir, '/GOANC_', n, "_average_abs.tab", sep='')
    write.table(go.l2fc.avg.abs, out.file, sep='\t')
    
}
 
 
l2fc.threshold <- 0
out.dir <- paste(folder, '/Deseq2/goanc')
dir.create(out.dir, showWarning=FALSE)

for (n in names(ds.data)) {
    cat("\nTracing GO ancestry for", n, '\n')
    res <- ds.data[[n]]$signif.annot
        
    annotate.go.ancestry(res, l2fc.threshold, out.dir)
    
    #ans <- readline("Press RETURN to continue: ")
    #if (ans == "q") break
}


# Analyze PFAM representation
# ---------------------------
#
# for PFAM, we can use
#
# Get PFAM database and use AC -> DE mapping
#library(PFAM.db)
db <- PFAMDE
mk <- PFAMDE[mappedkeys(PFAMDE)]
xx <- as.list(mk)
pfam.table <- toTable(PFAMDE)

out.dir <- paste(folder, '/DESeq2/pfam', sep='')
dir.create(out.dir, showWarning=FALSE)


#for (i in pfam.ids) print(xx[[i]])

# Get PFAM families and sort them by their representation in the dataset
options(width=200)
for (n in names(ds.data)) {
    sig <- as.data.frame(ds.data[[n]]$signif) 
    names(sig)
    sig$ensembl_gene_id <- rownames(sig)
    sig.a <- cbind(sig, 
           bm.fam.annot[ match(sig$ensembl_gene_id, bm.fam.annot$ensembl_gene_id), ])
    pfam.ids <- sig.a$pfam[ ! is.na(sig.a$pfam) ]
    aggregate(sig.a$log2FoldChange, by=list(Category=sig.a$pfam), FUN=sum)

    sig.a <- merge(sig, bm.fam.annot, by="ensembl_gene_id")
    head(sig.a, 10)
    
    sig.sum <- aggregate(sig.a$log2FoldChange, by=list(Category=sig.a$pfam), FUN=sum)
    sig.sum <- sig.sum[ order(sig.sum$x, decreasing=T), ]
    sig.sum <- sig.sum[ ! is.na(sig.sum$Category), ]
    sig.sum <- sig.sum[ sig.sum$Category != '', ]
    names(sig.sum) <- c("pfam", "sum.l2FC")
    # annotate
    for (i in 1:length(sig.sum$pfam)) { 
        pf <- sig.sum$pfam[i] ; 
        if (! is.null(xx[[pf]])) { 
            sig.sum$pfam.de[i] <- xx[[pf]] 
        } else { 
            sig.sum$pfam.de[i] <- '' 
        }
    }
    cat("\nMost represented PFAM families in", n, "\n")
    print(head(sig.sum, 20)) ; print(tail(sig.sum, 20))
    out.file <- paste(out.dir, '/PFAM_', n, "_sum.tab", sep='')
    write.table(sig.sum, out.file, sep='\t')

    sig.avg <- aggregate(sig.a$log2FoldChange, by=list(Category=sig.a$pfam), FUN=mean)
    sig.avg <- sig.avg[ order(sig.avg$x, decreasing=T), ]
    sig.avg <- sig.avg[ ! is.na(sig.avg$Category), ]
    sig.avg <- sig.avg[ sig.avg$Category != '', ]
    names(sig.avg) <- c("pfam", "avg.l2FC")
    # annotate
    for (i in 1:length(sig.avg$pfam)) { 
        pf <- sig.avg$pfam[i] ; 
        if (! is.null(xx[[pf]])) { 
            sig.avg$pfam.de[i] <- xx[[pf]] 
        } else { 
            sig.avg$pfam.de[i] <- '' 
        }
    }
    out.file <- paste(fout.dir, '/PFAM_', n, "_average.tab", sep='')
    write.table(sig.avg, out.file, sep='\t')
    
    # save also the files sorted by abs(l2fc)
    sig.sum.abs <- sig.sum[ order(abs(sig.sum$sum.l2FC), decreasing=T), ]
    sig.avg.abs <- sig.avg[ order(abs(sig.avg$avg.l2FC), decreasing=T), ]
    out.file <- paste(out.dir, '/PFAM_', n, "_sum_abs.tab", sep='')
    write.table(sig.sum.abs, out.file, sep='\t')
    out.file <- paste(out.dir, '/PFAM_', n, "_average_abs.tab", sep='')
    write.table(sig.avg.abs, out.file, sep='\t')


    #continue.on.key()
    ans <- continue.on.enter(prompt="Press return to continue ")
    if (ans == "q") break
}
cat('\n')
options(width=80)
