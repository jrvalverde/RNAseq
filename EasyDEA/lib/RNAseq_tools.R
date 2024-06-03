### FUNCTIONS FOR RNA-SEQ AND DIFFERENTIAL EXPRESSION ANALYSIS ###

#' create.output.hierarchy
#'
#' Create the hierarchy of output directories that we need for the analysis
#' 
#' @param	rnaseq.out	the folder where we want to store all our output
#' @param	use.both.reads	whether both reads of a paired-reads sequencing
#'			experiment should be required to map in the analysis
#' 
#' @return	the name of the output folder where results will be stored
#'
#' @usage	out.dir <- create.output.hierarchy(rnaseq.out, use.both.reads)
#' 
#' @examples	out.dir <- create.output.hierarchy('.', TRUE)
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
create.output.hierarchy <- function(rnaseq.out='.', use.both.reads=TRUE)
{
# set the name of the output folder and log file
    if (use.both.reads == T) {
        requireBothEnds <- T
        folder <- paste(rnaseq.out, "both_ends", sep='/')	# whether both ends sshould match or just any end
    } else {
        requireBothEnds <- F
        folder <- paste(rnaseq.out, "any_end", sep='/')
    }

    # create needed directory hierarchy
    # this could go into a separate function for simplicity
    dir.create(rnaseq.out, showWarnings=FALSE)
    dir.create(folder, showWarnings=FALSE)
    dir.create(file.path(folder, "log"), showWarning=FALSE)
    dir.create(file.path(folder, "img"), showWarning=FALSE)
    dir.create(file.path(folder, "annotation"), showWarning=FALSE)
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

    return(folder)
}

#' align.fastq
#'
#' Align all fastq files inside a folder against a reference genome
#' 
#' @param	fastq.dir	the path to the folder where the fastq files are stored
#'
#' @param	reference	the organism reference genome
#'
#' @param 	alignment.out	the path to a directory where we want to store the alignments
#' 
#' @return	nothing
#'
#' @usage	align.fastq(path, reference, aln.out)
#' 
#' @examples	align.fastq('./fastq', 'ref/Hsapiens.fna', 'alignments')
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#

# Align fastQ reads into the reference genome
align.fastq <- function(fastq.dir, reference, alignment.dir) {
   
	# Get the fastq file names
    R1.fastq.files <- list.files(path=fastq.dir, pattern='R1', full.names=TRUE)
    R2.fastq.files <- list.files(path=fastq.dir, pattern='R2', full.names=TRUE)
    print(R1.fastq.files)
    print(R2.fastq.files)
    
    # we'll check if the output files exist to avoid repeating
    # work already done
	
	ref.fasta <- reference
	ref.name <- sub("\\.[[:alnum:]]+$", "", basename(reference))

	# Go to the alignment directory   
	wd <- getwd()	
	setwd(alignment.dir)
	system(paste('ln -sfr', dirname(paste(wd, ref.fasta, sep = "/")), "ref"))
	
	# Check wether sorted bam files already exist.
	if (length(list.files(pattern='.sorted.bam$', ignore.case=T)) > 0){
		bam.files <- list.files(pattern='.sorted.bam$', ignore.case=T, full.name = F)
		cat('\nUSING EXISTING SORTED BAM FILES:')
		cat('\t', bam.files, sep='\n\t')
		setwd(wd)
		return(paste(alignment.dir, bam.files, sep="/"))
		
	# Check wether unsorted bam files exist.
	} else if (length(list.files(pattern='.bam$', ignore.case=T)) > 0){	
		bam.files <- list.files(pattern='.bam$', ignore.case=T, full.name = F)
		cat('\nUSING EXISTING BAM FILES:\n')
		cat(paste(bam.files), sep='\n\t')
		setwd(wd)
		return(paste(alignment.dir, bam.files, sep="/"))
	}

	# If Bam files do not exist, perform alignment with Rsubread
	if (! file.exists(paste('ref/', ref.name, '.00.b.tab', sep=''))) {

		# build the reference index inside the 'alignment.dir' directory
		cat('\nBUILDING INDEX\n')
		buildindex(basename = paste("ref", ref.name, sep = "/"),
					reference = paste(wd, ref.fasta, sep = "/"))

	} else {
		cat('\nUSING INDEXED REFERENCE\n')
	}    

	if (! file.exists(paste(basename(R1.fastq.files[1]), '.subread.BAM', sep=''))) {

		# Align the reads
		# IMPORTANT NOTE: WE NEED TWO FILES LISTING ALL THE FASTQ FILES TO ALIGN
		#	R1.fastq.files and R2.fastq.files

		cat('\nALIGNING USING R_SUBREAD\n')
		align(index = paste('ref', ref.name, sep='/'),
				#nthreads = nthreads,		##ADRIAN
				readfile1 = paste(wd, R1.fastq.files, sep = "/"),
				readfile2 = paste(wd, R2.fastq.files, sep = "/"))        

		# Align will generate the output in the fastq directory, we
		# so we move the alignment results to the output directory
		system(paste("mv ", paste(wd, fastq.dir, sep = "/"), "/*.BAM .", sep=""), ignore.stderr = T)
		system(paste("mv ", paste(wd, fastq.dir, sep = "/"), "/*.vcf .", sep=""), ignore.stderr = T)
		system(paste("mv ", paste(wd, fastq.dir, sep = "/"), "/*.summary .", sep=""), ignore.stderr = T)
	}

	bam.files <- list.files(pattern='.subread.bam$', ignore.case=T)

    
	# Alignment Statistics:
	# Get the bam file names and inspect them to count the number
	# of reads that map to each genome position
    
	if ( ! file.exists('bam_files_stats.txt') & ! length(list.files(pattern = '.idxstats$')) == length(bam.files)) {
        
		cat('\nGENERATING BAM FILES STATISTICS')
        props <- propmapped(files=bam.files)
        props
        write.table(props, file='bam_files_stats.txt', 
    	    row.names=T, col.names=T, sep='\t')
    }
	cat('READ ALIGNMENT - DONE')
	setwd(wd)
	return(paste(alignment.dir, bam.files, sep="/"))
}


#' compute.feature.counts
#'
#' Takes a list of aligned bam files and a reference and computes the times
#' each gene (feature) in the reference is matched by an aligned read. It
#' will try to use a file named 'reference'.gtf or, if one does not exist,
#' a file named 'reference'.gff as a source of annotation with the feature
#' (gene) coordinates. As such, the reference gtf/gff file must exist and
#' match the reference fasta sequence used in the alignment.
#' 
#' @param	bam.files	a list of bam files with reads mapped to the reference
#' @param	reference	the base name of the reference used for aligning and
#'				to be used to search for an annotation source
#' @param	requireBothEnds	whether both ends of a read pair are required to 
#'				match in order to be counted
#' @param	save.dir	the name of the directory where the feature
#'				counts will be saved for future reference
#' 
#' @return	the feature count table
#'
#' @usage	fc <- compute.feature.counts(
#'			bam.files=list.files(path = aln.dir, pattern = '.BAM$', full.names = TRUE),
#'			reference=ref.genome,
#'			requireBothEnds=T,
#'			save.dir='.')
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
compute.feature.counts <- function(bam.files, annotation, feature.count.dir, requireBothEnds = T) {

    ref.name  <- sub("\\.[[:alnum:]]+$", "", basename(annotation))
    ref.gtf   <- sub("\\.[[:alnum:]]+$", ".gtf", annotation)
    ref.gff   <- sub("\\.[[:alnum:]]+$", ".gff", annotation)
	
	# If feature count data already exists, load it
    if (file.exists(paste(feature.count.dir, "featureCounts.tab", sep = "/"))){
		
		fc <- read.delim(paste(feature.count.dir, "featureCounts.tab", sep = "/"), header = TRUE, row.names = 1, sep = "\t")
		cat(paste('\nUSING ALREADY EXISTING FEATURE COUNT FILE: ', feature.count.dir, "/featureCounts.tab\n", sep=""))
		return(fc)
	}
		
	# First, try to annotate counts using GTF file. If it is not found,
	# try with GFF. If neither of the annotation files exist, abort.
    if (! file.exists(ref.gtf) && ! file.exists(ref.gff)) {
        cat.err("Annotation files", ref.gtf, ".GTF or .GFF not found\n")    
	} else if (file.exists(ref.gtf)) {
		annot.ext <- ref.gtf
		isGTF <- TRUE
	} else {
		annot.ext <- ref.gff
		isGTF <- FALSE
	}
	
    fc <- featureCounts(files = bam.files, 
	    				annot.ext = annot.ext, 
        				isGTFAnnotationFile = isGTF,
	    				isPairedEnd = T, 
        				requireBothEndsMapped = requireBothEnds, 
        				primaryOnly = T, 
        				ignoreDup = T, 
        				useMetaFeatures = T)#, nthreads = nthreads) ##ADRIAN

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
	
	# Remove the xtension .bam from the feature count table header
	colnames(fc$counts) <- sub('\\..*$', '', colnames(fc$counts))
	fc$targets <- sub('\\..*$', '',fc$targets)
	
	# SAVE FEATURE COUNTS
	# -------------------
	# now the variable fc contains 4 columns: annotation, target, counts and stats
	# but they exist only in the RAM memory, they are not stored somewhere safe,
	# so, we save them and create new variables so as to make it easier for us to 
	# manipulate the data
	write.table(fc$counts, file=paste(feature.count.dir, 'featureCounts.tab', sep='/'),
				sep = "\t", row.names = TRUE, col.names = TRUE)
	write.csv(fc$counts, file=paste(feature.count.dir, 'featureCounts.csv', sep='/'))
	write.table(fc$stats, file=paste(feature.count.dir, 'featureCounts_stat.txt', sep='/'),
				sep = "\t", row.names = TRUE, col.names = TRUE)

	# save all the contents of 'fc' in an RDS file and in an Rdata file

	saveRDS(fc, file=paste(feature.count.dir, '/featureCounts.rds', sep=''))
	save(fc, file=paste(feature.count.dir, '/featureCounts.RData', sep=''))

	# 'fc' can later be recovered with:
	# 		fc <- readRDS(file='featureCounts.rds')  
	# 		fc <- load(file='featureCounts.RData')
		
	cat('\nFEATURE COUNT FINISHED\n')
	
	return(fc)
}


#' merge.count.files
#'
#' Takes an array of files containing individual count data and merge then by
#' the gene_id(leftmost column). Input files should be in tabular format and
#' each of them should correspond to a different sample. The merged output
#' table is saved in the same directory where the individual counts are found.
#' 
#' @param	count.files	Vector with the path of the individual feature count files
#' 
#' @return	table with merged count files
#'
#' @usage	merge.count.files(count.files)
#' 
#' @examples	merge.count.files(c("sample1.cnt", "sample2.cnt", "sample3.cnt"))
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
merge.count.files <- function(count.files) {

	if (length(count.files) <= 0) cat.err("Feature Count Files not found.\n", abort = TRUE)    
	
	labels <- sub('\\..*$', '', basename(count.files))
	
	cat('\nMerging count files:\n\n')
	cat(count.files, sep='\n')

	# Each file is saved as an element of the list. Each element is a table with two coulmns:
	#	Gene_ID and "Sample name" (counts).

	count.list <- list()
	for (i in 1:length(labels)) {
		label <- labels[i]
		filename <- count.files[i]
		count.list[[label]] <- read.table(filename, row.names = 'V1', header = F, stringsAsFactors = F)  #row.names='V1'
		colnames(count.list[[label]]) <- label
	}

	# Bind all htseq count files into a unique data frame by GeneID.
	all.counts <- bind_cols(count.list)
	
	# Save the merged dataframe into a file
	return(all.counts)
}

#' h.cluster.changes
#'
#' Cluster gene expression data using hierarchycal clustering
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
h.cluster.changes <- function(data.table, annotated.data, 
                              distance="euclidean",	# see ?dist()
                              clusters=1,
                              method="complete",	# see ?hclust()
                              estimate=TRUE,
                              interactive=TRUE
                              ) {
                              

    # Normalize the data
    by.row <- 1
    by.col <- 2
    dif <- data.table
    means <- apply(dif, by.col, mean)
    sds <- apply(dif, by.col, sd)		# standard deviations
    nor <- scale(dif,center=means,scale=sds)	# normalized values

    # Calculate distance matrix  
    distan = dist(nor, method=distance)

    # Hierarchical agglomerative clustering  
    cat.info("
    H I E R A R C H I C A L   C L U S T E R I N G
    =============================================
    \n\n")
    
    hc = hclust(distan)
#    if (verbose) {
#        as.png(plot(hc), 'hclust.png')
#        as.png(plot(hc,hang=-1), 'hclust.hang.png')
#        as.png(plot(hc,labels=rownames(data.table),main='Default from hclust')
#               "hclust.labelled.png")
#    }
    if (interactive) {
        print(plot(hc,labels=rownames(data.table),main='Default from hclust'))
        continue.on.enter("Press [ENTER] to continue ")
    } 
       
    # Cluster membership
    if ((clusters <= 1) & (estimate == T)) {
        ans <- continue.on.enter(prompt="How many clusters should we use? ")
        nclust <- as.numeric(ans)
    } else
        nclust <- clusters
    
    member <- cutree(hc, nclust)	# cut to 4 groups
    cat.info("Cluster membership counts\n")
    print(table(member))
    # Characterizing clusters 
    cat('\n')
    cat.info("Means by cluster in the normalized data\n")
    print(aggregate(nor, list(member), mean))
    cat('\n')
    cat.info("Means by cluster in the non-normalized data\n")
    print(aggregate(dif, list(member), mean))

    for (i in 1:nclust) {
        clus.i <- member[ member == i ]
        data.i <- annotated.data[annotated.data$ensembl.gene.id %in% names(clus.i), ]
        # sort by l2FC
        data.i <- data.i[order(data.i$log2FoldChange, decreasing=T), ]
        cat.info("Showing annotation for cluster", i, "(", length(clus.i), " elements)\n")
        
        #View(data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")],
        #     title=paste("Cluster no.", i, "(", length(clus.i), ") elements)"))
        
        # or, using tcltk and gWidget2
        #library(tcltk)
        #library(gWidgets2)
        data.to.show <- data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")]
        # clean up for showing
        data.to.show[ is.na(data.to.show) ] <- 'NA'
        window.name <- paste("Cluster no.", i, "(", length(clus.i), " elements)")
        #window <- gwindow(title=window.name, visible=TRUE)
        #tab <- gtable(data.to.show,
        #       container=window)
        window <- show.data.frame(data.to.show, window.name, visible=F)
        if (interactive) {
          visible(window) <- TRUE	# setting it to FALSE removes window
          keypress()
          visible(window) <- FALSE
        }
    }

    cat.info("
    
    Silhouette Plot for hierarchical clustering with normalized data
    ----------------------------------------------------------------
    Measure similarity of each object to its own cluster (cohesion) compared to
    other clusters (dispersion). Values range from -1 to +1. Large values
    indicate objects well matched to their own cluster and badly to neighboring
    clusters. If many points have low or negative value, the number of clusters
    is too low or too high.
    \n")
    #if (verbose)
    #    as.png(plot(silhouette(cutree(hc, nclust), distan)), "silhouette.png")
    print(plot(silhouette(cutree(hc, nclust), distan)))

    return(hc)
}



#' hcut.cluster changes
#'
#' Cluster gene expression data using hierarchycal clustering with 'hcut'
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
hcut.cluster.changes <- function(data.table, annotated.data,
                                    algorithm="hclust",	# see ?hcut()
                                    distance="euclidean", # see ?hcut()
                                    method="ward.D2",	# see ?hcut()
                                    clusters=1,	# n. of clusters
                                    nstart=1,	# ignored
                                    gap_bootstrap=500,
                                    estimate=FALSE
                                    ) {

    # Normalize the data
    by.row <- 1
    by.col <- 2
    dif <- data.table
    means <- apply(dif, by.col, mean)
    sds <- apply(dif, by.col, sd)
    nor <- scale(dif,center=means,scale=sds)

    cat.info("
    H c u t
    -------
    \n")

    if (estimate) {
        cat("Suggestions for the best number of clusters\n")

        cat("Plot by within-cluster sums of squares\n")
        cat("    Elbow method: look for a knee\n")
        print(fviz_nbclust(nor, hcut, method="wss"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Average Silhouette Method

        The average silhouette approach measures the quality of a clustering. It
        determines how well each observation lies within its cluster.

        A high average silhouette width indicates a good clustering. The average
        silhouette method computes the average silhouette of observations for
        different values of k.
        \n")
        print(fviz_nbclust(nor, hcut, method="silhouette"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Gap Statistic Method

        This approach can be utilized in any type of clustering method (i.e.
        K-means clustering, hierarchical clustering).

        The gap statistic compares the total intracluster variation for different
        values of k with their expected values under null reference distribution of
        the data.

        \n")
        gap_stat <- clusGap(nor, FUN=hcut, nstart=25,
                            K.max=15, B=gap_bootstrap)
        print(fviz_gap_stat(gap_stat))
        continue.on.enter("Press [ENTER] to continue ")
    }

    # Number of clusters
    if (clusters == 1) {
        ans <- continue.on.enter(prompt="How many clusters should we use? ")
        nclust <- as.numeric(ans)
    } else
        nclust <- clusters
    
    # clustering (N groups)
    set.seed(123)
    hc <- hcut(nor, k=nclust, 
               hc_func=algorithm, hc_method=method, hc_metric=distance, is_diss=FALSE)
    print(head(hc))
    member <- hc$cluster
    
    cat.info("Cluster membership counts\n")
    print(table(member))
    # Characterizing clusters 
    cat('\n')
    cat.info("Means by cluster in the normalized data\n")
    print(aggregate(nor, list(member), mean))
    cat('\n')
    cat.info("\nMeans by cluster in the non-normalized data\n")
    print(aggregate(dif, list(member), mean))

    for (i in 1:nclust) {
        clus.i <- member[ member == i ]
        data.i <- annotated.data[annotated.data$ensembl.gene.id %in% names(clus.i), ]
        # sort by l2FC
        data.i <- data.i[order(data.i$log2FoldChange, decreasing=T), ]
        cat("Showing annotation for cluster", i, "(", length(clus.i), ") elements)\n")

        # invoke a spreadsheet-style data viewer
        View(data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")],
             title=paste("Cluster no.", i, "(", length(clus.i), ") elements)"))
 #       keypress()
    }
    
    # plot a PCA of the clusters
    # If there are more than two dimensions (variables) fviz_cluster will perform
    # principal component analysis (PCA) and plot the data points according to
    # the first two principal components that explain the majority of the
    # variance.
    print(fviz_cluster(hc, data=nor, geom="point"))
    continue.on.enter("Press [ENTER] to continue ")
    # to plot gene names instead of ensembl.ids we need to change rownames to genes
    gor <- nor 
    rownames(gor) <- paste(annotated.data[ , "GENENAME"], 1:length(rownames(gor)))
    print(fviz_cluster(hc, data=gor))
    continue.on.enter("Press [ENTER] to continue ")
}



#' k.means.cluster changes
#'
#' Cluster gene expression data using K-means clustering
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
k.means.cluster.changes <- function(data.table, annotated.data,
                                    algorithm="Hartigan-Wong",	# see kmeans()
                                    clusters=1,	# n. of clusters
                                    nstart=1,	# n. of random start sets to choose
                                    estimate=F,
                                    gap_bootstrap=500
                                    ) {

    # Normalize the data
    by.row <- 1
    by.col <- 2
    dif <- data.table
    means <- apply(dif, by.col, mean)
    sds <- apply(dif, by.col, sd)
    nor <- scale(dif,center=means,scale=sds)

    cat.info("
    K - m e a n s
    -------------
    \n")
    
    if (estimate) {
        cat("Suggestions for the best number of clusters\n")

        #
        # Scree Plot
        cat("

        Scree plot using normalized data

        This allows us to evaluate how much variation we account for as we
        consider more clusters and decide what a reasonable number of clusters
        might be.
        We draw here the variances accounted for using up to 20 K-means clusters
        \n") 
        # compute variances by row
        wss <- (nrow(nor)-1)*sum(apply(nor, by.row, var))
        for (i in 2:20) wss[i] <- sum(kmeans(nor, centers=i)$withinss)
        print(plot(1:20, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares") )
        continue.on.enter("Press [ENTER] to continue ")

        # The scree plot will allow us to see the variabilities in clusters, 
        # we expect that if we increase the number of clusters, then the 
        # within-group sum of squares would come down. 
        #

        cat("Plot by within-cluster sums of squares\n")
        cat("    Elbow method: look for a knee\n")
        print(fviz_nbclust(nor, kmeans, method="wss"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Average Silhouette Method

        The average silhouette approach measures the quality of a clustering. It
        determines how well each observation lies within its cluster.

        A high average silhouette width indicates a good clustering. The average
        silhouette method computes the average silhouette of observations for
        different values of k.
        \n")
        print(fviz_nbclust(nor, kmeans, method="silhouette"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Gap Statistic Method

        This approach can be utilized in any type of clustering method (i.e.
        K-means clustering, hierarchical clustering).

        The gap statistic compares the total intracluster variation for different
        values of k with their expected values under null reference distribution of
        the data.

        \n")
        gap_stat <- clusGap(nor, FUN=kmeans, nstart=25,
                            K.max=15, B=gap_bootstrap)
        print(fviz_gap_stat(gap_stat))
        continue.on.enter("Press [ENTER] to continue ")
    }

    # Number of clusters
    if (clusters == 1) {
        ans <- continue.on.enter(prompt="How many clusters should we use? ")
        nclust <- as.numeric(ans)
    } else
        nclust <- clusters
    
    # K-means clustering (N groups)
    kc <- kmeans(nor, centers=nclust, nstart=nstart, algorithm=algorithm)
    print(head(kc))
    member <- kc$cluster
    
    cat("Cluster membership counts\n")
    print(table(member))
    # Characterizing clusters 
    cat("\nMeans by cluster in the normalized data\n")
    print(aggregate(nor, list(member), mean))
    cat("\nMeans by cluster in the non-normalized data\n")
    print(aggregate(dif, list(member), mean))

    for (i in (1:nclust)) {
        clus.i <- member[ member == i ]
        data.i <- annotated.data[annotated.data$ensembl.gene.id %in% names(clus.i), ]
        # sort by l2FC
        data.i <- data.i[order(data.i$log2FoldChange, decreasing=T), ]
        cat("Showing annotation for cluster", i, "(", length(clus.i), ") elements)\n")
        View(data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")],
             title=paste("Cluster no.", i, "(", length(clus.i), ") elements)"))
 #       keypress()
    }
    
    # plot a PCA of the clusters
    # If there are more than two dimensions (variables) fviz_cluster will perform
    # principal component analysis (PCA) and plot the data points according to
    # the first two principal components that explain the majority of the
    # variance.
    print(fviz_cluster(kc, data=nor, geom="point"))
    continue.on.enter("Press [ENTER] to continue ")
    # to plot gene names instead of ensembl.ids we need to change rownames to genes
    gor <- nor 
    rownames(gor) <- paste(annotated.data[ , "GENENAME"], 1:length(rownames(gor)))
    print(fviz_cluster(kc, data=gor))
    continue.on.enter("Press [ENTER] to continue ")
}



#'
#' pam.cluster.changes
#'
#'  Cluster gene expression data using Partition around medoids clustering
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
pam.cluster.changes <- function(data.table, annotated.data,
                                    distance="euclidean", 	# see dist()
                                    metric="euclidean",		# see pam()
                                    clusters=1,		# n. of clusters
                                    nstart=1,	# n. of random start sets to choose
                                    estimate=T,
                                    gap_bootstrap=500
                                    ) {
                                    
    if (nstart != 1) medoids="random" else medoids=NULL
        
    # Normalize the data
    by.row <- 1
    by.col <- 2
    dif <- data.table
    means <- apply(dif, by.col, mean)
    sds <- apply(dif, by.col, sd)
    nor <- scale(dif,center=means,scale=sds)

    cat("
    Partition around medoids
    ------------------------
    \n")

    if (estimate) {
        cat("Suggestions for the best number of clusters\n")

        cat("    Plot by within-cluster sums of squares\n\n")
        cat("    Elbow method: look at the knee\n")
        print(fviz_nbclust(nor, pam, method="wss"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Average Silhouette Method

        The average silhouette approach measures the quality of a clustering. It
        determines how well each observation lies within its cluster.

        A high average silhouette width indicates a good clustering. The average
        silhouette method computes the average silhouette of observations for
        different values of k.
        \n")
        print(fviz_nbclust(nor, pam, method="silhouette"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Gap Statistic Method

        This approach can be utilized in any type of clustering method (i.e.
        K-means clustering, hierarchical clustering).

        The gap statistic compares the total intracluster variation for different
        values of k with their expected values under null reference distribution of
        the data.

        \n")
        gap_stat <- clusGap(nor, FUN=pam, nstart=25,
                            K.max=15, B=gap_bootstrap)
        print(fviz_gap_stat(gap_stat))
        continue.on.enter("Press [ENTER] to continue ")
    }


    # Number of clusters
    if (clusters == 1) {
        ans <- continue.on.enter(prompt="How many clusters should we use? ")
        nclust <- as.numeric(ans)
    } else
        nclust <- clusters

    # cluster with partition around medioids for k=N clusters 
    # (data is not dissimilarity but distance)
    # compute distance (we have dissimilarity to a common reference but not
    # between the considered groups).
    cat("Computing clusters (may take some time)...\n")
    eucldist <- dist(nor, method=distance) 
    #cluster
    ### NOTE ### NOTE: may be worth trying to cluster separately with diss=TRUE)
    pam.clus <- pam(eucldist, k=nclust, diss=TRUE)
    print(head(pam.clus))
    member <- pam.clus$clustering

    cat("Cluster membership information\n")
    print(pam.clus$clusinfo)
    cat("\nMeans by cluster in the normalized data\n")
    print(aggregate(nor, list(member), mean))
    cat("\nMeans by cluster in the non-normalized data\n")
    print(aggregate(dif, list(member), mean))
    
    for (i in (1:nclust)) {
        clus.i <- member[ member == i ]
        data.i <- annotated.data[annotated.data$ensembl.gene.id %in% names(clus.i), ]
        # sort by l2FC
        data.i <- data.i[order(data.i$log2FoldChange, decreasing=T), ]
        cat("Showing annotation for cluster", i, "(", length(clus.i), ") elements)\n")
        View(data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")],
             title=paste("Cluster no.", i, "(", length(clus.i), ") elements)"))
 #       keypress()
    }
    
    # plot the partitioning
    print(clusplot(pam.clus, shade = FALSE,labels=F,
	    col.clus="blue",col.p="red",
            span=FALSE,
            main="PAM Cluster Mapping",cex=1.2))
    continue.on.enter("Press [ENTER] to continue ")

    # plot a PCA of the clusters
    # If there are more than two dimensions (variables) fviz_cluster will perform
    # principal component analysis (PCA) and plot the data points according to
    # the first two principal components that explain the majority of the
    # variance.
    print(fviz_cluster(pam.clus, data=nor, geom="point"))
    continue.on.enter("Press [ENTER] to continue ")
    # to plot gene names instead of ensembl.ids we need to change rownames to genes
    gor <- nor 
    rownames(gor) <- paste(annotated.data[ , "GENENAME"], 1:length(rownames(gor)))
    # this seemingly ignores the data argument!!!
    print(fviz_cluster(pam.clus, data=gor))
    continue.on.enter("Press [ENTER] to continue ")
}


#' cluster.changes
#'
#' Cluster gene expression data using any of a variety of methods for
#' clustering
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
cluster.changes <- function(data.table, annotated.data,
                            FUN=hcut,
                            clusters=1,			# n. of clusters
                            algorithm="default",	# see below	
                            distance="default", 	# see below
                            method="default",		# see see below
                            nstart=1,		# n. of ran dom start sets to choose
                            eps=1.0,		# epsilon for DBScan
                            gap_bootstrap=100,
                            normalize=TRUE,
                            estimate=TRUE,	# only if clusters > 1
                            output.folder=NULL	# NULL => no output desired
                                    ) {
# FUN -- one of c(hcut, kmeans, pam)
# algorithm -- for 'hcut' one of c(_"hclust"_, "agnes", "diana")
#              for 'kmeans' one of c(_"Hartigan-Wong"_, "Lloyd", "Forgy", "MacQueen")
#              for 'pam' one of c(_"original"_, "o_1", "o_2", "f_3", "f_4", "f_5", 
#                       "faster")
# distance -- for 'hcut' one of c(_"euclidean"_, "manhattan", "maximum", 
#                       "canberra", "binary", "minkowski", "pearson", "spearman", 
#                       "kendall")
#             for 'kmeans' it is ignored
#             for 'pam' one of c(_"euclidean"_, manhattan")
#             for 'dbscan' it is ignored
# method -- for 'hcut' one of c("ward.D"', "ward.D2", "single", "complete", 
#                       "average" (= UPGMA), "mcquitty" (= WPGMA), 
#                       "median" (= WPGMC), "centroid" (= UPGMC), 
#                       "weighted" (=WPGMA), "flexible", "gaverage")
#             for 'kmeans' ignored
#             for 'pam' ignored
#             for 'dbscan' ignored

    # Normalize the data
    by.row <- 1
    by.col <- 2
    if (normalize == TRUE) {
        means <- apply(data.table, by.col, mean)
        sds <- apply(data.table, by.col, sd)
        nor <- scale(data.table,center=means,scale=sds)
    } else {
        nor <- data.table
    }

    fun.name <- deparse(substitute(FUN))
    cat("
    C L U S T E R I N G    W I T H :   ", fun.name, "
    ---------------------------------------------
    \n")

    # re-activate next line to save computation time
    #if (clusters > 1) estimate <- FALSE		# we already know the number
    if (clusters < 1) estimate <- TRUE
    # we may have clusters == 1 and estimate == FALSE, e.g. to find outliers
    
    if (estimate == TRUE) {
        cat("Suggestions for the best number of clusters using", fun.name, "\n")

        if (normalize == TRUE) {
            cat("Plot by within-cluster sums of squares\n")
            cat("    Elbow method: look for a knee (normalized data)\n")
            print(fviz_nbclust(nor, FUN, method="wss"))
            out.png <- sprintf("%s/DESeq2_%s_wss.png",
                output.folder, fun.name)
            as.png(fviz_nbclust(nor, FUN, method="wss"), out.png)
            continue.on.enter("Press [ENTER] to continue ")
        }

        cat("Plot by within-cluster sums of squares\n")
        cat("    Elbow method: look for a knee (raw data)\n")
        print(fviz_nbclust(data.table, FUN, method="wss"))
        out.png <- sprintf("%s/DESeq2_%s_wss.png",
            output.folder, fun.name)
        as.png(fviz_nbclust(data.table, FUN, method="wss"), out.png)
        continue.on.enter("Press [ENTER] to continue ")

        if (fun.name != "dbscan") {
            cat("
            Average Silhouette Method

            The average silhouette approach measures the quality of a clustering. It
            determines how well each observation lies within its cluster.

            A high average silhouette width indicates a good clustering. The average
            silhouette method computes the average silhouette of observations for
            different values of k.
            \n")
            print(fviz_nbclust(nor, FUN, method="silhouette"))
            out.png <- sprintf("%s/DESeq2_%s_silhouette.png",
                output.folder, fun.name)
            as.png(fviz_nbclust(nor, FUN, method="silhouette"), out.png)

            continue.on.enter("Press [ENTER] to continue ")
        }

        cat("
        Gap Statistic Method

        This approach can be utilized in any type of clustering method (e.g.
        K-means, hierarchical clustering, partition around medoids).

        The gap statistic compares the total intracluster variation for different
        values of k with their expected values under null reference distribution of
        the data.

        We are looking for the highest peak identified.

        \n")
        gap_stat <- clusGap(nor, FUN=FUN, K.max=15, B=gap_bootstrap)
        print(fviz_gap_stat(gap_stat))
        out.png <- sprintf("%s/DESeq2_%s_gap_stat.png",
            output.folder, fun.name)
        as.png(fviz_gap_stat(gap_stat), out.png)
        continue.on.enter("Press [ENTER] to continue ")
    }
    
    # Select the number of clusters
    if ((clusters < 1) & (fun.name != 'dbscan')) {
        ans <- continue.on.enter(prompt="How many clusters should we use? ")
        nclust <- as.numeric(ans)
    } else {
        nclust <- clusters
    }
    
    # Do the clustering (N groups). This is necessarily method-specific
    if (fun.name == "hcut") {
        if (algorithm == "default") algorithm <- 'hclust'
        if (distance == "default") distance <- "euclidean"
        if (method == "default") method <- "complete"
        
        cl <- hcut(nor, k=nclust, 
               hc_func=algorithm, hc_method=method, hc_metric=distance, is_diss=FALSE)
        print(head(cl))
        member <- cl$cluster
    } else if (fun.name == "kmeans") {
        if (algorithm == "default") algorithm <- "Hartigan-Wong"
        if (distance == "default") distance <- "euclidean"
        
        cl <- kmeans(nor, centers=nclust, nstart=nstart, algorithm=algorithm)
        print(head(cl))
        member <- cl$cluster
    } else if (fun.name == "pam") {
        if (algorithm == "default") algorithm <- "faster"
        if (distance == "default") distance <- "euclidean"
        
        # cluster with partition around medioids for k=N clusters 
        # (data is not dissimilarity but distance)
        # compute distance (we have dissimilarity to a common reference but not
        # between the considered groups).
        cat("Computing clusters (may take some time)...\n")
        use_dist_matrix <- FALSE	### NOTE ### fviz_cluster fails, why?
        if (use_dist_matrix == TRUE) {
            ### NOTE ### It may be worth trying to cluster with diss=TRUE)
            # calculate dissimilarity matrix and cluster it
            #distm <- get_dist(nor, method=distance, stand=TRUE) 
            cl <- pam(get_dist(nor, method=distance, stand=TRUE),
                      k=nclust, diss=TRUE, 
                      nstart=nstart, metric=distance, variant=algorithm)
        } else {
            cl <- pam(nor, k=nclust, diss=FALSE, 
                      nstart=nstart, metric=distance, variant=algorithm)
        }

        print(head(cl))
        member <- cl$cluster
        cat("Summary\n")
        print(cl$clusinfo)
    }  else if (fun.name == "dbscan") {
        if (algorithm == "defaut") algorithm <- "hybrid"
        if (distance == "default") distance <- "manhattan"
        cl <- dbscan(nor, eps=eps, MinPts=4, showplot=1)
        # showplot=1 makes it produce a movie plot
        continue.on.enter("Press [ENTER] to continue ")
        print(head(cl))
        member <- cl$cluster
        names(member) <- annotated.data$ensembl.gene.id
        nclust <- max(member)
        plot(cl, nor, main="DBScan")
        out.png <- sprintf("%s/DESeq2_%s_plot.png",
            output.folder, fun.name)
        as.png(plot(cl, nor, main="DBScan"), out.png)
        continue.on.enter("Press [ENTER] to continue ")
    }
        
    # this may take very long and is of little use for now
    if (FALSE) {
        cat("Distance plot\n")
        # plot distances
        distm <- get_dist(nor, distance, stand=TRUE)
        print(fviz_dist(distm,
                  gradient=list(low="blue", mid="white", high="red")))
        out.png <- sprintf("%s/DESeq2_%s_distances.png",
            output.folder, fun.name)
         as.png(fviz_dist(distm,
                    gradient=list(low="blue", mid="white", high="red")),
                    out.png)
        continue.on.enter("Press [ENTER] to continue ")
    }
    
    cat("Cluster membership counts\n")
    print(table(member))
    # Characterizing clusters 
    cat("\nMeans by cluster in the normalized data\n")
    print(aggregate(nor, list(member), mean))
    cat("\nMeans by cluster in the non-normalized data\n")
    print(aggregate(data.table, list(member), mean))

    for (i in (min(member):max(member))) {
        clus.i <- member[ member == i ]
        data.i <- annotated.data[annotated.data$ensembl.gene.id 
                                 %in% names(clus.i), ]
        # sort by l2FC
        data.i <- data.i[order(data.i$log2FoldChange, decreasing=T), ]
        cat("Showing annotation for cluster", i, "(", length(clus.i), ") elements)\n")
        #View(data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")],
        #     title=paste("Cluster no.", i, "(", length(clus.i), ") elements)"))
        d <- data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")]
        d[ is.na(d) ] <- 'NA'
        t <- paste("Cluster no.", i, "(", length(clus.i), ") elements")
        w <- show.data.frame(d, t, TRUE)
        visible(w) <- FALSE
        #visible(w) <- TRUE
#        keypress()
         clus.file <- sprintf("%s/DESeq2_%s_nc=%03d_c=%03d.tab",
             output.folder, fun.name, max(member), i)
         write.table(data.i, file=clus.file,
             row.names=T, col.names=T, sep='\t')
    }
    
    # plot a PCA of the clusters
    # If there are more than two dimensions (variables) fviz_cluster will perform
    # principal component analysis (PCA) and plot the data points according to
    # the first two principal components that explain the majority of the
    # variance.
    print(fviz_cluster(cl, data=nor, geom="point", ellipse.type="convex"))
    out.png <- sprintf("%s/DESeq2_%s_nc=%03d_PCA.png",
            output.folder, fun.name, max(member))
    as.png(fviz_cluster(cl, data=nor, geom="point", ellipse.type="convex"),
           out.png)
    continue.on.enter("Press [ENTER] to continue ")
    # to plot gene names instead of ensembl.ids we need to change rownames to genes
    gor <- nor 
    rownames(gor) <- paste(annotated.data[ , "GENENAME"], 1:length(rownames(gor)))
    rcl <- cl
    if (fun.name == 'hcut') {
        cat("Setting names to gene names\n")
        dimnames(rcl$data)[[1]] <- rownames(gor)
    } else if (fun.name == 'pam') {
        names(cl$cluster <- rownames(gor))
    }
    print(fviz_cluster(rcl, data=gor, 
              ellipse.type="convex", show.clust.cen=FALSE, labelsize=6))
    out.png <- sprintf("%s/DESeq2_%s_nc=%03d_PCA_genes.png",
            output.folder, fun.name, max(member))
    as.png(fviz_cluster(rcl, data=gor, 
              ellipse.type="convex", show.clust.cen=FALSE, labelsize=6),
           out.png)
    continue.on.enter("Press [ENTER] to continue ")
    
    return(cl)
}



# # # # edgeR functions

#
# SAVE TOP 'N' GENES
# ------------------
# save the table of the n.genes top expressed genes sorted in various orders
eR.save.top <- function(fit, folder, n.genes=500, sort.by='p', coef=1, p.value=0.01)
{
    for (by in sort.by) {
        # default adjustment is BH
        tt <- topTable(fit, 
                       coef=coef, 
                       number=n.genes, 
                       sort.by=by, 
                       p.value=p.value)
        
		file <- paste(folder, 
                      '/edgeR/cmp_coef=', coef, 
                      '_top_', n.genes, 
                      '_by', by, 
                      '.tab', 
                      sep='')
        write.table(tt, file, row.names=T, col.names=T, sep='\t')
    }
}


#
# SAVE TOP 'N' GENES ANNOTATED
# ----------------------------
# we can now run again the topTable and we will have all the annotation 
# information linked 
#	n.genes <- 500 ALREADY DEFINED ABOVE
eR.save.top.annotated <- function(fit, folder, n.genes=500, sort.by='p', coef=1, p.value=0.01) {

    for (by in sort.by) {
        tt <- topTable(fit, 
                       coef=coef, 
                       number=n.genes, 
                       sort.by=by, 
                       p.value=p.value)
        file <- paste(folder, 
                      '/edgeR/cmp=', coef, 
                      '_top_', n.genes, 
                      '_by_', by, 
                      '_annotated.tab', 
                      sep='')
        write.table(tt, file, row.names=T, col.names=T, sep='\t')
    }
}



# compare to a threshold and save
# SAVE TOP RESULTS RELATIVE TO THRESHOLD
# --------------------------------------
eR.save.top.treat.ann <- function(fit, folder, n.genes=500, sort.by='p', coef=1, p.value=0.01, threshold=1) {
    for (by in sort.by) {
        #     vvvvvvvv		here we use topTreat instead of topTable
        tt <- topTreat(fit, 
                       coef=coef, 
                       number=n.genes, 
                       sort.by=by, 
                       p.value=p.value)
        file <- paste(folder, 
                      '/edgeR/cmp=', coef, 
                      '_top_', n.genes, 
                      '_by_', by, 
                      '_lfc>=', threshold, 
                      '_annotated.tab', 
                      sep='')
        write.table(tt, file, row.names=T, col.names=T, sep='\t')
    }
}



# save top N genes, annotated with biomart
eR.save.top.ann.thresh <- function(fit, folder, n.genes=500, sort.by='p', coef=1, p.value=0.01, threshold=1) {
    for (by in sort.by) {
        file <- paste(folder, 
                      '/edgeR/cmp=', coef, 
                      '_top_', n.genes, 
                      '_by_', by, 
                      '_lfc>=', threshold, 
                      '_threshold_annotated.tab', sep='')
        tt <- topTable(fit, 
                       coef=coef, 
                       number=n.genes, 
                       sort.by=by, 
                       p.value=p.value)
        write.table(tt, file, row.names=T, col.names=T, sep='\t')
    }
}



# annotation functions for edgeR

eR.fit.annotate <- function(fit, ens.db, biomart.db) {

	if (all(str_sub(rownames(fit), 1, 3) == "ENS")){
        by='GENEID'
    }else{
    	by='GENENAME'}
        
     ens.ann <- ensembldb::select(ens.db, 
                  keytype=by, keys=rownames(fit), 
                  columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION',
                             'GENENAME', 'GENEID', 'ENTREZID',
                             'TXNAME', 'TXID', 'TXBIOTYPE',
                             'PROTEINID', 'UNIPROTID'
                            ))

#    ann <- ensembldb::select(ens.db, 
#              keytype= 'GENEID', keys=rownames(fit), 
#              columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION',
#                         'GENENAME', 'GENEID', 'ENTREZID', 
#                          'TXNAME', 'TXID', 'TXBIOTYPE', 
#                          'PROTEINID', 'UNIPROTID'))

    # use only first annotation
    ens.ann.1 <- ens.ann[ ! duplicated(ens.ann$GENEID), ]

    # As long as it works perfectly we continue by connecting the 
    # annotation information to the fit using a new bucket named "genes"
    # (fit is an object of class "MArrayLM" (package "limma"), DGELRT or
    # DGELM (package edgeR) which is a list. I.e. we assign the annotation 
    # to a new element named 'genes' of this list.
    #	This works if they are in the same order
    #fit$genes <- ens.ann.1
    #	This works by matching rownames and gene IDs
    fit$genes <- ens.ann.1[ match(rownames(fit), ens.ann.1[,by], nomatch = NA), ]
    rownames(fit$genes) <- rownames(fit)
    #
    # add biomart annotation
    fit$genes <- merge(fit$genes, biomart.db, by.x="GENEID", by.y="ensembl_gene_id", all.x = TRUE)
    #head(fit$genes)


    return (fit)

}

eR.save.fit <- function(fit, name) {
    # SAVE THE ANNOTATED FIT
    # ----------------------
    # save all the contents of 'fit' in an RDS file
    saveRDS(fit, file=paste(name,'.rds', sep=''))
    #	'fit' can later be recovered with: 
    #		fit <- readRDS(file=paste(folder, 'edgeR/', name, '.rds', sep=''))
    # and save as well as Rdata file
    save(fit, file=paste(name, '.RData', sep=''))
    #	'fit' can later be recovered with: 
    #		fc <- load(file=paste(folder, '/edgeR/', name, '.RData', sep=''))
}

eR.save.top.fit <- function(fit, file, n.genes=500, sort.by='PValue', p.value=0.01) {    # default adjust.method is BH
    # default p.value is 1 (all genes)
    tt <- topTags(fit, n=n.genes, sort.by=sort.by, p.value=p.value)
    n <- dim(tt)[1]
    # cap at the maximum number of genes requested
    if (n > n.genes) n <- n.genes
    file.name <- paste(file, '_top_', n, '_by_', sort.by, '.tab', sep='')
    write.table(tt$table[1:n, ], file.name, 
	row.names=T, col.names=T, sep='\t')
    # if fit is annotated the annotation will also be saved
}

#' Plot counts-counts per million (cpm) correlation in order to define a count
#' threshold.
#'
counts.cpm.plot <- function(counts, cpm, out.png = NULL) {

	as.png({	plot(x = as.matrix(counts), y = as.matrix(cpm), xlim = c(0,20), ylim = c(0,5),
						pch = 1, cex = 0.8, xlab = 'Counts', ylab = 'Counts per Million (CPM)')
				abline(v = 10, col = 2, lwd = 2)
				abline(v = 15, col = 2, lwd = 2)
				#abline(h = 0.3, col = 2)
				arrows(10, 3, 15, 3, angle = 20, code = 3, col = 1, length = 0.2, lwd = 3, lty = 1)
				text(x=12.5, y=3.2, 'Expression\nthreshold', cex= 1)
	}, out.png, overwrite = TRUE)
}



eR.differential.gene.expression<- function(counts, 
                                           metadata,
										   design.column,
                                           cpm.threshold,
										   n.genes,
                                           ens.db,
                                           folder)
{

    # then we proceed to the analysis.. for this, we will need other packages
    #library(edgeR)
    #library(limma)
    #library(RColorBrewer)
    #library(gplots)

    # The next step is to distinguish the genes whose expression is significant
    # from  the ones that have an 'insignificant' expression that could be by
    # chance or irrelevant to the  situation of the cell. So we set up the
    # threshold of expression for each gene to a minimum 10-15 counts. Since
    # we will work with normalized CPM (counts per million) data, we need
    # to check which number of CPM corresponds to  to this amount of counts in
    # order to filter the data

    # convert counts to CPM
    cpm.counts <- cpm(counts)

    # we'll plot the correlation of cpm and counts to see which is the number of
    # cpm that corresponds to 10-15 counts minimum. we see this information
    # graphically we use for example column 1... we could check every file
    # independently but since they are all similar in number of reads there 
    # should be no need

    out.png <- paste(folder, '/edgeR/img/edgeR_CPMcounts.png', sep='')
    counts.cpm.plot(counts = counts, cpm = cpm.counts, out.png = out.png)

    # we will set up the threshold to 0.25 according on the plot we have drawn
    # this command will return a table of trues and falses, then, we want to keep
    # only the rows that exceed the threshold at least in three different
    # cases(experiments .bam files)

    cat(paste("\n\tFiltering genes by CPM threshold", cpm.threshold, "\n\n"))
    thres <- cpm.counts > cpm.threshold		# 0.5 for Gallus
    keep <- rowSums(thres) >=3
    counts.keep <- counts[keep,]

	if(VERBOSE) print(table(keep))

    # then we store in a different variable the genes whose counts exceed the
    # threshold and visualise the content of the new variable to see the amount of
    # remaining genes
   	
	if ( VERBOSE ) {
        # some paranoid manual checks
        # we are using counts instead of cpm
        cpmavg <- data.frame(vd000.0=apply(counts.keep[,13:15],1,mean), 
                             vd000.1=apply(counts.keep[,1:3],1,mean),
                             vd001.0=apply(counts.keep[,10:12],1,mean), 
                             vd010.0=apply(counts.keep[,7:9],1,mean), 
                             vd100.0=apply(counts.keep[,4:6],1,mean) 
                             )

        f.c000.0 <- data.frame(vd000.0_vd000.1=(cpmavg$vd000.0 / cpmavg$vd000.1),
                               vd000.0_vd001.0=(cpmavg$vd000.0 / cpmavg$vd001.0),
                               vd000.0_vd010.0=(cpmavg$vd000.0 / cpmavg$vd010.0),
                               vd000.0_vd100.0=(cpmavg$vd000.0 / cpmavg$vd100.0)
                              ) 
        row.names(f.c000.0) <- row.names(cpmavg)
        l.f.c000.0 <- log2(f.c000.0)
        write.table(f.c000.0, file=paste(folder, '/edgeR/hand.fc_000.0.tab', sep=''), sep='\t')
        write.table(l.f.c000.0, file=paste(folder, '/edgeR/hand.lfc_000.0.tab', sep=''), sep='\t')
        write.table(cpmavg, file=paste(folder, '/edgeR/hand.cpmavg.tab', sep=''), sep='\t')
    }

    # We have manipulated the data discarding whatever is not of high interest
    # and now we need to see the differencial expression and highlight the
    # differences among the cells

    # Convert the counts.keep to a DGEList
    dge <- DGEList(counts.keep)
	dge$samples$group <- as.factor(metadata[ ,design.column])
    
	# do TMM normalization
    cat("\n\tCalculating TMM normalization factors\n\n")
    dge <- calcNormFactors(dge)
    if (VERBOSE) print(dge$samples)

    # Plot the library size of the different samples.
    out.png <- paste(folder, '/edgeR/img/edgeR_sample_lib_size.png', sep='')
    as.png(barplot(dge$samples$lib.size, cex.names= 0.8,
					main = "Library Size",
					col = dge$samples$group,
					names.arg=dge$samples$group,
					ylab = "Reads"), out.png, overwrite=TRUE)


    # Now, do some quality control plots, barplots and boxplots
    # we need normalized data counts so we take the logarithm
    # we need to group together all the experiment data that correspond to the
    # same  sample (target) and asign a different color to each one so as to
    # distinguish them graphically.
	
    # We first set up the colours we will use for the different plots 
    # and then we create all the colors in between in the palette
    mypalette <- brewer.pal(11, 'RdYlBu')
    morecolors <- colorRampPalette(mypalette)
	
	logcpm <- cpm(dge$counts, log=TRUE)
	colors <- morecolors(length(levels(dge$samples$group)))
    group.col <- colors[dge$samples$group] 

    out.png <- paste(folder, '/edgeR/img/edgeR_log2_cpm.png', sep='')
    as.png( {
            par(mfrow=c(1,1))
            boxplot(logcpm, xlab='', ylab=' Log2 counts per million', 
	        col= group.col, las=2, outline=FALSE)
            abline(h=median(logcpm), col='red')
            graphics::title('Boxplots of logCPMs unnormalised')
        }, out.png )

    # Now we produce an MDS plot to see any significant difference between the
    # groups
    out.png <- paste(folder, '/edgeR/img/edgeR_mds_plot.png', sep='')
    as.png( {
            par(mfrow= c(1,1))
            plotMDS(dge, col = group.col)
        }, out.png)


    # Next, we will identify the to differentially expressed genes (higher variance).
    # We apply a funcion that calculates the variance by rows (genes) and then
    # retrieve the n.genes most DE genes

    cat(paste("\n\tCalculating top", n.genes, "genes with the highest variance\n\n"))   
	var_genes <- apply(logcpm, by.rows, var)
    select_var <- names(sort(var_genes, decreasing=TRUE))[1:n.genes]
    highly_var <- logcpm[select_var,]
    #dim(highly_var)

    # we plot the heatmap.2 (gplots) without a line (trace), scale by row
    # (difference in color) margins (something about the labels used), also we
    # reverse the colors because by default the red is associated with low
    # expression and we are more used to it meaning "hot"
    out.png <- paste(folder, '/edgeR/img/edgeR_heatmap.png', sep='')
    as.png( {
            #margins <- par("mar")
            #par(mar=c(25, 5, 5, 10))
            heatmap.2(highly_var, 
                      col= rev(morecolors(50)), 
                      trace='none', 
                      ColSideColors = group.col,
                      scale='row', 
                      margins= c(15,5))
            #par(mar=margins)
        }, out.png)

    # This plot is of limited use. We'd better have other names for rows and 
    # columns and plot the top 50 (most variable) to see them well
	n <- 50
	high_var <- highly_var
	colnames(high_var) <- gsub("_R1.[Bb][Aa][Mm]", "", colnames(highly_var))
    
	## NOTE: we should add the annotation here, before doing the next plot
    if (all(str_sub(rownames(high_var), 1, 3) == "ENS")){
        by='GENEID'
    }else{
    	by='GENENAME'}
	genenames <- ensembldb::select(ens.db, keys = rownames(high_var), 
                	    keytype = by, 
                	    columns=c('GENENAME', 'GENEID'))
	out.png <- paste(folder, '/edgeR/img/edgeR_heatmap.', n,'.png', sep='')
	as.png( {
        	#margins <- par("mar")
        	#par(mar=c(25, 5, 5, 10))
        	heatmap.2(high_var[1:n,], 
                    	  col=rev(morecolors(n)), 
                    	  trace='none', 
                    	  ColSideColors = group.col,
                    	  scale='row', 
                    	  margins= c(15,5),
                    	  labRow=genenames[1:n, 1])
        	#par(mar=margins)
			}, out.png )

    # We can automate the estimation of the dispersion and add it to the dge object
    cat("\n\tEstimating Dispersion\n")
	cat("\n\tGenerating DGEList object\n\n")

	dge <- estimateCommonDisp(dge)

    # And now we can estimate gene-wise dispersion estimates allowing for a
    # possible trend with average count size. These will allow us to use a GLM
    # instead of a plain LM.
    # This gives us the BCV (Biological Coefficient of Variation) between samples
    dge <- estimateGLMTrendedDisp(dge)
    dge <- estimateTagwiseDisp(dge)

    out.png <- paste(folder, '/edgeR/img/edgeR_BCV_dispersions.png', sep='')
    as.png(plotBCV(dge), out.png)

    return(dge)
}

eR.voom.variation.analysis <- function(dge, folder)
{
     ############## JR #########################
    # 
    # edgeR user's guide
    # Chapter 3 
    # 
    # Specific experimental designs 3.1 Introduction In this chapter, we outline
    # the principles for setting up the design matrix and forming contrasts for
    # some typical experimental designs.
    # 
    # Throughout this chapter we will assume that the read alignment, normalization
    # and dispersion estimation steps described in the previous chapter have
    # already been completed. We will assume that a DGEList object y has been
    # created containing the read counts, library sizes, normalization factors and
    # dispersion estimates.
    # 
    # 3.2 Two or more groups
    # 
    # 3.2.1 Introduction
    # 
    # The simplest and most common type of experimental design is that in which a
    # number of experimental conditions are compared on the basis of independent
    # biological replicates of each condition. Suppose that there are three
    # experimental conditions to be compared, treatments A, B and C, say. The
    # samples component of the DGEList data object might look like:
    # 
    # > y$samples
    # group lib.size norm.factors
    # Sample1 A 100001 1
    # Sample2 A 100002 1
    # Sample3 B 100003 1
    # Sample4 B 100004 1
    # Sample5 C 100005 1
    # 
    # Note that it is not necessary to have multiple replicates for all the
    # conditions, although it is usually desirable to do so. By default, the
    # conditions will be listed in alphabetical order, regardless of the order that
    # the data were read:
    # 
    # > levels(y$samples$group)
    # [1] "A" "B" "C"
    # 29
    # 
    # 3.2.2 Classic approach
    # 
    # The classic edgeR approach is to make pairwise comparisons between the
    # groups. For example,
    # 
    # > et <- exactTest(y, pair=c("A","B"))
    # > topTags(et)
    # 
    # will find genes differentially expressed (DE) in B vs A. Similarly
    # 
    # > et <- exactTest(y, pair=c("A","C"))
    # 
    # for C vs A, or
    # 
    # > et <- exactTest(y, pair=c("C","B"))
    # 
    # for B vs C.
    # 
    # Alternatively, the conditions to be compared can be specified by number, so
    # that
    # 
    # > et <- exactTest(y, pair=c(3,2))
    # 
    # is equivalent to pair=c("C","B"), given that the second and third levels of
    # group are B and C respectively.
    # 
    # Note that the levels of group are in alphabetical order by default, but can
    # be easily changed.
    # 
    # Suppose for example that C is a control or reference level to which
    # conditions A and B are to be compared. Then one might redefine the group
    # levels, in a new data object, so that C is the first level:
    # 
    # > y2 <- y
    # > y2$samples$group <- relevel(y2$samples$group, ref="C")
    # > levels(y2$samples$group)
    # [1] "C" "A" "B"
    # 
    # Now
    # 
    # > et <- exactTest(y2, pair=c("A","B"))
    # 
    # would still compare B to A, but
    # 
    # > et <- exactTest(y2, pair=c(1,2))
    # 
    # would now compare A to C.
    # 
    # When pair is not specified, the default is to compare the first two group
    # levels, so
    # 
    # > et <- exactTest(y)
    # 
    # compares B to A, whereas
    # 
    # > et <- exactTest(y2)
    # 
    # compares A to C.
    # 
    # 
    # 
    # 3.2.3 GLM approach
    # 
    # The glm approach to multiple groups is similar to the classic approach, but
    # permits more general comparisons to be made. The glm approach requires a
    # design matrix to describe the treatment conditions. We will usually use the
    # model.matrix function to construct the design matrix, although it could be
    # constructed manually. There are always many equivalent ways to define this
    # matrix. Perhaps the simplest way is to define a coefficient for the
    # expression level of each group:
    # 
    # > design <- model.matrix(~0+group, data=y$samples)
    # > colnames(design) <- levels(y$samples$group)
    # > design
    # A B C
    # Sample1 1 0 0
    # Sample2 1 0 0
    # Sample3 0 1 0
    # Sample4 0 1 0
    # Sample5 0 0 1
    # attr(,"assign")
    # [1] 1 1 1
    # attr(,"contrasts")
    # attr(,"contrasts")$group
    # [1] "contr.treatment"
    # 
    # Here, the 0+ in the model formula is an instruction not to include an
    # intercept column and instead to include a column for each group.
    # 
    # One can compare any of the treatment groups using the contrast argument of
    # the glmQLFTest or glmLRT function. For example,
    # 
    # > fit <- glmQLFit(y, design)
    # > qlf <- glmQLFTest(fit, contrast=c(-1,1,0))
    # > topTags(qlf)
    # 
    # will compare B to A. The meaning of the contrast is to make the comparison
    # -1*A + 1*B + 0*C, which is of course is simply B-A.
    # 
    # The contrast vector can be constructed using makeContrasts if that is
    # convenient. The above comparison could have been made by
    # 
    # > BvsA <- makeContrasts(B-A, levels=design)
    # > qlf <- glmQLFTest(fit, contrast=BvsA)
    # 
    # One could make three pairwise comparisons between the groups by
    # 
    # > my.contrasts <- makeContrasts(BvsA=B-A, CvsB=C-B, CvsA=C-A, levels=design)
    # > qlf.BvsA <- glmQLFTest(fit, contrast=my.contrasts[,"BvsA"])
    # > topTags(qlf.BvsA)
    # > qlf.CvsB <- glmQLFTest(fit, contrast=my.contrasts[,"CvsB"])
    # > topTags(qlf.CvsB)
    # > qlf.CvsA <- glmQLFTest(fit, contrast=my.contrasts[,"CvsA"])
    # > topTags(qlf.CvsA)
    # 
    # which would compare B to A, C to B and C to A respectively.
    # 
    # 
    # Any comparison can be made. For example,
    # 
    # > qlf <- glmQLFTest(fit, contrast=c(-0.5,-0.5,1))
    # 
    # would compare C to the average of A and B. Alternatively, this same contrast
    # could have been specified by
    # 
    # > my.contrast <- makeContrasts(C-(A+B)/2, levels=design)
    # > qlf <- glmQLFTest(fit, contrast=my.contrast)
    # 
    # with the same results.
    # 
    ############## NOTE END #########################

    # we know that the variability we see in the expression depends on infection so
    # we have to take this into account and create a model
    # 0+ forces the design to include all groups and not have
    # an intercept (reference) column

	#target <- as.factor(metadata[, design.column])
    #design <- model.matrix(~  target)
    #colnames(design) <- levels(target)
    cat("\n\tGenerating Design Matrix\n\n")

    design <- model.matrix(~ dge$samples$group)
    colnames(design) <- c("Intercept", levels(dge$samples$group)[-1])
	rownames(design) <- rownames(dge$samples)
    if (VERBOSE) print(design)

    # IF WE DO NOT INCLUDE THE "~ 0 + " IN THE FORMULA, MODEL.MATRIX()
    # WILL USE AS REFERENCE THE FIRST ALPHABETICAL ORDER LEVEL !!!
    # 
    # Another problem is that we are limited to the comparisons defined by
    # coefficents, if we want more control we need to use makeContrasts
    #
	
	### Let us test for differential expression with the DGE data
	# 
	# 1. EDGER Classical Approach: Negative Bionomial Linear Models
    # ==============================================================
	# First, fit genewise NB GLMs and compute the coefficients for
	# the different levels of our design column.
    cat("\n\tComputing Model Coefficients\n\n")
	gfit <- glmFit(dge, design)
    if (VERBOSE) {
        #print(names(gfit))
        print(head(coef(gfit)))
    }
    
    # We now conduct Likelihood Ratio tests and show the top genes
    # for the selected comparison coefficients (reference vs. coeff)
    cat("\n\tPerforming Likelihood Ratio Tests\n\n")
    for (i in 1:ncol(gfit)){
		lrt <- glmLRT(gfit, coef = i)	# coef = 1... length(gfit$coefficients)
    	print(topTags(lrt))
		cat("\n")
    }
	# the problem here us that we are limited to comparisons defined by
    # the fitting coefficients (ref vs. variable-in-coeff).
    # If we want more control we need to use makeContrasts()
    # This is due to the formula used:  ~ table[ , design.column ]
    # If we used ~ 0 + table[ , design.column] we would have diffs for
    # all levels, but no reference.
    
	# We can perform a similar analysis to glmFit, but also estimating
	# genewise QuasiLikelihood (QL) dispersion values. Likelihood Ratio
	# Test (LRT) are replaced by Bayes Quasilikelihood F-Tests.
	cat("\n\tPerforming Quasi-likelihood F-Tests\n\n")
	qfit <- glmQLFit(dge, design)
    if (VERBOSE) {
        #print(names(gfit))
        print(head(coef(qfit)))
    }
    
    # We now conduct Likelihood Ratio tests and show the top genes
    # for the selected comparison coefficients (reference vs. coeff)
    for (i in 1:ncol(qfit)){
		qlt <- glmQLFTest(qfit, coef = i)	# coef = 1... length(gfit$coefficients)
    	print(topTags(qlt))
		cat("\n")
    }

	
	
	# 2. Limma - Voom Approach: Negative Bionomial Linear Models
    # ==============================================================
	# First, estimate Mean - Variance relationship if the count data. Use trend
	# to assign wheights to each observation, adjusting for heteroscedasticity.
    #
    cat("\n\tComputing Mean - Variance Trend\n\n")
    v <- voom(dge, design, plot = FALSE)
    out.png <- paste(folder, '/edgeR/img/edgeR_voom.png', sep='')
    as.png( {
            par(mfrow= c(1,1))
            voom(dge, design, plot = TRUE)
        }, out.png )

    # Carry on a differential expression analysis using the VOOM transformed data
    #	
	# With Limma, we use the weighted data and the design matrix to fit
	# gewewise linear models and compute the coefficients.
    cat("\n\tFitting Weighted Linear Models\n\n")
    fit <- lmFit(v, design)
	if (VERBOSE) cat('Coefficients:\n') ; print(head(fit$coefficients))

    # EBayes uses empirical Bayes to squeze genewise dispersion towards
	# the global trend. Then, computes moderated t-statistics and
	# moderate F-statistics to test differential expression
	cat("\n\tComputing Empirical Bayes Statistics\n\n")
	fit <- eBayes(fit)

	# Next, the results of the different testing strategies are contrasted
	# to determina differentially expressed genes.
	results <- decideTests(fit)
    if (VERBOSE) {
		print(summary(results))
    	for (i in 1:ncol(fit)){
			cat("\n\n", colnames(fit)[i], ":\n")
			print(topTable(fit, coef = i, sort.by='p'))
		}
    }

    eR.save.top(fit, folder, n.genes, sort.by = c('p', 'B', 'logFC', 'AveExpr'))

    # the problem here is that we are limited to the comparisons defined by
    # coefficents, if we want more control we need to use makeContrasts
    #
    # but for now we will leave it here, as we will do it in more detail
    # with DESeq2

    return(fit)

}

eR.fit.annotate.save <- function(fit,
                            ens.db,
                            biomart.ann,
                            folder,
                            n.genes = 1000	# top N genes to save in tables
                            )
{
    #---------------------------------------------------------------------------
    #Now is time to connect all the results we have with the existing
    #information  we know from the literature, so we will retrieve infos from
    #ENSEMBL and connect with the genes  we have kept. We are interested only in
    #genes and transcriptomes (non-characterized genes) from the organism
    #of interest

    # First we need to identify wether the gene identifiers used for
	# differential analysis correspond to GENEID or GENENAME fields in
	# Ensembl notation. 
	
	if (all(str_sub(rownames(fit), 1, 3) == "ENS")){
        by='GENEID'
    }else{
    	by='GENENAME'}
		   
    # Once we have prepared the database and identified the Key type
	# we extract the annotation for the genes in our Differentially
	# Expressed gene list ('fit').
    ens.ann <- ensembldb::select(ens.db, 
                      column = by, keytype = by, keys = rownames(fit), 
                      columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION', 
                                 'GENENAME', 'GENEID', 'ENTREZID',
                                 'TXID', 'TXBIOTYPE',
                                 'PROTEINID', 'UNIPROTID'
                                ))

    #ann <- ensembldb::select(ens.db, 
    #              keytype= 'GENEID', keys=rownames(fit), 
    #              columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION',
    #                         'GENENAME', 'GENEID', 'ENTREZID', 
    #                          'TXNAME', 'TXBIOTYPE', 
    #                          'PROTEINID', 'UNIPROTID'))



    # SAVE ANNOTATION
    # ---------------
    # this is all the annotation for all the genes in 'fit'
    write.table(ens.ann, file=paste(folder, '/annotation/ensembl.annotation.txt', sep=''), 
	    sep='\t', row.names=T, col.names=T)

    # check if the amount of genes we have is the same as the number of 
    # the annotations that we have extracted
    if ( ! length(ens.ann[, by])==length(rownames(fit)) ) {
        cat("\nNumber of annotations does not match gene number")
        cat("\nUsing only one entry per gene\n")
        ens.ann.1 <- ens.ann[ ! duplicated(ens.ann$GENEID), ]
    } else {
        ens.ann.1 <- ens.ann
    }
    write.table(ens.ann.1, file=paste(folder, '/annotation/ensembl.annotation.1st.txt', sep=''), 
	    sep='\t', row.names=T, col.names=T)


    # As long as it works perfectly we continue by connecting the 
    # annotation information to the fit using a new bucket named "genes"
    # (fit is an object of class "MArrayLM" (package "limma") which is a
    # list. I.e. we assign the annotation to a new element named 'genes'  
    # of this list.
    #
    #	This checks that genes and annotation go in the same order)
    fit$genes <- ens.ann.1[ match(rownames(fit), ens.ann.1[, by], nomatch = NA), ]

    
    # and now we will also add to 'fit$genes' the biomaRt annotation
    fit$genes <- merge(fit$genes, bm.annot.1, by.x = "GENEID", by.y = "ensembl_gene_id")
    if (VERBOSE) print(head(fit$genes))

    # SAVE THE FULL ANNOTATED FIT
    # ---------------------------
    # save all the contents of 'fit' in an RDS file
    saveRDS(fit, file = paste(folder, '/edgeR/annotatedVOOMfit.rds', sep=''))
    #	'fit' can later be recovered with: 
    #		fit <- readRDS(file=paste(folder, '/annotatedVOOMfit.rds', sep=''))
    # and save as well as Rdata file
    save(fit, file = paste(folder, '/edgeR/annotatedVOOMfit.RData', sep=''))
    #	'fit' can later be recovered with: 
    #		fc <- load(file=paste(folder, '/annotatedVOOMfit.RData', sep=''))


    # save top annotated genes
    eR.save.top.annotated(fit, folder, n.genes, 
                          sort.by = c('p', 'logFC', 'AveExpr'))
    

    return(fit)

}


eR.fit.treat <- function(fit,
                         threshold = 1,
                         folder,
                         verbose = T )
{
    # Similar to Empitical Bayes but Testing relative to a LFC threshold.
	# LFC = 1 means a 2x fold change. This is more strict than eBayes.
    fit.thres <- treat(fit, lfc = threshold)
    if (VERBOSE) {
        res.thres <- decideTests(fit.thres)
        print(summary(res.thres))
    }
    if (VERBOSE)
        print(topTreat(fit.thres, coef = 1, sort.by='p'))
    
    # save topTreat data
    eR.save.top.treat.ann(fit.thres, folder, n.genes, 
                          sort.by=c('p', 'logFC', 'AveExpr'))

    # also save fitted to a threshold as tables
    eR.save.top.ann.thresh(fit.thres, folder, n.genes, 
                           sort.by=c('p', 'logFC', 'AveExpr'))
    
    return(fit.thres)
}

eR.dge.all.comparisons <- function(dge, 
                               design.column,
                               ens.db,
                               biomart.db,
                               folder)
{
	cat("\n\tPerforming all pairwise comparisons by:", toupper(design.column), "\n")
    # --------------------------------------------------------------------------------
    # Do all possible pairwise comparisons.
	# We need to define a new design without intercept, since we do not want to
	# fix a reference level.

    design <- model.matrix(~0 + as.factor(dge$samples$group))
    grps <- levels(as.factor(dge$samples$group))
	colnames(design) <- grps
    n.grps <- length(grps)
	if (VERBOSE) cat("\nGroups:", grps, "\n")
	
	## Fit new genewise multiple linear models with quasilikelihood varianc
	## estimation.
	cat("\n\tComputing Quasi-Likelihood Dispersion\n\n")
    qlfit <- glmQLFit(dge, design = design, robust = TRUE, abundance.trend = TRUE)
    out.png <- paste(folder, '/edgeR/img/edgeR_QLFit_', design.column, '.png', sep='')
    as.png(plotQLDisp(qlfit), out.png)

    # if we wanted to apply a log-fold-change theshold, we could do
    # the calculations using qlfit here before testing for contrasts 
    # using e.g.
    # qlf.cmp <- glmTreat(glfit, cef=1..ncol(glfit$design), lfc=threshold)
    # or
    # qlf.cmp <- glmTreat(glfit, contrast=cmp, lfc=threshold)
    #

    eR.data <- list()
    for (i in grps) {
        for (j in grps) {
            if (i == j) next	# ignore self-comparisons
            formula <- paste(i, '-', j)
            cat("\n\n\tComputing DGE:", i, '-', j, '\n\n')
            cmp <- makeContrasts(formula, levels = design)
            # glmQLFTest is similar to glmLRT except it uses Bayes quasi-likelihood
            # the P-values are always >= those produced by glmLRT
            qlf.cmp <- glmQLFTest(qlfit, contrast = cmp)
            # or if a threshold has been defined
            #qlf.cmp <- glmTreat(glfit, contrast=cmp, lfc=threshold)

            # get summary of up/down regulated genes
            if (VERBOSE) {
                print(summary(decideTests(qlf.cmp)))
            }
            png.file <- paste(folder, "/edgeR/img/edgeR_QLF_MD_", 
                              design.column, '_', i, '-', j, '.png', sep='')
            as.png( { 
                    plotMD(qlf.cmp)
                    abline(h = c(-1, 1), col="darkgreen")
                    }, png.file)
                    
            if (INTERACTIVE) {
                plotMD(qlf.cmp)
                abline(h=c(-1, 1), col="darkgreen")
                continue.on.enter("Press [RETURN] to continue: ")
            }
            
            # ANNOTATE the results without saving them
            ### NOTE consider using eR.fit.annotate.ensembl.biomart
            qlf.cmp <- eR.fit.annotate(qlf.cmp, ens.db = ens.db, biomart.db = biomart.db)

            name <- paste(folder, '/edgeR/fit_', design.column, '_', i, '_-_', j, '_annot', sep='')
            # qlf.cmp is a list of tables, if we want to save it,
            # we'll need to save the whole object
            # will add .rds and .RData to the files created
            #
            eR.save.fit(qlf.cmp, name)
            # defaults: n.genes=500, sort.by='PValue', p.valu=0.01
            # we'll save all (<=100.000) significant genes

            name <- paste(folder, '/edgeR/comp_', design.column, '_', i, '_-_', j, '_annot', sep='')
	    	# qlf.cmp$table is a table with logFC, logCPM, F and PValue
            # that is what weill be saved when using topTags and write.table
            # will add "_top_" n "_by_" sort.by
            # defaults: n.genes=500, sort.by='PValue', p.value=0.01
            eR.save.top.fit(qlf.cmp, file=name, n.genes=100000)

            # we can use limma to test for over-representation of gene
            # ontology (GO) terms or KEGG pathways with goana() or kegga()
            # using the entrez.gene.ids of DE genes, to an FDR of 0.05 (default)
	    # we can use species.KEGG="gga" or "cjo"
	    # ( see https://www.kegg.jp/kegg/catalog/org_list.html )
	    #
            # eR.go <- goana(qlf.cmp, species="Cj")
	    	# eR.go <- goana(qlf.cmp, species="Gg")
	    	# eR.kegg <- kegga(qlf.cmp, species.KEGG="cjo")
	    	# eR.kegg <- kegga(qlf.cmp, species.KEGG="gga")
	    	# topGO(go, sort="up", number=n.genes)
	    	# topKEGG(keg, sort="up", number=n.genes)

            qlf.result <- list(
	                      eR.cmp = qlf.cmp
                          #, eR.go=eR.go
			      			#, eR.kegg=eR.kegg
			      			)
            eR.data[[formula]] <- qlf.result 
            #print(topTags(qlf.cmp))
        }
    }
    
    return(eR.data)
}



##### DESEQ2



ds2.annotate.results <- function(dds, ens.db, org.db)
{
    if (all(str_sub(rownames(dds), 1, 3) == "ENS")){
        by='GENEID'
    }else{
    	by='GENENAME'}

    # annotate with ensembl ens.db

    ens.ann <- ensembldb::select(ens.db, 
                  keytype= by, keys=rownames(dds), 
                  columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION',
                             'GENENAME', 'GENEID', 'ENTREZID', 
                              'TXNAME', 'TXID', 'TXBIOTYPE', 
                              'PROTEINID', 'UNIPROTID'))

    ens.ann.1 <- ens.ann[ ! duplicated(ens.ann$GENEID), ]

    # Add annotation to dds to keep everything in one place
    #dds$ens.annot.1 <- ens.ann.1
    #dds$bm.annot.1 <- b,.annot.1

    #
    # Extract annotation using Org object if avaiable
    #
    if ( ! is.null(org.db) ) {

        # and now we should be able to use the Org package if we successfully built it
        # at the beginning.
        geneSymbols <- ens.ann$ENTREZID

        # retrieve go ids
        go.ann <- AnnotationDbi::select(org.db, 
	        	keys = as.character(ens.ann$ENTREZID),
                columns=c("ENTREZID", "GO", "GOALL", "ONTOLOGY","ONTOLOGYALL"), 
                keytype="ENTREZID", 
                multiVals="CharacterList")

        # retrieve corresponding descriptions
        GOdescription <- AnnotationDbi::select(GO.db, keys = go.ann$GO, 
                         columns= c("GOID", "TERM", "DEFINITION", "ONTOLOGY"), 
                         keytype= "GOID")

    }
    
    return(list(ens.ann, ens.ann.1, geneSymbols, go.ann, GOdescription))

}



ds2.dds.compare.annotate.save <- function(dds, 
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
                      ensembl.ann[ match(cmp.df.a$ensembl.gene.id, ensembl.ann$GENEID, nomatch = NA), ])
    # bm.annot fails to annotate some entries (why?)
    cmp.df.a <- cbind(cmp.df.a, 
                      biomart.ann[ match(cmp.df.a$ensembl.gene.id, biomart.ann$ensembl_gene_id, nomatch = NA), ])

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




ds2.dds.plot.and.save <- function(dds, 
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
                ensembl.ann[ match(raw.df$ensembl.gene.id, ensembl.ann$GENEID, nomatch = NA), ]
                )
    # bm.annot fails to annotate some entries (why?)
    raw.df <- cbind(
                raw.df, 
                biomart.ann[ 
                  match(raw.df$ensembl.gene.id, biomart.ann$ensembl_gene_id, nomatch = NA),
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
    out.png <- paste(folder, '/DESeq2/img/DESeq2_raw', out.base, '_hist.png', sep='')
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
                ensembl.ann[ match(srt.df$ensembl.gene.id, ensembl.ann$GENEID, nomatch = NA), ])
    srt.df <- cbind(srt.df, 
                biomart.ann[ match(srt.df$ensembl.gene.id, biomart.ann$ensembl_gene_id, nomatch = NA), ])
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


ds2.interactively.print.n.most.significant.de.genes <- function(dds.results, n=10) 
{
    continue.on.enter("You may want to maximize your terminal before continuing ")

    options(width = 200)
    n <- 10
    for (i in names(dds.results)) {
        cat("\n\n\tMost significant", n, "genes for", i, "\n\n")
        print(head(dds.results[[i]]$signif[ , c("log2FoldChange", "entrezgene_description")] ), n)
        continue.on.enter("Press [ENTER] to continue ")
    }
    options(width = 80)
    continue.on.enter("Done, you can restore your terminal now ")
}


ds2.interactively.plot.top.up.down.regulated.gene <- function(dds,
                                                              dds.results, 
                                                              design.column)
{
    # Plot counts of the gene with maximal l2FC (up or down)
    threshold <- 0
    for (n in names(dds.results)) {
        res <- dds.results[[n]]$signif
        up <- res[ res$log2FoldChange > threshold, ]	# up regulated
        down <- res[ res$log2FoldChange < -threshold, ]	# down regulated

        most.up <- rownames(up)[which.max(up$log2FoldChange)]	# max up
        cat("\ncounts of most overexpressed gene in", n, ":\n",
            most.up, 
            up$log2FoldChange[which.max(up$log2FoldChange)], 
            up$entrezgene_description[which.max(up$log2FoldChange)], 
            "\n")
        #plotCounts(dds, gene=rownames(res)[which.max(res$log2FoldChange)], intgroup="PFU")
        print(plotCounts(dds, gene=most.up, intgroup=design.column))
        k <- continue.on.key()
        if (k == "q") break

        most.down <- rownames(down)[which.min(down$log2FoldChange)]	# min down
        cat("\ncounts of most underexpressed gene in", n, ":\n",
            most.down, 
            down$log2FoldChange[which.min(down$log2FoldChange)], 
            down$entrezgene_description[which.min(down$log2FoldChange)],
            "\n")
        #plotCounts(dds, gene=rownames(res)[which.min(res$log2FoldChange)], intgroup="PFU")
        print(plotCounts(dds, gene=most.down, intgroup=design.column))
        continue.on.key()
        if (k == "q") break
    }
    cat("\n")

}



# 1) ANNOTATE RESULTS

annotateDESeqResults <- function(results,
				ensembl.db,
				biomart.db,
				org.db = NULL,
				annotation.dir,		# For saving annotation dataset
				save = TRUE,
				out.file.base,
				out.dir ){
	
	cat("\n	Annotating genes ...\n\n")
	
	## Get gene annotation from ENSEMBL databsase and BioMart database

	if ( ! file.exists(paste(annotation.dir, 'ensembl.annotation.1st.txt', sep = '/'))) {
		ensembl.ann <- annotateGenesFromENSEMBL(ensembl.db = ensembl.db,
							gene.ids = rownames(results),
							save = save,
							out.dir = annotation.dir)
	} else {
		ensembl.ann <- read.table(paste(annotation.dir, 'ensembl.annotation.1st.txt', sep = '/'))
	}
	
	if ( ! is.null(org.db) ){
		## Get gene annotation from ENSEMBL databsase and BioMart database
		if ( ! file.exists(paste(annotation.dir, 'orgDb.GO.annotation.1st.txt', sep = '/'))) {
			org.ann <- annotateGenesFromORG(org.db = org.db,
								gene.ids = rownames(results),
								save = save,
								out.dir = annotation.dir)
		} else {
			org.ann <- read.table(paste(annotation.dir, 'orgDb.GO.annotation.1st.txt', sep = '/'))
		}
	}
	
	# Remove duplicated records preserving only first entry

	ensembl.ann <- ensembl.ann[ ! duplicated(ensembl.ann$GENEID), ]
	biomart.ann <- biomart.db[ ! duplicated(biomart.db$ensembl_gene_id), ]
	if ( ! is.null(org.db) ) org.ann <- org.ann[ ! duplicated(org.ann$ENTREZID), ]


	## Add annotation to DESeq results data.frame to keep everything in one place

	# Add row names (ensembl.gene.id) as an additional column
	res.ann <- results
	res.ann$gene.id <- rownames(res.ann)
	
	### Identify keys for ENSEMBL ###
	if (all(str_sub(res.ann$gene.id, 1, 3) == "ENS")){
        by='GENEID'
    }else{
    	by='GENENAME'}

	# ENSEMBL GENENAMES match the GENE.IDs from our results dataframe
	#res.ann <- cbind(res.ann, ensembl.ann[ match(res.ann$ensembl.gene.id, ensembl.ann$GENENAME), ])
	for (i in colnames(ensembl.ann)){
		res.ann[i] <- ensembl.ann[match(res.ann$gene.id, ensembl.ann[,by], nomatch = NA), i]
	}	
	
	# BIOMART ENTREZ ACCESSION match the GENE.IDs from our results dataframe
	#res.ann <- cbind(res.ann, biomart.ann[ match(res.ann$ensembl.gene.id, biomart.ann$entrezgene_accession), ])
	for (i in colnames(biomart.ann)){
		res.ann[i] <- biomart.ann[match(res.ann$GENEID, biomart.ann$ensembl_gene_id, nomatch = NA), i]
	}

	if ( ! is.null(org.db) ){

		### Identify keys for NCBI (OrgDb Package) ###
		if (all(str_sub(res.ann$gene.id, 1, 3) == "ENS")){
        	by = 'ENSEMBL'
    	} else if (! any(is.na(as.numeric(res.ann$gene.id)))){
    		by = 'ENTREZID'
		} else {
			by = 'SYMBOL'}

		# NCBI SYMBOL match the GENE.IDs from our results dataframe
		for (i in colnames(org.ann)){
			res.ann[i] <- org.ann[match(res.ann$gene.id, org.ann[,by], nomatch = NA), i]
		}
	}

	## Save annotated results as data frame (table)
	if (save){
		write.table(data.frame(res.ann), sep = "\t",
					file = paste(out.dir, '/DESeq2/', out.file.base, '_ann_results.tab', sep = ''),
					row.names = TRUE, col.names = TRUE)

		saveRDS(res.ann, paste(out.dir, "/DESeq2/", out.file.base, "_ann_results.rds", sep=""))
	}
	return(res.ann)
}

# 2) COMPARE AND ANNOTATE RESULTS FROM DESEQ2 OBJECT
ddsCompareAnnotate <- function(	dds, 
								contrast,			# c(factor, numerator, denominator) 
								filterFun = ihw,	# IHW increases statistical power
								alpha = 0.01,
								annotate = TRUE, 
								ensembl.db,
								biomart.db,
								org.db = NULL,
								annotation.dir,		# For saving annotation tables
								save = TRUE,
								out.file.base,
								out.dir ) {

	results <- results(dds,
					contrast = contrast, 
					pAdjustMethod="BH",
					filterFun = filterFun,
					alpha = alpha)
	
	print(head(results))
	print(mcols(results)$description)	# Description of the Results fields

    ## Provide default output base name
    if (missing(out.file.base)) {
		out.file.base <- paste(contrast[1], contrast[2], 'vs', contrast[3], sep = '_')
	}

	if (save){
	
		## Save comparison results
		write.table(data.frame(results),
				file = paste(out.dir, '/DESeq2/raw/', out.file.base, '_raw_results.tab', sep = ''),
				row.names = TRUE, col.names = TRUE)

		saveRDS(results, paste(out.dir, "/DESeq2/raw/", out.file.base, "_raw_results.rds", sep=""))
	
		## Save summary
		sink(paste(out.dir, '/DESeq2/', out.file.base, '_summary.txt', sep=''), split=T)
		summary(results)
		sink()
	}
	
	if (annotate){
		results <- annotateDESeqResults(results = results,
										ensembl.db = ensembl.db,
										biomart.db = biomart.db,
										org.db = org.db,
										annotation.dir = annotation.dir,
										save = save,
										out.file.base = out.file.base,
										out.dir = out.dir )
	}
	return(results)
}



# 3) IDENTIFY TOP N DIFFERENTIALLY EXPRESSED GENES
top_n <- function(ds, n, col = "padj", decreasing = F) {
	
	df <- as.data.frame(ds)
	ordered.idx <- order(abs(df[,col]), decreasing = T)
	ordered.df <- df[ordered.idx, ]
	
	if (dim(ordered.df)[1] < n) n <- dim(ordered.df)[1]
	top.n <- ordered.df[1:n, ]
	
	return(top.n)
}


# 4) PERFORM DESEQ2 COMPARISON, ANNOTATE AND PLOT THE RESULTS

analysePlotDESeq <- function(	#results,
								dds, 
								contrast,			# c(factor, numerator, denominator) 
								filterFun = ihw,	# IHW increases statistical power
								alpha = 0.01,
								shrnk.type,
								annotate = TRUE,
								ensembl.db,
								biomart.db,
								org.db = NULL,
								annotation.dir,		# For saving annotation data
								save = TRUE,
								overwrite = FALSE,
								out.file.base,
								out.dir ) {

    ## Provide default output base name
    if (missing(out.file.base)) {
		out.file.base <- paste(contrast[1], contrast[2], 'vs', contrast[3], sep = '_')
	}

    ## Check if final comparison results already exist. If so 
	if ( ! file.exists(paste(out.dir, "/DESeq2/", out.file.base, "_final_result_list.rds", sep = "")) || (overwrite == TRUE)){

		## Check if DESeqResults object is provided
		
		results <- ddsCompareAnnotate(	
							dds = dds, 
							contrast = contrast,
							filterFun = ihw,
							alpha = alpha,
							annotate = annotate, 
							ensembl.db = ensembl.db,
							biomart.db = biomart.db,
							org.db = org.db,
							annotation.dir = annotation.dir,	
							save = save,
							out.file.base = out.file.base,
							out.dir = out.dir)

		## Histogram plot p-value distribution
		if (save) {
			out.png <- paste(out.dir, '/DESeq2/img/DESeq2_', out.file.base, '_padj_hist.png', sep='')
		} else out.png <- NULL

		as.png( {
				margins <- par("mar")
				par(mar = c(5, 5, 5, 5))
				hist(results$padj, main=paste("Adj. p-value distribution", "\n", out.file.base, sep = ''),
				breaks = 1/alpha, xlab = "Adj. p-value")
				par(mar = margins)
				},
		out.png, overwrite = TRUE)
		
		
		## Shrinkage of data to improve visualization and ranking
		cat("\n\tShrinking Log 2 Fold Change ...\n\n")    
		coef <- paste(contrast[1], contrast[2], 'vs', contrast[3], sep = '_')
		shrunk.lfc <- lfcShrink(dds = dds, coef = coef, type = shrnk.type)
		
		if (save) {
		    out.png <- paste(out.dir, "/DESeq2/img/DESeq2_", out.file.base, "_shrunk_", shrnk.type, "_MA.png", sep='') 
		} else out.png <- NULL
		
		as.png( {

			dim <- par("mfrow")
			par(mfrow = c(2,1))
		   
			# plotMA shows the log2 fold changes attributable to a given variable
			# over the mean of normalized counts for all the samples in the dataset
			# Points above alpha are colored, outliers are shown as directed triangles    
			DESeq2::plotMA(results, alpha = alpha, main = "Log2 Fold Change")
			
			# it is useful to look at the shrunk l2fc values, which removes the noise
			# associated with l2fc changes from low-count genes
			DESeq2::plotMA(shrunk.lfc, alpha = alpha, main = paste("Shrunken LFC (", shrnk.type, ")", sep = ""))
			
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
			par(mfrow = dim)
			},
			
		out.png, overwrite = TRUE)
		
		## Save shrunken results
		if (save) {
		
		    out.file <- paste(out.dir, "/DESeq2/shrunk/", out.file.base, "_shrunk_", shrnk.type, "_raw_results.rds", sep='')
		    cat("\n\tSaving", out.file, '\n\n')
		    saveRDS(shrunk.lfc, file=out.file)
		
			out.file <- paste(out.dir, "/DESeq2/shrunk/", out.file.base, "_shrunk_", shrnk.type, "_raw_results.tab", sep='')
			cat("\tSaving", out.file, '\n')
			write.table(shrunk.lfc, 
				        file=out.file,
				        row.names=T, col.names=T, sep='\t')
		}
		
		## Annotate shrunken data for visualization in MA plot	
		if (annotate) {
			
			ann.shrunk <- annotateDESeqResults(results = shrunk.lfc,
									ensembl.db = ensembl.db,
									biomart.db = biomart.db,
									org.db = org.db,
									annotation.dir = annotation.dir,
									save = FALSE)	
			if (save) {

			## Save shrunken annotated results		
				out.file <- paste(out.dir, "/DESeq2/shrunk/", out.file.base, "_shrunk_", shrnk.type, "_ann_results.rds", sep='')
				cat("\n\tSaving", out.file, '\n\n')
				saveRDS(ann.shrunk, file=out.file)
			
				out.file <- paste(out.dir, "/DESeq2/shrunk/", out.file.base, "_shrunk_", shrnk.type, "_ann_results.tab", sep='')
				cat("\tSaving", out.file, '\n')
				write.table(ann.shrunk, 
					        file=out.file,
					        row.names=T, col.names=T, sep='\t')
			}
		# Else, unnannotated shrunk.lfc will be represented in 
		} else ann.shrunk <- shrunk.lfc
		
		if (FALSE) {
		    if (save) {
		        out.file <- paste(out.dir, "/DESeq2/img/DESeq2_", out.file.base, "_shrunk_", shrnk.type, "_l2FC.png", sep='')
		        cat("\n\tPlotting", out.file, '\n\n')
		    } else out.file <- NULL
		    
		    as.png( {
					ggplot2::ggplot(data.frame(ann.shrunk), 
					ggplot2::aes(x = log2(baseMean), y=log2FoldChange),
					environment=environment()) + 	# this is supposed to make it work in a local env
					
					ggplot2::geom_point(ggplot2::aes(colour = padj < alpha), shape = 10, size = 1) +
					
					#ggplot2::geom_text(data=~top_n(.x, 10, wt=-padj),
					ggplot2::geom_text(data = top_n(ann.shrunk, 10, "padj"),
					ggplot2::aes(label = rownames(ann.shrunk))) +
					ggplot2::labs(x="Log2 Mean of Normalised Counts", y="Log2 Fold Change")
		    
		    }, out.file , overwrite = TRUE)
		}
		
		# Find statistically significant changes. Ignore genes with NAN p-value.
		signif <- results[ (results$padj < alpha) & ! is.na(results$padj) , ]
		signif$abs_lfc <- abs(signif$log2FoldChange)
		
		if (save) {
		    out.file <- paste(out.dir, '/DESeq2/signif/signif_', out.file.base, "_α<", alpha, ".tab", sep='')
		    cat("\n\tSaving", out.file, '\n')
		    write.table(data.frame(signif), file = out.file, sep = '\t', row.names = T, col.names = T)
		}
		
		# Order the data, firstly by the decreasing absolute 
		# value of log2FoldChange and secondly by the increasing pvalue....
		srt.signif <- signif[ order( signif$abs_lfc,
								signif$padj,
								decreasing=c(T, F)), ]
		if (save) {
		    out.file <- paste(out.dir, "/DESeq2/signif/signif_sorted_", out.file.base, ".tab", sep='')
		    cat("\tSaving", out.file,'\n\n')
		    write.table(srt.signif, file = out.file, sep = '\t', row.names = T, col.names = T)
		}
		
		## If shrunken data is unnanotated, return null.
		if ( ! annotate) ann.shrunk <- NULL
		
		# Generate output list containing results    
		res.list <- list(	result.ann = results, 
		            		shrunk = shrunk.lfc, 
		            		shrunk.ann = ann.shrunk,
							signif = srt.signif[,1:7],
		            		signif.ann = srt.signif)
		
		# Save results as an object to avoid redundant analyses.
		cat("\n\tSaving final results ...\n\n")		
		saveRDS(res.list, file=paste(out.dir, "/DESeq2/", out.file.base, "_final_result_list.rds", sep=''))
		save(res.list, file=paste(out.dir, "/DESeq2/", out.file.base, "_final_result_list.RData", sep=''))
		
		cat("\n =========================== \n")
		cat(  "| A N A L Y S I S   D O N E |\n")
		cat(  " =========================== \n")
        
	} else {
    
		cat("\n\t Loading from file ", basename(out.file.base), "_final_result_list.rds", "\n\n", sep = "")
    	res.list <- readRDS(paste(out.dir, "/DESeq2/", out.file.base, "_final_result_list.rds", sep = ""))
		#load(paste(out.dir, "/DESeq2/", out.file.base, "_final_result_list.RData", sep = ""))
		cat("\n\t DONE! \n\n")
        }
        return(res.list)
}
