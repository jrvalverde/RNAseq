### MOST RECENT

#
# We start from data that has already been cleaned by performing a quality
# check with FastQC and subsequent edge trimming.
#
# The next step is to align the reads and calculate the counts of reads
# that map to each gene. This has been done as well previously, obtaining
# alignments in BAM format: the .bam files (binary files containing the
# reads aligned to the reference genome sequence) and the .bam.bai files
# containing the indexes for each one). This step is mostly a matter of
# time and has already been done:
# 
#   Paired-end Illumina short-reads were aligned against Coturnix japonica
#   genome (v2.0 primary assembly) using RNA-STAR (1) (--outReadsUnmapped
#   Fastx; --alignIntronMax 10000; -- alignMatesGapMax 10000). PCR and optical
#   duplicates were marked using the MarkDuplicates function of Picard-Tools 
#   (GATK)(2) (TAGGING_POLICY=ALL). Alignment results, saved as BAM files, 
#   were sorted and indexed using samtools (3).
# 
# We take the alignments produced by the Bioinformatics for Proteomics and
# Genomics Service of CNB directly, which are stored in the folder
# 'Alignments_Coturnix', and the reference genome used by them which is in
# the folder 'refGenomes/Coturnix_Japonica'. The reference files correspond
# to the 2.0 primary assembly. This is important for we will need their
# indexes for the next step.
#
# Now we need to calculate the counts of the reads that map to gene regions,
# summing up the reads that match each gene. For this, we use each .bam file
# and process it with the function featureCounts from R package Rsubread, and
# the reference genome data, which is in file Cjaponica.gtf. This last file
# contains the information about the features annotated in the genome of
# Coturnix japonica, including the coordinates of each feature in the reference
# genome. The function featureCounts will use these coordinates to know to
# which feature each read maps.
#
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
#reference <- 'Gg_GRGc7b'
#release <- "GRCg6a"
#target.organism <- 'Gallus_gallus'
#ncbi.taxid <- '9031'
#genus <- 'Gallus'
#species <- 'gallus'
#version <- '6a'
#alignment.dir <- 'gg-aln-Rsubread-6a'
#ens.db.pkg <- "EnsDb.Ggallus"
#ens.version <- '106'
#mart.name <- 'ggallus_gene_ensembl'
#rnaseq.out <- 'rnaseq-6a'
#org.package <- "org.Ggallus.eg.db"
#org.package <- "org.Gg.eg.db"
#KEGG_org <- 'gga'


n.genes <- 1000		# number of top genes to revise

fastq.data <- 'fastq-qc-paired'
#data <- fastq.data

out <-  alignment.dir

if (BOTH == T) {
    requireBothEnds <- T
    folder <- paste(rnaseq.out, "both_ends", sep='/')	# whether both ends sshould match or just any end
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


openlog <- function(name=logfile) {
    i <- 0
    name <- sprintf("%s.%03d", logfile, i)
    # check the existence of increasing numbers starting at 0 of log file versions
    while (file.exists(name)) {
        i <- i + 1
        name <- sprintf("%s.%03d", logfile, i)
    }
    # when we find a free number, use it
    # split=T means that the output will go to the screen and to a file
    sink(name, split=T, append=T)
}

closelog <- function() {
    # close ALL sinks. CAUTION CAUTION CAUTION
    while (sink.number()) { sink() }
}

#library(keypress)
continue.on.key <- function() {
    if (! INTERACTIVE)
        return('')

    cat("\nPress any key to continue: ")
    key <- keypress()
    return(key)
}

continue.on.enter <- function(prompt='Press ENTER to continue: ') {
    if (! INTERACTIVE)
        return('')
    
    return(readline(prompt))
}

#library(keypress)
more.columns <- function (data, columns=c(1:dim(data)[2]), lines=20, header='', prompt="More? ") {
    maxln <- dim(data)[1]
#print(dim(data))
    if (maxln < lines) lines <- maxln
    cont <- TRUE
    st <- 1; en <- lines
#print(paste("start",st, "end", en, sep=' '))
    ans <- ''
    while (cont) {
#print(paste(">>> start",st, "end", en, sep=' '))
        #cat(" ", data[ st:en, columns], '\n', sep='')
        #print(data[ st:en, columns])
        if (header != '') cat(header, '\n', sep='')
        for (i in st:en) {
            cat("<", i, "> ", sep="")
            for (c in columns) {
                cat('\t"', data[i, c], '"', sep='')
            }
            cat('\n')  
        }
        cat(prompt)
        ans <- keypress()
#print(paste("|", ans, "|", sep=''))
        if (ans == ' ') {
            st <- st + lines
            en <- en + lines
#print(paste("start",st, "end", en, sep=' '))
        } else if (ans == 'b' | ans == 'B') {
            st <- st - lines
            en <- en - lines
        } else if (ans == '>') {
            en <- maxln
            st <- en - (lines - 1)
        } else if (ans == '<') {
            st <- 1
            en <- lines
        } else if (ans == '' | ans == 'down') {
            st <- st + 1
            en <- en + 1
        } else if (ans == 'up') {
            st <- st - 1
            en <- en - 1
        } else if (ans == 'q' | ans == 'Q') {
            cont <- FALSE
        } else cont <- FALSE
        # check limits
        if (en > maxln) {
            en <- maxln
            st <- en - (lines - 1)
        }
        if (st < 1) {
            st <- 1
            en <- lines
        }
        cat('\n')
#print(paste("start",st, "end", en, sep=' '))
    }
    return (ans)
}


show.data.frame <- function(df, window.name='data.frame', visible=TRUE) {
    #library(tcltk)
    #library(gWidgets2)
    window.name <- window.name
    window <- gwindow(title=window.name, visible=FALSE)
    tab <- gtable(df,
           container=window)
    if (INTERACTIVE)
        visible(window) <- TRUE	# setting it to FALSE removes window
    return(window)
}



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
    if ( ! file.exists(paste(save,dir, 'bam_files_stats.txt', sep='/'))) {
        bam.files <- list.files(path=aln.out, pattern='.BAM$', full.names = TRUE)
        bam.files
        props <- propmapped(files=bam.files)
        props
        write.table(props, file=paste(save,dir, 'bam_files_stats.txt', sep='/'), 
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



h.cluster.changes <- function(data.table, annotated.data, 
                              distance="euclidean",	# see dist()
                              clusters=1,
                              method="complete",	# see hclust()
                              estimate=T
                              ) {
                              

    # Normalize the data
    by.row <- 1
    by.col <- 2
    dif <- data.table
    means <- apply(dif, by.col, mean)
    sds <- apply(dif, by.col, sd)
    nor <- scale(dif,center=means,scale=sds)

    # Calculate distance matrix  
    distan = dist(nor, method=distance)

    # we will try several clustering strategies

    # Hierarchical agglomerative clustering  
    cat("
    H I E R A R C H I C A L   C L U S T E R I N G
    =============================================
    \n\n")
    
    # we should consider using hcut() instead which cuts and allows fviz use
    
    hc = hclust(distan)
#    plot(mydata.hclust)
#    plot(mydata.hclust,hang=-1)
    print(plot(hc,labels=rownames(data.table),main='Default from hclust'))
    continue.on.enter("Press [ENTER] to continue ")
    
    # Cluster membership
    if ((clusters <= 1) & (estimate == T)) {
        ans <- continue.on.enter(prompt="How many clusters should we use? ")
        nclust <- as.numeric(ans)
    } else
        nclust <- clusters
    
    member <- cutree(hc, nclust)	# cut to 4 groups
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
        cat("Showing annotation for cluster", i, "(", length(clus.i), " elements)\n")
        
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
        window <- show.data.frame(data.to.show, window.name)
        visible(window) <- F
        #visible(window) <- TRUE	# setting it to FALSE removes window
#       keypress()
    }

    cat("
    
    Silhouette Plot for hierarchical clustering with normalized data
    ----------------------------------------------------------------
    Measure similarity of each object to its own cluster (cohesion) compared to
    other clusters (dispersion). Values range from -1 to +1. Large values
    indicate objects well matched to their own cluster and badly to neighboring
    clusters. If many points have low or negative value, the number of clusters
    is too low or too high.
    \n")
    print(plot(silhouette(cutree(hc, nclust), distan)))

    return(hc)
}



hcut.cluster.changes <- function(data.table, annotated.data,
                                    algorithm="hclust",	# see hcut()
                                    distance="euclidean", # see hcut()
                                    method="ward.D2",	# see hcut()
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

    cat("
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
    print(fviz_cluster(hc, data=nor, geom="point"))
    continue.on.enter("Press [ENTER] to continue ")
    # to plot gene names instead of ensembl.ids we need to change rownames to genes
    gor <- nor 
    rownames(gor) <- paste(annotated.data[ , "GENENAME"], 1:length(rownames(gor)))
    print(fviz_cluster(hc, data=gor))
    continue.on.enter("Press [ENTER] to continue ")
}



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

    cat("
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



pam.cluster.changes <- function(data.table, annotated.data,
                                    distance="euclidean", # see dist()
                                    metric="euclidean",	# see pam()
                                    clusters=1,	# n. of clusters
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
    Partition aroung medoids
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
    ### JR ### NOTE: may be worth trying to cluster separately with diss=TRUE)
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
        use_dist_matrix <- FALSE	### JR ### fviz_cluster fails, why?
        if (use_dist_matrix == TRUE) {
            ### JR ### NOTE: may be worth trying to cluster with diss=TRUE)
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


##############################################################################
##############################################################################
##############################################################################

							S T A R T

##############################################################################
##############################################################################
##############################################################################




openlog(logfile)
# when a script ends or when stop() is called, there are a number
# of housekeeping tasks to do. on.exit() allows us to add additional
# tasks so that they, too, are executed at the end of the script
# this way we need not keep track of all the sink()s
on.exit(closelog, add=T)





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



###################################################
# SELECT ENSEMBL ANNOTATION DATABASE TO USE
###################################################

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

# check it
print(ens.db)
head(keys(ens.db, 'GENEID'))
columns(ens.db)


###################################################
###   B I O M A R T  A N N O T A T I O N   ###
###################################################

# ---------------------------------------------------------------
# M A K E     O R G     P A C K A G E
# ---------------------------------------------------------------

org.db <- NULL

# use AnnotationHub to seek a suitable package
qo <- query(ah, "Orgdb")
if ( last.ref %in% names(qo) ) {
    org.db <- query(ah, "Orgdb")[[ last.ref ]]
} 

if (is.null(org.db)) {
	library(org.package)
	# org.package has already been loaded
	# convert the packahe name to variable name and save the variable
	org.db <- get(org.package)
}

# at this point ens.db should contain the ENSEMBL data and 
# *.org.db the Org type data.





# ---------------------------------------------------------------
# O B T A I N   B I O M A R T   A N N O T A T I O N
# ---------------------------------------------------------------
if ( file.exists(paste(folder, 'biomaRt.annotation.1st.txt', sep='/')) ) 
    bm.annot.1 <- read.table(
    		file=paste(folder, '/biomaRt.annotation.1st.txt', sep=''), 
	        sep='\t', 
                header=T)
} else {
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
                                    "entrezgene_id", 
									"entrezgene_accession", 
									"entrezgene_description"),
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
                                    "go_id", 
									"name_1006", 
									"definition_1006", 
									"go_linkage_type", 
									"namespace_1003"), 
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
                                    "goslim_goa_accession", 
									"goslim_goa_description"),
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
    #out <- 'gg-aln-Rsubread-6a'
    out <- 'Alignments_Coturnix'
    align.fastq(path=path, reference=reference, aln.out=out, save.dir=folder) 
} else {
    #out <- 'gg-aln-Rsubread-6a'
    out <- 'Alignments_Coturnix'
}

bam.files <- list.files(path = out, pattern = '.bam$', full.names = TRUE)

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

# Now we have the feature counts and we can do differential expression


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
# A N A L Y S I S     W I T H     E D G E R
# ---------------------------------------------------------------


# use convenience names
counts <- fc$counts
annotation <- fc$annotation
stats <- fc$stats
#target <- read.delim("data/SampleInfo.txt")
#target <- read.delim("data/SampleInfo.tab")
target <- read.delim(paste(rnaseq.out, 'SampleInfo.tab', sep='/'))
print(target)

# then we proceed to the analysis.. for this, we will need other packages:
#library(edgeR)
#library(limma)
#library(RColorBrewer)
#library(gplots)

# The next steps are to distinguish the genes whose expression is significant
# from  the ones that have an 'insignificant' expression that could be by
# chance or irrelevant to the  situation of the cell. So we set up the
# threshold of expression for each gene to a minimum 10-15 counts. Since
# we will work with normalized CPM (counts per million) data, we need
# to check which number of CPM corresponds to  to this amount of counts in
# order to filter the data

# convert counts to CPM
myCPM <- cpm(counts)


# we'll plot the correlation of cpm and counts to see which is the number of
# cpm that corresponds to 10-15 counts minimum. we see this information
# graphically we use for example column 1... we could check every file
# independently but since they are all similar in number of reads there 
# should be no need

out.png <- paste(folder, '/edgeR/img/edgeR_CPMcounts.png', sep='')
as.png( {
        plot(counts[,1:15], myCPM[,1:15], xlim=c(0,20), ylim=c(0,5))
        abline(v=10, col=1)
    }, out.png )

# we will set up the threshold to 0.5 according on the plot we have drawn
# this command will return a table of trues and falses, then, we want to keep
# only the rows that exceed the threshold at least in three different
# cases(experiments .bam files)
thres <- myCPM > 0.5
keep <- rowSums(thres) >=3
if (INTERACTIVE)
    table(keep)

# Note that we could as well have used the counts directly with
# keep <- rowSums(cnt > 10) >= 3
# table(keep)
# which would be less approximate and yield a different selection


# then we store in a different variable the genes whose counts exceed the
# threshold and visualise the content of the new variable to see the amount of
# remaining genes
counts.keep <- counts[keep,]
dim(counts.keep)


# We have manipulated the data discarding whatever is not of high interest
# and now we need to see the differencial expression and highlight the
# differences among the cells

# convert the counts.keep to a DGEList
dge.keep <- DGEList(counts.keep)


# do TMM normalization
dge.norm <- calcNormFactors(dge.keep)
dhe <- dge.norm
dge$samples


# Now, do some quality control plots, barplots and boxplots
# we need normalized data counts so we take the logarithm
# we need to group together all the experiment data that correspond to the
# same  cell type (target) and asign a different color to each one so as to
# distinguish them graphically; we know we have 4 groups so we have to specify
# 4 different colors
out.png <- paste(folder, '/edgeR/img/edgeR_sample_lib_size.png', sep='')
as.png(barplot(dge$samples$lib.size), out.png)

logcpm <- cpm(dge$counts, log=TRUE)
group.col <- c('red', 'blue', 'green', 'yellow')[target$PFU] 

out.png <- paste(folder, '/edgeR/img/edgeR_log2_cpm.png', sep='')
as.png( {
        par(mfrow=c(1,1))
        boxplot(logcpm, xlab='', ylab=' Log2 counts per million', 
	        col= group.col, las=2, outline=FALSE)
        abline(h=median(logcpm), col='red')
        title('Boxplots of logCPMs unnormalised')
    }, out.png )

# then we check the MDS plot to see any significant difference between the
# groups
out.png <- paste(folder, '/edgeR/img/edgeR_mds_plot.png', sep='')
as.png( {
        par(mfrow= c(1,1))
        plotMDS(dge, col=group.col)
    }, out.png)


# now is time to see the differences in the expression (variance). We have
# to apply a funcion that calculates the variance by rows (genes) and then
# retrieve the n.genes most DE genes

logcounts <- cpm(dge, log=TRUE)
var_genes <- apply(logcounts, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:n.genes]
highly_var <- logcounts[select_var,]
if (INTERACTIVE)
   dim(highly_var)

# we now set up the colours we will use for the heatmap plot 
# and then we create all the colors in between in the palette

mypalette <- brewer.pal(11, 'RdYlBu')
morecolors <- colorRampPalette(mypalette)

# we plot the heatmap.2 (gplots) without a line (trace), scale by row
# (difference in color) margins (something about the labels used), also we
# reverse the colors because by default the red is associated with low
# expression and we are not familiar with it
out.png <- paste(folder, '/edgeR/img/edgeR_heatmap.png', sep='')
as.png( {
        #margins <- par("mar")
        #par(mar=c(25, 5, 5, 10))
        heatmap.2(highly_var, 
                  col= rev(morecolors(50)), 
                  trace='none', 
                  ColSideColors=group.col,
                  scale='row', 
                  margins= c(15,5))
        #par(mar=margins)
    }, out.png)

# This plot is of limited use. We'd better have other names for rows and 
# columns and plot the n first (most variable) to see them well
n <- 50
high_var <- highly_var
colnames(high_var) <- gsub("_R1.bam", "", colnames(highly_var))
name <- ensembldb::select(ens.db, keys=rownames(high_var), 
                   column='GENEID', keytype='GENEID', 
                   columns=c('GENENAME'))
out.png <- paste(folder, '/edgeR/img/edgeR_heatmap.', n,'.png', sep='')
as.png( {
        #margins <- par("mar")
        #par(mar=c(25, 5, 5, 10))
        heatmap.2(high_var[1:n,1:15], 
                      col=rev(morecolors(50)), 
                      trace='none', 
                      ColSideColors=group.col,
                      scale='row', 
                      margins= c(15,5),
                      labRow=name[1:n, 2])
        #par(mar=margins)
    }, out.png )


# We can automate the estimation of the dispersion and add it to the dge object
dge <- estimateCommonDisp(dge)

# And now we can estimate gene-wise dispersion estimates allowing for a
# possible trend with average count size. These will allow us to use a GLM
# instead of a plain LM.
# This gives us the BCV (Biological Coefficient of Variation) between samples
dge <- estimateGLMTrendedDisp(dge)
dge <- estimateTagwiseDisp(dge)

out.png <- paste(folder, '/edgeR/img/edgeR_BCV_dispersions.png', sep='')
as.png(plotBCV(dge), out.png)

############## JR #########################
# 
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
# edgeR User’s Guide
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
############## JR #########################

# we know that the variability we see in the expression depends on infection so
# we have to take this into account and create a model
# 0+ forces the design to include all groups and not have
# an intercep (reference) column
design.column <- "viral.dose"

design <- model.matrix(~0 + target[ , design.column])	# no reference
#design
#design <- model.matrix(~ target[ , design.column])
colnames(design) <- levels(as.factor(target[ , design.column]))
design

# This is intended to be used by the Gallus gallus experiment
#design.persist.inf <- model.matrix(~ 0 + target$Persistent * target$Infected)
#design.persist.inf <- model.matrix(~ target$Persistent * target$Infected)
#design.persist.inf

# IF WE DO NOT INCLUDE THE "~ 0 + " IN THE FORMULA, MODEL.MATRIX()
# WILL USE AS REFERENCE THE FIRST ALPHABETICAL ORDER LEVEL !!!
# WHICH HERE WOULD BE 0.1 !!!
# Another problem is that we are limited to the comparisons defined by
# coefficents, if we want more control we need to use limma::makeContrasts
#
# If we use as formula "~ table[ , design.column ]" then we are limited to 
# comparisons defined by the fitting coefficients (ref vs. variable-in-coeff).
# If we use "~ 0 + table[ , design.column]" we would have diffs for
# all levels, but no reference.


# Let us test for differential expression with the DGE data using a GLM:
# first, fit genewise GLMs
gfit <- glmFit(dge, design)
if (INTERACTIVE) {
	names(gfit)
	head(coef(gfit))
}

fit <- gfit

# We can now conduct Likelihood Rato tests and show the top genes
# for the selected comparison coefficients (reference vs. coeff)
lrt <- glmLRT(gfit, coef=1)	# coef = 1... length(gfit$coefficients)
topTags(lrt)


# Before going forward, let us prepare the annotation we will use:
# we now want to extract some annotation 
# to add useful information to our genes
ens.ann <- ensembldb::select(ens.db, 
                  column='GENID', keytype= 'GENEID', keys=rownames(fit), 
                  columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION', # no longer available
                             'GENENAME', 'GENEID', 'ENTREZID', # empty
                             'TXID', 'TXBIOTYPE', # these make the call fail
                             'PROTEINID', 'UNIPROTID' # no longer available
                            ))



# SAVE ANNOTATION
# ---------------
write.table(ens.ann, file=paste(folder, '/ensembl.annotation.txt', sep=''), 
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


# let's do a VOOM analysis with the formula employed
# in the current design

DO_VOOM=TRUE
if (DO_VOOM) {
    #voom transformation of the data 
    v <- voom(dge, design, plot=FALSE)
    out.png <- paste(folder, '/edgeR/img/edgeR_voom.png', sep='')
    as.png( {
            par(mfrow= c(1,1))
            voom(dge, design, plot=TRUE)
        }, out.png )



    # Carry on a summary variation analysis using the VOOM transformed data
    #
    # we fit the results of the voom depending on the model we created with our
    # design, we overwrite the fit with the use of eBayes (stat model), and
    # depending on that we decide abut which tests fit the best and summarise the
    # results: topTable gives us helpful information about everything

    v.fit <- lmFit(v, design)
    v.fit <- eBayes(v.fit)
    results <- decideTests(v.fit)
    summary(results)
    #topTable(fit, coef= 1, sort.by='p')

    # SAVE TOP 'N' GENES
    # ------------------
    # save the table of the n.genes top expressed genes sorted in various orders
    eR.save.top <- function(fit, folder, n.genes=500, sort.by='p', coef=1, p.value=0.01) {
        # default adjustment is BH
        tt <- topTable(fit, coef=coef, number=n.genes, sort.by='p', p.value=p.value)
        write.table(tt, paste(folder, '/edgeR/cmp_coef=', coef, '_top_', n.genes, '_by', sort.by, '.tab', sep=''), 
	    row.names=T, col.names=T, sep='\t')

    }

    eR.save.top(v.fit, folder, n.genes, 'p')
    eR.save.top(v.fit, folder, n.genes, 'B')
    eR.save.top(v.fit, folder, n.genes, 'logFC')
    eR.save.top(v.fit, folder, n.genes, 'AveExpr')


    # the problem here is that we are limited to the comparisons defined by
    # coefficents, if we want more control we need to use makeContrasts


    #---------------------------------------------------------------------------

    #Now is time to connect all the results we have with the existing
    #information  we know from the literature, so we will retrieve infos from
    #ENSEMBL and connect with the genes  we have kept. We are interested only in
    #genes and transcriptomes (non-characterized genes) from the organism
    #Coturnix japonica. there is no database available for this organism so we
    #have to create it by ourselves

    # we created two different databases one from the gtf file (small archive)
    # and one from the data we extracted from ENSEMBL and stored it in a sqlite 
    # file

    # we'll use ens.db from above

    # As long as it works perfectly, we continue by connecting the 
    # annotation information to the fit using a new bucket named "genes"
    # (fit is an object of class "MArrayLM" (package "limma") which is a
    # list. I.e. we assign the annotation to a new element named 'genes'  
    # of this list.
    #	This works if they are in the same order
    #fit$genes <- ens.ann.1
    #	This checks they go in the same order)
    v.fit$genes <- ens.ann.1[ match(rownames(fit), ens.ann.1$GENEID), ]


    # SAVE THE ANNOTATED FIT
    # ----------------------
    # save all the contents of 'fit' in an RDS file
    saveRDS(v.fit, file=paste(folder, '/edgeR/annotatedVOOMfit.rds', sep=''))
    #	'fit' can later be recovered with: 
    #		fit <- readRDS(file=paste(folder, '/annotatedVOOMfit.rds', sep=''))
    # and save as well as Rdata file
    save(v.fit, file=paste(folder, '/edgeR/annotatedVOOMfit.RData', sep=''))
    #	'fit' can later be recovered with: 
    #		fc <- load(file=paste(folder, '/annotatedVOOMfit.RData', sep=''))


    # SAVE TOP 'N' GENES ANNOTATED
    # ----------------------------
    # we can now run again the topTable and we will have all the annotation 
    # information linked 
    #	n.genes <- 500 ALREADY DEFINED ABOVE
    eR.save.top.annotated <- function(fit, folder, n.genes=500, sort.by='p', coef=1, p.value=0.01) {
        # default adjustment is BH
        tt <- topTable(fit, coef=coef, number=n.genes, sort.by='p', p.value=p.value)
        write.table(tt, paste(folder, '/edgeR/cmp=', coef, '_top_', n.genes, '_annotated_by_', sort.by, '.tab', sep=''), 
	    row.names=T, col.names=T, sep='\t')

    }

    eR.save.top.annotated(v.fit, folder, n.genes, sort.by='p')
    eR.save.top.annotated(v.fit, folder, n.genes, sort.by='logFC')
    eR.save.top.annotated(v.fit, folder, n.genes, sort.by='AveExpr')


    #now create the volcano plot for only the top 1/2 genes
    out.png <- paste(folder, '/edgeR/img/edgeR_volcanoplot.png', sep='')
    as.png(volcanoplot(v.fit, highlight=n.genes/2, coef=1, names=fit$genes$SYMBOL),
        out.png)

    #Testing relative to a threshold: 1 means a 2x fold change
    threshold=1
    v.fit.thres <- treat(v.fit, lfc=threshold)
    res.thres <- decideTests(v.fit.thres)
    summary(res.thres)
    topTreat(v.fit.thres, coef=1, sort.by='p')

    # SAVE TOP RESULTS RELATIVE TO THRESHOLD
    # --------------------------------------
    #	n.genes <- 500 ALREADY DEFINED ABOVE
    eR.save.top.ann.thresh <- function(fit, folder, n.genes=500, sort.by='p', coef=1, p.value=0.01, threshold=1) {
        # default adjustment is BH
        tt <- topTreat(fit, coef=coef, number=n.genes, sort.by='p', p.value=p.value)
        write.table(tt, paste(folder, '/edgeR/cmp=', coef, '_top_', n.genes, '_annotated_lfc>=', threshold, '_by_', sort.by, '.tab', sep=''), 
	    row.names=T, col.names=T, sep='\t')

    }

    eR.save.top.ann.thresh(v.fit.thres, folder, n.genes, 'p')
    eR.save.top.ann.thresh(v.fit.thres, folder, n.genes, 'logFC')
    eR.save.top.ann.thresh(v.fit.thres, folder, n.genes, 'AveExpr')


    # and now we need to annotate 'fit$genes' with the biomaRt data
    v.fit$genes <- merge(v.fit$genes, bm.annot.1, by.x="GENEID", by.y="ensembl_gene_id")
    head(v.fit$genes)

    v.fit.thres$genes <- merge(v.fit.thres$genes, bm.annot.1, by.x="GENEID", by.y="ensembl_gene_id")
    head(v.fit.thres$genes)


    # save fit data with the extended version containing an extra
    # annotated bucket:
    # save all the contents of 'fit' in an RDS file
    saveRDS(v.fit, file=paste(folder, '/edgeR/annotatedVOOMfit+.rds', sep=''))
    #	'fit' can later be recovered with: 
    #		fit <- readRDS(file=paste(folder, '/annotatedVOOMfit+.rds', sep=''))
    # and save as well as Rdata file
    save(v.fit, v.fit.thres, file=paste(folder, '/edgeR/annotatedVOOMfit+.RData', sep=''))
    #	'fit' can later be recovered with: 
    #		fc <- load(file=paste(folder, '/annotatedVOOMfit+.RData', sep=''))

    # save top N genes, annotated-
    eR.save.top.bm.ann.thresh <- function(fit, folder, n.genes=500, sort.by='p', coef=1, p.value=0.01, threshold=1) {
        # default adjustment is BH
        tt <- topTable(fit, coef=coef, number=n.genes, sort.by='p', p.value=p.value)
        write.table(tt, paste(folder, '/edgeR/cmp=', coef, '_top_', n.genes, '_bm.annotated_lfc>=', threshold, '_by_', sort.by, '.tab', sep=''), 
	    row.names=T, col.names=T, sep='\t')

    }

    eR.save.top.bm.ann.thresh(v.fit.thres, folder, n.genes, 'p')
    eR.save.top.bm.ann.thresh(v.fit.thres, folder, n.genes, 'B')
    eR.save.top.bm.ann.thresh(v.fit.thres, folder, n.genes, 'logFC')
    eR.save.top.bm.ann.thresh(v.fit.thres, folder, n.genes, 'AveExpr')

}  # if (VOOM == TRUE)

# ---------------------------------------------------------------------
# Do all comparisons at once

fit <- v.fit

eR.annotate.fit <- function(fit, end.db, biomart) {
    ens.ann <- ensembldb::select(ens.db, 
                  column='GENEID', 
                  keytype= 'GENEID', keys=rownames(fit), 
                  columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION', # no longer available
                             'GENENAME', 'GENEID', 'ENTREZID', # empty
                             'TXNAME', 'TXID', 'TXBIOTYPE', # these make the call fail
                             'PROTEINID', 'UNIPROTID' # no longer available
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
    fit$genes <- ens.ann.1[ match(rownames(fit), ens.ann.1$GENEID), ]
    rownames(fit$genes) <- rownames(fit)
    #
    # add biomart annotation
    fit$genes <- merge(fit$genes, biomart, by.x="GENEID", by.y="ensembl_gene_id")
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
	row.names=F, col.names=T, sep='\t')
    # if fit is annotated the annotation will also be saved
}

# for Coturnix we use 'viral.dose' as basis for comparison
#design.column <- 'viral.dose'
# for Gallus we will use 'Src'
# For Coturnix we will use 'PFU' or 'viral.dose'
design.column <- 'viral.dose'
design <- model.matrix(~0 + target[ , design.column])
colnames(design) <- levels(as.factor(target[ , design.column]))
grps <- levels(as.factor(target[ , design.column]))
colnames(design) <- grps
n.grps <- length(grps)

qlfit <- glmQLFit(dge, design, robust=TRUE, abundance.trend=TRUE)
png.file <- paste(folder, '/edgeR/img/edgeR_QLFit_', design.column, '.png', sep='')
as.png(plotQLDisp(qlfit), png.file, overwrite=T)

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
        cat("\nComputing DGE:", i, '-', j, '\n')
        cmp <- makeContrasts(formula, levels=design)
        # glmQLFTest is similar to glmLRT except it uses Bayes quasi-likelihood
        # the P-values are always >= those produced by glmLRT
        qlf.cmp <- glmQLFTest(qlfit, contrast=cmp)
        # or if a threshold has been defined
        #qlf.cmp <- glmTreat(glfit, contrast=cmp, lfc=threshold)
        
        # get summary of up/down regulated genes
        if (VERBOSE)
			print(summary(decideTests(qlf.cmp)))
        png.file <- paste(folder, "/edgeR/img/edgeR_QLF_MD_", 
                          design.column, '_', i, '-', j, '.png', sep='')
        as.png( { 
                plotMD(qlf.cmp)
                abline(h=c(-1, 1), col="darkgreen")
                }, png.file)
        #continue.on.enter("Press [RETURN] to continue: ")
        
        # ANNOTATE the results
        qlf.cmp <- eR.annotate.fit(qlf.cmp, ens.db, bm.annot.1)

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
        qlf.result <- list(
	                  eR.cmp=qlf.cmp
                          #, eR.go=eR.go
			  #, eR.kegg=eR.kegg
			  )
        eR.data[[formula]] <- qlf.result 
        #print(topTags(qlf.cmp))
    }
}

# now eR.data is a list where each element is a comparison A-B (A minus B),
# i.e. each A-B is a list of tables, one of them named "table" and containing
# logFC, logCPM, F and PValue
#
# We would like to have FDR-corrected p-values as well, which can be got by
# defining a threshold. But that implies we know of a meaningful one, which
# we don't here.





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
		row.names=1) 
countData <- as.matrix(countData)
countData <- countData[rowSums(countData)>1, ]
head(countData)

colData <- read.delim(paste(rnaseq.out, "SampleInfo.tab", sep='/'))

colData$PFU <- as.factor(colData$PFU)
colData$Src <- as.factor(colData$Src)
colData$viral.dose <- as.factor(colData$viral.dose)

# SAVE NORMALIZED COUNTS FOR AI ANALYSIS
#
# We do not need normalized counts for differential gene expression
# because DESeq2 calculates them automatically.
#
# But we do need the normalized counts for the AI analysis, so we 
# will get and save them.
#
# edgeR-type TMM normalization depends on the comparison design
# try both and see if there are differences
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
# and we should make sure that in the design, wt is the first element in
# countData.

# We need to have metadata and feature counts in the appropriate
# name order for DESeq2 to work correctly.
# Here we reorder columns of metadata and counts to be in the proper
# order
metadata <- rbind(colData[13:15,], colData[1:3,], colData[10:12,], colData[7:9,], colData[4:6,])
rownames(metadata) <- 1:15

# sort feature counts from countdata
counts <- cbind(countData[,13:15], countData[,1:3], countData[,10:12], countData[,7:9], countData[,4:6])

dds.pfu <- DESeqDataSetFromMatrix(counts, metadata,  design=~PFU)
dds.pfu <- estimateSizeFactors(dds.pfu)	# add norm factors to dds.pfu
sizeFactors(dds.pfu)
normalized_counts <- counts(dds.pfu, normalized=TRUE)	# get normalized counts
write.table(normalized_counts, file="normalized_counts.tab",
	        sep='\t', row.names=T, col.names=T) # better if row.names=F



# RUN DGEseq2 DGE analysis

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

# use ens.ann from above

ens.ann.1 <- ens.ann[ ! duplicated(ens.ann$GENEID), ]

# Add annotation to dds to keep everything in one place
#dds$ens.annot.1 <- ens.ann.1
#dds$bm.annot.1 <- b,.annot.1



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

dds_compare_select_plot_annotate_save <- function(dds, 
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
        res <- dds_compare_select_plot_annotate_save( 
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


if (INTERACTIVE) {
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
}



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


# This should go up above all 
references <- c('wt')		# Cj
#references <- c('wt', 'PC')	# Gg
#contrasts.column <- "Src"
contrasts.column <- "PFU"

for (ref in references) {
    # Find genes that change w.r.t. the reference strain

    # create a convenience text variable to simplify/unify filenames below
    ccol_ref <- paste(contrasts.column, ref, sep='_')

    name <- paste(contrasts.column, ref, levels(as.factor(target[ , contrasts.column]))[1], sep='_')
    common <- rownames(ds.data[[name]]$signif)
    for (i in levels(as.factor(target[ , contrasts.column]))) {
        if (i == ref) next
        name <- paste(ccol_ref, i, sep='_')
        print(name)
        common <- intersect(common, rownames(ds.data[[name]]$signif))
    }
    length(common)	# 681 in Coturnix, 1670 in Gallus

    data.table <- data.frame(genes=common)
    #for (i in c("wt", "PC", "P", "1.0", "0.1")) {
    # compare ref (wt) against the different PFU-infected samples
    for (i in levels(as.factor(target[ , contrasts.column]))) {
    # we do not want to include PC this time
    #for (i in c("P", "1.0", "0.1")) {
        if (i == ref) next
        name <- paste(ccol_ref, i, sep='_')
        print(name)
        data.table[name] <- ds.data[[name]]$signif[common, "log2FoldChange"]
    }
    rownames(data.table) <- common

    data.annot <- ds.data[[name]]$signif.annot[common, ]

    par(mfrow=c(1,1))

    set.seed(1963)

    # First have a general look at the methods to get a feeling for the
    # best number of clusters


    dif <- data.table[ , -1]
    by.row <- 1
    by.col <- 2
    means <- apply(dif, by.col, mean)
    sds <- apply(dif, by.col, sd)
    nor <- scale(dif,center=means,scale=sds)
    # Do a scatterplot matrix
    car::scatterplotMatrix(dif)
    out.png <- paste(folder, "/DESeq2/img/", 
               'scatterplot_matrix_', ccol_ref, ".png", 
               sep='')
    as.png( {
            print(car::scatterplotMatrix(dif))
            }, out.png)


    # Try to guess the optimum number of K-means clusters
    # NBClust
    out.log <- paste(folder, "/DESeq2/cluster/NBClust_", ccol_ref, ".log", sep='')
    sink(out.log, split=T)
    library("NbClust")
    # predict best number of clusters for hierarchical and k-means clustering
    nbc <- NbClust(dif, diss=NULL, 
            distance="euclidean", method="complete", 
            min.nc=3, max.nc=10, 
            index="all", alphaBeale=0.1)
    # 4 for the C. japonica wt vs others data
    # 5 for G. gallus wt vs infected.data
    sink()


    # Let the user see the various clusters and make a decision
    #
    out.log <- paste(folder, "/DESeq2/cluster/kmeans/DESeq2_kmeans_all_", 
                     ccol_ref, ".log", sep='')
    sink(out.log, split=T)
    for (i in 2:10) {
        cat("Clustering with K-means (", i, " clusters)\n")
        cl <- kmeans(nor, i)				# NOTE: nor
        print(table(cl$cluster))
        print(fviz_cluster(cl, geom = "point", data=nor))
        out.png <- sprintf("%s/DESeq2/cluster/kmeans/DESeq2_kmeans_%s_nc=%03d.png", 
        	folder, ccol_ref, i)
        as.png(fviz_cluster(cl, geom = "point", data=nor,
                    main=paste("K-means clsutering nc =", i) ),
                out.png)
        #continue.on.enter("Press [ENTER] to continue ")
    }
    sink()
    
    cat("The log and plots for K-means clustering have been saved to\n", 
        "\t", folder, "/DESeq2/cluster/kmeans/\n", 
        "please, inspect them and select the best number of clusters\n",
        sep='')
    
    #kmeans.nc <- readline("How many clusters should I use for k-means? ")
    #kmeans.nc <- as.numeric(kmeans.nc)
    # for gg
    #kmeans.nc <- 5
    # for Cj
    kmeans.nc <- 4


    out.log <- paste(folder, "/DESeq2/cluster/pam/DESeq2_pam_all_", 
                     ccol_ref, ".log", sep='')
    sink(out.log, split=T)
    for (i in 2:10) {
        cat("\nClustering with PAM (", i, " clusters )\n")
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
    pam.nc <- 5

    out.log <- paste(folder, "/DESeq2/cluster/dbscan/DESeq2_dbscan_all_", 
                     ccol_ref, ".log", sep='')
    sink(out.log, split=T)
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
    # for Cj
    dbscan.eps <- 
	# 0.8


    hclust.nc <- 4	# we'll set it by hand for now

    # Now, proceed in detail, method by method, with a deeper and more 
    # detailed analysis

    # possible distances are c("euclidean", "maximum", "manhattan", "canberra",
    #	"binary", "minkowski", "pearson", "spearman", "kendall")
    # possible methods are c("ward.D", "ward.D2", "single", "complete",
    #	"average" (UPGMA), "mcquitty" (WPGMA), "median" (WPGMC), "centroid" (UPGMC))
    out.log <- paste(folder, "/DESeq2/cluster/hcluster/DESeq2_hcluster_", contrasts.column, "_", ref, ".txt", sep='')
    sink(out.log, split=T)
    h.cluster.changes(
        #data.table=pca.table[ , -1],
        #annotated.data=ds.data[[name]]$signif.annot[common, ],
        clusters=hclust.nc,
        data.table=dif,
        annotated.data=data.annot
    )
    sink()


    boot <- 100
 
    hclust_folder <- paste(folder, "/DESeq2/cluster/hclust_", ccol_ref, sep='')
    dir.create(hclust_folder, showWarnings=FALSE)
    
    hcut_folder <- paste(folder, "/DESeq2/cluster/hcut_", ccol_ref, sep='')
    dir.create(hcut_folder, showWarnings=FALSE)
    
    out.log <- paste(hcut_folder, "/DESeq2_hcut_", 
                     ccol_ref, ".txt", sep='')
    sink(out.log, split=T)
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


    kmeans_folder <- paste(folder, "/DESeq2/cluster/kmeans_", ccol_ref, sep='')
    dir.create(kmeans_folder, showWarnings=FALSE)

    out.log <- paste(kmeans_folder, "/DESeq2_kmeans_", 
                     ccol_ref, ".txt", sep='')
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


    pam_folder <- paste(folder, "/DESeq2/cluster/pam_", ccol_ref, sep='')
    dir.create(pam_folder, showWarnings=FALSE)

    out.log <- paste(pam_folder, "/DESeq2_pam_", 
                     ccol_ref, ".txt", sep='')
    sink(out.log, split=T)
    pam.cl <- cluster.changes(
        #data.table=pca.table[ , -1],
        #annotated.data=ds.data[[name]]$signif.annot[common, ],
        data.table=dif,
        annotated.data=data.annot,
        FUN=pam,
        distance="euclidean", # see dist()
        clusters=pam.nc,	# n. of clusters
        nstart=pam.nc*10,	# n. of random start sets to choose
        estimate=T,
        gap_bootstrap=boot,
        output.folder=pam_folder
        )
    sink()


    dbscan_folder <- paste(folder, "/DESeq2/cluster/dbscan_", ccol_ref, sep='')
    dir.create(dbscan_folder, showWarnings=FALSE)

    out.log <- paste(dbscan_folder, "/DESeq2_dbscan_", 
                     ccol_ref, ".txt", sep='')
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



    # Now cluster by experiment
    # -------------------------
	
    fid <- t(dif)
    means <- apply(fid, by.col, mean)
    sds <- apply(fid, by.col, sd)
    ron <- scale(fid,center=means,scale=sds)

    # here we have a small number of rows and can set a maximum number of clusters
    maxclust <- nrow(fid) - 1
    
    # hclust
    distan = dist(ron, method="euclidean")
    hcl <- hclust(distan)
    plot(hcl,labels=rownames(fid),main='Default from hclust')
    out.png <- sprintf(
            "%s/DESeq2_hcluster_grps_%s.png",
	    hclust_folder, ccol_ref)
    as.png(plot(hcl,labels=rownames(fid),main='Default from hclust'), out.png)

    # kmeans
    fviz_nbclust(ron, kmeans, method="silhouette", k.max=maxclust)
    out.png <- sprintf(
            "%s/DESeq2_kmeans_grps_silhouette_%s.png",
	    kmeans_folder, ccol_ref)
    as.png(fviz_nbclust(ron, kmeans, method="silhouette", k.max=maxclust), out.png)
    kcl <- kmeans(ron, centers=3, nstart=100)
    fviz_cluster(kcl, data=ron)
    out.png <- sprintf(
            "%s/DESeq2_kmeans_grps_%s_nc=%03d.png",
	    kmeans_folder, ccol_ref, 3)
    as.png(fviz_cluster(kcl, data=ron), out.png)

    # pam
    fviz_nbclust(ron, pam, method="silhouette", k.max=maxclust)
    out.png <- sprintf(
            "%s/DESeq2_pam_grps_silhouette_%s.png",
	    pam_folder, ccol_ref)
    as.png(fviz_nbclust(ron, pam, method="silhouette", k.max=maxclust), out.png)
    pcl <- pam(ron, k=3, diss=F)
    fviz_cluster(pcl, data=ron)
    out.png <- sprintf(
            "%s/DESeq2_pam_grps_%s_nc=%03d.png",
	    pam_folder, ccol_ref, 3)
    as.png(fviz_cluster(pcl, data=ron) , out.png)

    # dbscan
    # this results in the same two PCs but one cluster
    fviz_nbclust(ron, dbscan, method="silhouette", k.max=maxclust)
    out.png <- sprintf(
            "%s/DESeq2_dbscan_grps_silhouette_%s.png",
	    dbscan_folder, ccol_ref)
    as.png(fviz_nbclust(ron, dbscan, method="silhouette", k.max=maxclust), out.png)
    ### JR ### eps should be tuned for each experiment
    dcl <- dbscan(ron, eps=0.2, MinPts=2, showplot=1)
    fviz_cluster(dcl, data=ron)
    out.png <- sprintf(
            "%s/DESeq2_dbscan_grps_%s_eps=%03.2f.png",
	    dbscan_folder, ccol_ref, 0.2)
    as.png(fviz_cluster(dcl, data=ron), out.png)

    # this fails
    tryCatch( {
        NbClust(ron, diss=NULL, 
                distance="euclidean", method="complete", 
                min.nc=3, max.nc=10, 
                index="all", alphaBeale=0.1)
    } )
}

while (sink.number() > 0) sink()

