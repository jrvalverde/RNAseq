library(DESeq2)
library(tibble)			# general tools
library(tidyr)
library(stringr)
library(dplyr)
library(readr)

for (column in c("ko", "persistent", "sample", "status", "type", "wt") ) {
    rnaseq.out <- paste("rnaseq.by", column, "cnb", sep='.')
    folder <- paste(rnaseq.out, "any_end", sep='/')
    cat("Processing", rnaseq.out, "\n")

    countDataFile <- paste(folder, "/feature_counts/featureCounts.csv", sep='')

    countData <- read.csv(countDataFile, row.names=1)
	countData <- as.matrix(countData)
    countData <- countData[rowSums(countData)>1, ]
    head(countData)

    colData <- read.delim(paste(rnaseq.out, "SampleInfo.tab", sep='/'), stringsAsFactors=TRUE)


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
    #
    metadata <- colData

    counts <- countData

    dds.col <- DESeqDataSetFromMatrix(counts, 
                metadata,  
                design=as.formula(paste("~", column)))
                
    dds.col <- estimateSizeFactors(dds.col)	# add norm factors to dds.col
    sizeFactors(dds.col)
    normalized_counts <- counts(dds.col, normalized=TRUE)	# get normalized counts

    norm_counts_file <- paste(folder, "/feature_counts/normalizedCounts_by_", column, ".csv", sep='')
    cat("Writing", norm_counts_file, "\n")
    write.table(normalized_counts, 
                file=norm_counts_file,
	            sep='\t', row.names=T, col.names=T) # better if row.names=F

}
