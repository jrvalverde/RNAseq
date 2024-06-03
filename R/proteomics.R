# load packages
library("MSnbase")
library(MSnSet.utils)
library("RColorBrewer") ## Color palettes
library("ggplot2")  ## Convenient and nice plotting
library("reshape2") ## Flexibly reshape data
library("vsn")
library("MALDIquant")
library("MALDIquantForeign")
library("BRAIN")
library("Rdisop")
library("OrgMassSpecR")
library("isobar")
library("DEP")

## Install missing packages
cran_packages <- c("remotes", "dplyr", "ggplot2")
for (pkg_i in cran_packages) {
  if (!require(pkg_i, quietly = T, character.only = T))
    install.packages(pkg_i)
}
if (!require("MSnSet.utils", quietly = T))
  remotes::install_github("PNNL-Comp-Mass-Spec/MSnSet.utils")
## ------------------------
library(MSnSet.utils)
library(dplyr)
library(ggplot2)

library("msmsTests")

if (!require("synapter")) {
    install.packages("BiocManager")
    BiocManager::install("synapter")
}

library("synapter")
library("qvalue")



# for multifile
#abundanceFile <- "Abundances_Scaled.tab"
#abundanceFile <- "Abundances_Normalized.tab"

# for single file
###nexprColsPrefix="SclAbun"
exprColsPrefix="NormAbun"
#exprColsPrefix="No.PSMs"

log2cutoff <- 1
p.cutoff <- 0.05




pos.left <- 1
pos.above <- 2
pos.below <- 3
pos.right <- 4
out.text <- function( labels, ... ) {
    n.lines <- length(labels)
    plot.new()
    text(x=0, y=seq(1, 1-(0.01*(n.lines-1)), -0.01), 
         labels=labels, 
         pos=pos.right, 
         ...)
}

out <- paste('Analysis', exprColsPrefix, 'pdf', sep='.')
pdf(out, paper="a4")


#
#
# https://www.bioconductor.org/packages/release/data/experiment/vignettes/RforProteomics/inst/doc/RforProteomics.html#4_Quantitative_proteomics
#
#

# read in the data
# using separate files
# (left here for reference just in case we need it in the future)
#abun <- readMSnSet(exprsFile=abundanceFile, 
#		  featureDataFile="Features.tab", 
#                  phenoDataFile="Phenotypes.tab", 
#                  sep='\t' )
                  
# using a single file (plus a metadata file)
# data.tab should have been saved as a TAB file with " string delimiter
#
# 1. get the metadata
phenoData <- read.table('metadata.tab', header=T)
# set column names to use for counts
exprCols <- sub("Count", exprColsPrefix, phenoData$label)
phenoData$label <- sub("Count", exprColsPrefix, exprCols)
rownames(phenoData) <- phenoData$label		#rownames(phenoData) <- rownames(pData(abun))

# 2. get all the data
fnames="accession"

abun <- readMSnSet2("data.tab", ecol=exprCols, fnames=fnames, sep='\t')
head(abun)
head(exprs(abun))
head(fData(abun))
head(pData(abun))	# not present in the data file, we need to fill it

# 3. populate pData
#pData(abun) <- cbind(pData(abun), read.table('metadata.tab', header=T)
#	or the equivalent, safer and more general
pData(abun) <- phenoData
pData(abun)

# now we have the data in a MSnSet structure


# process
# get filtered quantitation data
qnt <- filterNA(abun)
# process it
processingData(qnt)
# get protein quantitation data
#    we will group quantitation data by protein accession code
#	this should have no effect as we are already working at the
#	protein level
protqnt <- combineFeatures(qnt,
                           groupBy = fData(qnt)$accession,
                           method = sum)
# and indeed, the values remain untouched when working at the protein level

# plot N random proteins 
n.prots <- 8          

# allocate 1/4 of the page for text
#layout(matrix(c(1, 1, 2, 2, 2, 2, 2, 2), 4, 2, byrow = TRUE))
#layout.show(2)
#out.text(labels=paste("Protein intensity for", n.prots, "last proteins", sep=' '))

set.seed(16160423)		# Shakespeare, Cervantes, Inca Garcilaso die
n.random.prots <- sample( 1:dim(protqnt)[1], n.prots, replace=F )
exprs.to.plot <- exprs(protqnt)[n.random.prots, ]
acc.to.plot <- rownames(protqnt)[n.random.prots]

cls <- brewer.pal(n.prots, "Set1")
matplot(t(exprs.to.plot), type = "b",
        lty = 1, col = cls,
        ylab = "Protein intensity (summed peptides)",
        xlab = "Measure group")
legend("topright", tail(featureNames(protqnt), n=n.prots),
       lty = 1, bty = "n", cex = .8, col = cls)
       
       
# normalyze
qntS <- normalise(qnt, "sum")
qntV <- normalise(qntS, "vsn")
qntV2 <- normalise(qnt, "vsn")
# plot SD vs means (by row)
#meanSdPlot(qntV)
#meanSdPlot(qntV2)

# Pick a few datums to display


# fData(qnt)$accession gets the accession column from the data in qnt
idx <- sapply(acc.to.plot, grep, fData(qnt)$accession)
idxs <- sapply(idx, head, 3)
small <- qntS[unlist(idxs), ]	# unlist idxs and extract them from qntS

idxm <- sapply(idx, head, 10)
medium <- qntV[unlist(idxm), ]

s <- exprs(small)	# from qntS and head(3)
m <- exprs(medium)	# from qntV and head(10)
head(s)
head(m)

# use more convenient identifiers
colnames(s) <- c("WT1", "WT2",
                 "MT1", "MT2")
colnames(m) <- c("WT1", "WT2",
                 "MT1", "MT2")

# we must use rownames to identify unique proteins, so we use the accession
rownames(s) <- fData(small)$accession
rownames(m) <- fData(medium)$accession

# change chosen gene names to something more readable
#rownames(m)[grep("CYC", rownames(m))] <- "CYT"
#rownames(m)[grep("ENO", rownames(m))] <- "ENO"
#rownames(m)[grep("ALB", rownames(m))] <- "BSA"
#rownames(m)[grep("PYGM", rownames(m))] <- "PHO"
#rownames(m)[grep("ECA", rownames(m))] <- "Background"


# draw the heatmaps
cls <- c(brewer.pal(length(unique(rownames(m)))-1, "Set1"),
         "grey")
names(cls) <- unique(rownames(m))
wbcol <- colorRampPalette(c("white", "darkblue"))(256)

msg=paste('Heatmap of Sum-normalized', n.prots, 'random proteins')
heatmap(s, col = wbcol, RowSideColors=cls[rownames(s)], main=msg)
msg=paste('Heatmap of Sum+VSN-normalized', n.prots, 'random proteins')
heatmap(m, col = wbcol, RowSideColors=cls[rownames(m)], main=msg)


# draw spikes plots
dfr.s <- data.frame(exprs(small),
                  Protein = as.character(fData(small)$accession),
                  Feature = featureNames(small),
                  stringsAsFactors = FALSE)

colnames(dfr.s) <- c("WT1", "WT2", "MT1", "MT2",
                   "Protein", "Feature")
# set easier names
#dfr$Protein[dfr$Protein == "sp|P00924|ENO1_YEAST"] <- "ENO"
#dfr$Protein[dfr$Protein == "sp|P62894|CYC_BOVIN"]  <- "CYT"
#dfr$Protein[dfr$Protein == "sp|P02769|ALBU_BOVIN"] <- "BSA"
#dfr$Protein[dfr$Protein == "sp|P00489|PYGM_RABIT"] <- "PHO"
#dfr$Protein[grep("ECA", dfr$Protein)] <- "Background"
dfr.s.2 <- melt(dfr.s)
## Using Protein, Feature as id variables

ggplot(aes(x = variable, y = value, colour = Protein),
       data = dfr.s.2) +
  geom_point() +
  geom_line(aes(group=as.factor(Feature)), alpha = 0.5) +
  facet_grid(. ~ Protein) + theme(legend.position="none") +
  labs(x = "Reporters", y = "Normalised intensity") +
  ggtitle("Sum-normalized plots")

#ggsave(plot = myplot, filename = "myplot.pdf", device = "pdf")

dfr.m <- data.frame(exprs(medium),
                  Protein = as.character(fData(small)$accession),
                  Feature = featureNames(small),
                  stringsAsFactors = FALSE)

colnames(dfr.m) <- c("WT1", "WT2", "MT1", "MT2",
                   "Protein", "Feature")
# set easier names
#dfr$Protein[dfr$Protein == "sp|P00924|ENO1_YEAST"] <- "ENO"
#dfr$Protein[dfr$Protein == "sp|P62894|CYC_BOVIN"]  <- "CYT"
#dfr$Protein[dfr$Protein == "sp|P02769|ALBU_BOVIN"] <- "BSA"
#dfr$Protein[dfr$Protein == "sp|P00489|PYGM_RABIT"] <- "PHO"
#dfr$Protein[grep("ECA", dfr$Protein)] <- "Background"
dfr.m.2 <- melt(dfr.m)
## Using Protein, Feature as id variables

ggplot(aes(x = variable, y = value, colour = Protein),
       data = dfr.m.2) +
  geom_point() +
  geom_line(aes(group=as.factor(Feature)), alpha = 0.5) +
  facet_grid(. ~ Protein) + theme(legend.position="none") +
  labs(x = "Reporters", y = "Normalised intensity") +
  ggtitle("Sum+VSN-normalized plots")



#############################################################################
#			ANALYZE DIFFERENTIAL EXPRESSION
#############################################################################


# DEA
#
# https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/DEA.html

#m <- abun

#m <- filterNA(m)

# ...



# synapter
#
# library(synapter)
# synapterGuide()
#
# Trivial t-test based analysis
#
# we start from an MSnSet like protqnt, qnt, qntS, qntV or qntV2
# ...
#
# the TOP3 approach uses the 3 most intense peptides to compute protein
# intensity.
# As we have protein data already, we can skip it and initial filtering
# and normalisation.

if (exprColsPrefix == "No.PSMs") {
    syn_data <- qntV
} else if (exprColsPrefix == "NormAbun") {
    syn_data <- qntV2		# data is already sum-normalized
}

# convert expressions to log2
exprs(syn_data) <- log2(exprs(syn_data))

## apply a t-test and extract the p-value
pv <- apply(exprs(syn_data), 1 ,
            function(x)t.test(x[1:2], x[3:4])$p.value)

## calculate q-values
qv <- qvalue(pv)$qvalues

## calculate log2 fold-changes
lfc <- apply(exprs(syn_data), 1 ,
             function(x) mean(x[1:3], na.rm=TRUE)-mean(x[4:6], na.rm=TRUE))
## create a summary table
syn_res <- data.frame(cbind(exprs(syn_data), pv, qv, lfc))
## reorder based on q-values
syn_res <- syn_res[order(syn_res$q.value), ]
colnames(syn_res) <- c("wt1", "wt2", "mut1", "mut2", "p.value", "q.value", "log2FC")
knitr::kable(head(round(syn_res, 3)))

# save
out=paste("synapter", exprColsPrefix, "all.tab", sep='.')
write.table(syn_res, file=out, row.names=T, col.names=T)

# volcano plot
plot(syn_res$log2FC, -log10(syn_res$q.value),
     col = ifelse(grepl("mut", colnames(syn_res)),
       "#4582B3AA",
       "#A1A1A180"),
     pch = 19,
     xlab = expression(log[2]~fold-change),
     ylab = expression(-log[10]~q-value))
grid()
abline(v = -1, lty = "dotted")
abline(h = -log10(0.1), lty = "dotted")
legend("topright", c("mut", "wt"),
       col = c("#4582B3AA", "#A1A1A1AA"),
       pch = 19, bty = "n")
       
# heatmap
msg <- paste("synapter trivial analysis of", exprColsPrefix)
heatmap(as.matrix(syn_res[1:length(exprCols)]), main=msg)


#
#
#	DEP 
#
#

# We can use any of the protqnt/qntS/qntV/qntV2 datasets 
# calculated above or just process the data from scratch
######

# For use with MSnData instead of a table:
# see https://bioconductor.org/packages/3.18/bioc/html/DEP.html
# convert to SummarizedExperiment
#se <- as(abun, "SummarizedExperiment")
# to convert back, use as(se, "MSnSet")

# default is to start from raw data
from.abundances <- F
from.sum.norm <- F
from.sum.vsn.norm <- F
from.vsn.norm <- F

#if (exprColsPrefix == "NormAbun") {
#    from.abundances <- T
#    # data is already Sum normalized, we will do additional VSN below,
#    # so we do not need Sum nor VNS nor Sum+VSN
#}
# we'll prefer the default (start from scratch) so we can keep meta-information

if (from.abundances == T) {
    data_se <- as(protqnt, "SummarizedExperiment")
    exprs(data_se) <- log2(exprs(data_se))
} else if (from.sum.norm == T) {
    data_se <- as(qntS, "SummarizedExperiment")
    exprs(data_se) <- log2(exprs(data_se))
} else if (from.sum.vsn.norm == T) {
    data_se <- as(qntV, "SummarizedExperiment")
    exprs(data_se) <- log2(exprs(data_se))
} else if (from.vsn.norm == T) {
    data_se <- as(qntV2, "SummarizedExperiment")
    exprs(data_se) <- log2(exprs(data_se))
} else {		
    # work from scratch
    #
    # from data table
    #
    # https://bioconductor.org/packages/3.18/bioc/vignettes/DEP/inst/doc/DEP.html

    # read data
    data <- read.table('data.tab', header=T, sep='\t', quote='"')
    colnames(data)
    #
    # remove rows with NA
    # data <- na.omit(data)
    #

    # check there are no duplicate proteins
    data$accession %>% duplicated() %>% any()

    #data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
    data_unique <- make_unique(data, "accession", "Description", delim = ";")

    LFQ_columns <- grep(exprColsPrefix, colnames(data_unique)) # get LFQ column numbers

    # get experimental design (first four rows are Scaled Abundances)
    experimental_design <- read.table('metadata.tab', header=T)
    # fix labels to correct name
    experimental_design$label <- sub("Count", exprColsPrefix, experimental_design$label)
    experimental_design

    # copy "accession" to "ID" and "name" is we didn't use make_unique()
    #data_unique$name <- data_unique$accession
    #data_unique$ID <- data_unique$accession

    # convert to summarized experiment
    data_se <- make_se(data_unique, LFQ_columns, experimental_design)
    # assay data is log2 transformed and rownames correspond to protein names
    # name and ID should have been generated by make_unique
    # colData contains the experimental design in label, condition and replicate
    #	as well as a new ID column
    # IMPORTANT: it applies a log2 transform!!!
    # we do not need to add it later, but if we skip this, we need
    # log2 transformed expressions.
}



# filter missing values (proteins not identified in all replicas)
# first explore:
plot_frequency(data_se)
# we need to impute missing values, but only for proteins with not too many


# Filter for proteins that are identified in all replicates of at least 
# one condition
data_filt <- filter_missval(data_se, thr = 0)

# Less stringent filtering:
# Filter for proteins that are identified in 2 out of 3 replicates of at 
# least one condition
#data_filt2 <- filter_missval(data_se, thr = 1)

# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

# The data is background corrected and normalized by variance stabilizing
# transformation (vsn).
# for Sum-normalized data (NormAbun) this will result in Sum+VSN
# Normalize the data
data_norm <- normalize_vsn(data_filt)
meanSdPlot(data_norm)

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)

# The remaining missing values in the dataset need to be imputed. The data can
# be missing at random (MAR), for example if proteins are quantified in some
# replicates but not in others. Data can also be missing not at random (MNAR),
# for example if proteins are not quantified in specific conditions (e.g. in
# the control samples). MNAR can indicate that proteins are below the detection
# limit in specific samples, which could be very well the case in proteomics
# experiments. For these different conditions, different imputation methods
# have to be used, as described in the MSnbase vignette and more specifically
# in the impute function descriptions.
# 
# To explore the pattern of missing values in the data, a heatmap is plotted
# indicating whether values are missing (0) or not (1). Only proteins with at
# least one missing value are visualized.

# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

# To check whether missing values are biased to lower intense proteins, the
# densities and cumulative fractions are plotted for proteins with and without
# missing values.

# Plot intensity distributions and cumulative fraction of proteins with and 
# without missing values
plot_detect(data_filt)

# In this case, missing values do not seem to be biased to low-density when
# using SclAbun.
# With PSMs and NormAbun it does look MNAR.

# Indeed the proteins with missing values have on average low intensities. This
# data (MNAR and close to the detection limit) should be imputed by a
# left-censored imputation method, such as the quantile regression-based
# left-censored function (\u201cQRILC\u201d) or random draws from a
# left-shifted distribution (\u201cMinProb\u201d and \u201cman\u201d). In
# contrast, MAR data should be imputed with methods such as k-nearest neighbor
# (\u201cknn\u201d) or maximum likelihood (\u201cMLE\u201d) functions. See the
# MSnbase vignette and more specifically the impute function description for
# more information.

# All possible imputation methods are printed in an error, if an 
# invalid function name is given. The error will tell us which functions
# are available
#impute(data_norm, fun = "")

# Impute missing data using random draws from a Gaussian distribution 
# centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Impute missing data using random draws from a manually defined 
# left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)

# The effect of the imputation on the distributions can be visualized.
# Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)
plot_imputation(data_norm, data_imp_man)
plot_imputation(data_norm, data_imp_knn)

# we will use MinProb imputation (data_imp)
####
# if no imputation, then use filtered data
#data_imp <- data_filt
####


# Differential enrichment analysis
# 
# Protein-wise linear models combined with empirical Bayes statistics are used
# for the differential enrichment analysis (or differential expression
# analysis). The test_diff function introduced here uses limma and
# automatically generates the contrasts to be tested. For the contrasts
# generation, the control sample has to be specified. Additionally, the types
# of contrasts to be produced need to be indicated, allowing the generation of
# all possible comparisons (\u201call\u201d) or the generation of contrasts of
# every sample versus control (\u201ccontrol\u201d). Alternatively, the user
# can manually specify the contrasts to be tested (type = \u201cmanual\u201d),
# which need to be specified in the argument test.

# Differential enrichment analysis  based on linear models and empherical 
# Bayes statistics

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "wt")

# Test all possible comparisons of samples
data_diff_all_contrasts <- test_diff(data_imp, type = "all")

# Test manually defined comparisons
data_diff_manual <- test_diff(data_imp, type = "manual", 
                              test = c("wt_vs_mut", "mut_vs_wt"))

# in our simple setup we can do with just mut_vs_wt (data_diff)

# Finally, significant proteins are defined by user-defined cutoffs using 
# add_rejections.

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = p.cutoff, lfc = log2(log2cutoff))


# Visualization of the results
# 
# The results from the previous analysis can be easily visualized by a number
# of functions. These visualizations assist in the determination of the optimal
# cutoffs to be used, highlight the most interesting samples and contrasts, and
# pinpoint differentially enriched/expressed proteins.

# PCA plot
# 
# The PCA plot can be used to get a high-level overview of the data. This can
# be very useful to observe batch effects, such as clear differences between
# replicates.
# Plot the first and second principal components
DEP::plot_pca(dep, x = 1, y = 2, n = dim(dep)[1], point_size = 4)
# there is another plot_pca in MSnSet.utils which won't work with an SE


# Correlation matrix
# 
# A correlation matrix can be plotted as a heatmap, to visualize the Pearson
# correlations between the different samples.
# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")


# Heatmap of all significant proteins
# 
# The heatmap representation gives an overview of all significant proteins
# (rows) in all samples (columns). This allows to see general trends, for
# example if one sample or replicate is really different compared to the
# others. Additionally, the clustering of samples (columns) can indicate closer
# related samples and clustering of proteins (rows) indicates similarly
# behaving proteins. The proteins can be clustered by k-means clustering
# (kmeans argument) and the number of clusters can be defined by argument k.
# Plot a heatmap of all significant proteins with the data centered per 
# protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))

# Alternatively, a heatmap can be plotted using the contrasts, i.e. the direct
# sample comparisons, as columns. Here, this emphasises the enrichment of
# mutant compared to the control sample.

# Plot a heatmap of all significant proteins (rows) and the tested 
#contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)


# Volcano plots of specific contrasts
# 
# Volcano plots can be used to visualize a specific contrast (comparison
# between two samples). This allows to inspect the enrichment of proteins
# between the two samples (x axis) and their corresponding adjusted p value (y
# axis). The add_names argument can be set to FALSE if the protein labels
# should be omitted, for example if there are too many names.
# Plot a volcano plot for the contrast "mut vs wt""
DEP::plot_volcano(dep, contrast = "mut_vs_wt", label_size = 4, 
	add_names = TRUE)
## Warning: ggrepel: 58 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps

# Barplots of a protein of interest
# 
# It can also be useful to plot the data of a single protein, for example if
# this protein is of special interest.
# # PSMs
# extreme <- c('Q9S746', 'Q93VG5', 'Q9LJN1',
# 	     'O04311', 'A0A178W0D3')
# # SclAbun
# extreme <- c('Q9LUD4', 'Q93VG5', 'F4KC80', 'P49107', 'O80653', 'Q8GYW0', 'A0A1P8B2M2',
# 	     'Q9SR37', 'O04310', 'O04311')
# # NormAbun
# extreme <- c('F4KC80', 'Q9LUD4', 'P49107', 'Q93VG5', 'O80653', 'Q8GYW0', 'A0A1P8B2M2',
# 	     'Q9SR37', 'O04310', 'O4311')
# general: get significant results IDs
extreme <- (get_results(dep) %>% filter(significant))$name

# Plot a barplot for the extreme proteins
plot_single(dep, proteins = extreme)
# Plot a barplot for the proteins with the data centered
plot_single(dep, proteins = extreme, type = "centered")

# Frequency plot of significant proteins and overlap of conditions
# 
# Proteins can be differentially enriched/expressed in multiple comparisons. To
# visualize the distribution of significant conditions per protein and the
# overlap between conditions, the plot_cond function can be used.
#
# Plot a frequency plot of significant proteins for the different conditions
#plot_cond(dep)
# this is of little use as here we only have one condition

# Results table
# 
# To extract a table containing the essential results, the get_results function
# can be used.
# Generate a results table
data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()
colnames(data_results)

# Of these columns, the p.val and p.adj columns contain the raw and adjusted p
# values, respectively, for the contrast as depicted in the column name. The
# ratio columns contain the average log2 fold changes. The significant columns
# indicate whether the protein is differentially enriched/expressed, as defined
# by the chosen cutoffs. The centered columns contain the average log2 fold
# changes scaled by protein-wise centering.

# Generate a data.frame from the resulting SummarizedExperiment object
# 
# You might want to obtain an ordinary data.frame of the results. For this purpose, the package provides functions to convert SummarizedExperiment objects to data.frames. get_df_wide will generate a wide table, whereas get_df_long will generate a long table.

# Generate a wide data.frame
df_wide <- get_df_wide(dep)
# Generate a long data.frame
df_long <- get_df_long(dep)

# Save our results object for reuse
# 
# To facilitate future analysis and/or visualization of our current data,
# saving our analyzed data is highly recommended. We save the final data
# object (dep) as well as intermediates of the analysis, i.e. the initial
# SummarizedExperiment object (data_se), normalized data (data_norm), imputed
# data (data_imp) and differentially expression analyzed data (data_diff). This
# allows us to easily change parameters in future analysis.

# Save analyzed data
out <- paste("DEP", exprColsPrefix, "RData", sep='.')
save(data_se, data_norm, data_imp, data_diff, dep, file=out)
# These data can be loaded in future R sessions using this command
load(out)

# save table
out <- paste("DEP", exprColsPrefix, "diff.tab", sep='.')
write.table( df_wide[df_wide$significant==T, ], file=out, row.names=T, col.names=T)

# simplify
out <- paste("DEP", exprColsPrefix, "diff.simple.tab", sep='.')
write.table( df_wide[df_wide$significant==T, c("name", "Description", "mut_vs_wt_diff", "mut_vs_wt_p.adj")], 
	file=out, row.names=F, col.names=T)


#
#
#	msmsTests using edgeR
#
#

# https://support.bioconductor.org/p/104999/
# we can also use package edgeR through "msmsTests"
# https://bioconductor.org/packages/release/bioc/manuals/msmsTests/man/msmsTests.pdfBiocManager::install("msmsTests")

# we start from the raw data as edgeR will do all the processing for us
res <- abun

res <- filterNA(res)

null.f <- "y~replica"
alt.f <- "y~status+replica"
div <- apply(exprs(res), 2, sum)

poisson=F
if (poisson == T) {
    # Differential expression tests on spectral counts
    # 
    # Spectral counts (SpC) is an integer measure which requires of tests suited to
    # compare counts, as the GML methods based in the Poisson distribution, the
    # negative-binomial, or the quasilikelihood [11]
    # 
    # Generally speaking no test is superior to the other. They are just more or
    # less indicated in some cases. The Poisson regression requires the estimation
    # of just one parameter, and is indicated when the number of replicates is low,
    # two or three. A drawback of the Poisson distribution is that its variance
    # equals its mean, and it is not able to explain extra sources of variability a
    # part of the sampling. Quasi-likelihood is a distribution independent GLM, but
    # requires of the estimation of two parameters and hence needs a higher number
    # of replicates, i.e. not less than four. The negative-binomial requires two
    # parameters too, but we may use the implementation of this GLM in the edgeR
    # package [12] which uses an empirical Bayes method to share information across
    # features and may be employed with a restricted number of replicates; in the
    # worst case it limits with the Poisson solution.

    # Poisson GLM regression 
    # 
    # When using the Poisson distribution we implicitly accept a model not
    # sensitive to biolog- ical variability [11]. So it is just recommended in
    # cases where we have very few replicates, if any, and we do not expect a
    # signicant biological variability between samples.

    ### Remove all zero rows
    e <- pp.msms.data(res)
    dim(e)

    ### Null and alternative model
    #null.f <- "y~replica"
    #alt.f <- "y~status+replica"
    null.f <- "y~1"
    alt.f <- "y~condition"
    ### Normalizing condition
    div <- apply(exprs(e),2,sum)
    ### Poisson GLM
    pois.res <- msms.glm.pois(e,alt.f,null.f,div=div)
    str(pois.res)

    ### DEPs on unadjusted p-values
    sum(pois.res$p.value<=0.01)

    ### DEPs on multitest adjusted p-values
    adjp <- p.adjust(pois.res$p.value,method="BH")
    sum(adjp<=0.01)

    ### The top features
    o <- order(pois.res$p.value)
    head(pois.res[o,],20)
}

qlglm=F
if (qlglm == T) {
    # Quasi-likelihood GLM regression
    # 
    # The quasi-likelihood is a distribution free model that allows for overdispersion, and could
    # be indicated where an appreciable source of biological variability is expected. In this
    # model, instead of specifying a probability distribution for the data, we just provide a
    # relationship between mean and variance. This relationship takes the form of a function,
    # with a multiplicative factor known as the overdispersion, which has to be estimated from
    # the data [11]. Its use in proteomics has been documented by Li et al. (2010)[13].

    ### Quasi-likelihood GLM
    ql.res <- msms.glm.qlll(res,alt.f,null.f,div=div)
    str(ql.res)

    ### DEPs on unadjusted p-values
    sum(ql.res$p.value<=0.01)

    ### DEPs on multitest adjusted p-values
    adjp <- p.adjust(ql.res$p.value,method="BH")
    sum(adjp<=0.01)

    ### The top features
    o <- order(ql.res$p.value)
    head(ql.res[o,],20)

}



# edgeR: negative binomial GLM regression example
# 
# The negative-binomial provides another model that allows for overdispersion.
# The im- plementation adopted in this package is entirely based in the
# solution provided by the package edgeR [12] which includes empirical Bayes
# methods to share information among features, and thus may be employed even
# when the number of replicates is as low as two. The negative-binomial is
# downward limited, when no overdispersion is observed, by the Poisson
# distribution.

eR.res <- msms.edgeR(res,alt.f,null.f,div=div,fnm="status")
str(eR.res)

### DEPs on unadjusted p-values
sum(eR.res$p.value <= 0.01)	# see how many have p <= 0.01

### DEPs on multitest adjusted p-values
adjp <- p.adjust(eR.res$p.value,method="BH")
sum(adjp <= 0.01)
eR.res$p.adjust <- adjp

### The top N features
o <- order(eR.res$p.value)
head(eR.res[o,], n.prots)

# Let's now convert the data to counts if needed by multyplying by a factor:
#exprs(res) <- round(exprs(res) * 10)	# we got ScaledAbun with one decimal
#x <- msms.edgeR(res, alt.f, null.f, div = div, fnm = "status")


# see fature data for significant differences
head( rownames(eR.res[ eR.res$p.value <= p.cutoff, ]) )	# significant diffs
head( rownames(fData(res)) %in%  rownames(eR.res[ eR.res$p.value < p.cutoff, ]) )
head( fData(res)[rownames(fData(res)) %in%  rownames(eR.res[ eR.res$p.value <= p.cutoff, ]), ] )

# save p <= p.cutoff
sigp <- fData(res)[rownames(fData(res)) %in%  rownames(eR.res[ eR.res$p.value <= p.cutoff, ]), ] 
# add statistics
sigp <- merge(sigp, eR.res, by='row.names', all.x=T, all.y=F)[ , -1]
rownames(sigp) <- sigp$accession
out <- paste("eR.signif", exprColsPrefix, "p", p.cutoff, "tab", sep='.')
write.table( sigp, file=out, row.names=F, col.names=T)

# save padj <= 0.01
sigpa <- fData(res)[rownames(fData(res)) %in%  rownames(eR.res[ eR.res$p.adjust <= 0.01, ]), ] 
sigpa <- merge(sigpa, eR.res, by='row.names', all.x=T, all.y=F)[ , -1]
rownames(sigpa) <- sigpa$accession
out <- paste("eR.signif", exprColsPrefix, "padj=0.01.tab", sep='.')
write.table( sigpa, file=out, row.names=F, col.names=T)

# simplify
out <- paste("eR.signif", exprColsPrefix, "p", p.cutoff, "simple.tab", sep='.')
write.table(sigp[ , c("accession", "Description", "p.value", "p.adjust", "LogFC")],
	   file=out, row.names=F, col.names=T)

out <- paste("eR.signif", exprColsPrefix, "padj=0.01.simple.tab", sep='.')
write.table( sigpa[ , c("accession", "Description", "p.value", "p.adjust", "LogFC")], 
	   file=out, row.names=F, col.names=T)

# Reproducibility
# 
# In the omics field, reproducibility is of biggest concern. Very low p-values
# for a protein in an experiment are not enough to declare that protein as of
# interest as a biomarker. A good biomarker should give as well a reproducible
# signal, and posses a biologically signicant eect size. According to our
# experience, a protein giving less than three counts in the most abundant
# condition results of poor reproducibility, to be declared as statistically
# signicant, experiment after experiment. On the other hand most of the false
# positives in an spiking experiment show log fold changes below 1. These two
# observations [6] allow to improve the results obtained in the previous
# sections. The trick is to flag as relevant those proteins which have low
# p-values, high enough signal, and good effect size. In performing this
# relevance filter we may even accept higher adjusted p-values than usual. This
# flagging is provided by the function test.results.

### Cut-off values for a relevant protein as biomarker
alpha.cut <- 0.05		# now we do not need to use 0.01
SpC.cut <- 2
lFC.cut <- 1
### Relevant proteins according to previous adjustments
eR.tbl <- test.results(eR.res,			# the dataframe with the stats
			res, 			# the MSnSet with the data
                        pData(res)$status,	# the factor used in the tests
                        "mut",			# treatment level name
                        "wt",			# contol level name
                        div,			# weights used as divisors
			alpha=alpha.cut,minSpC=SpC.cut,minLFC=lFC.cut,
			method="BH")$tres

(eR.nms <- rownames(eR.tbl)[eR.tbl$DEP])

# without the post-test filter a relatively low adjusted p-value cut-off is
# required to keep an acceptable number of false positives. The post-test filter
# allows to relax the p-value cut-off improving at the same time both the number
# of true positives and false positives. 

# save the results
eR.wide <- merge(fData(res), eR.tbl, by='row.names', all=F)
rownames(eR.wide) <- eR.wide$accession
eR.wide <- eR.wide[ , -1]

out <- paste("eR.mark", exprColsPrefix, "p", alpha.cut, "tab", sep='.')
write.table(eR.wide, file=out, row.names=F, col.names=T)

out <- paste("eR.mark", exprColsPrefix, "p", alpha.cut, "simple.tab", sep='.')
write.table(eR.wide[ , c("accession", "Description", "p.value", "p.adjust", "LogFC")],
	   file=out, row.names=F, col.names=T)

# A useful tool to visualize the global results of dierential expression tests
# is a table of accumulated frequencies of features by p-values in bins of log
# fold changes. It may help in nding the most appropriate post-test lter cut-o
# values in a given experiment.

pval.by.fc(eR.tbl$adjp,eR.tbl$LogFC)

### Filtering by minimal signal
flt <- eR.tbl$wt > 2
pval.by.fc(eR.tbl$adjp[flt], eR.tbl$LogFC[flt])


# Another usual tool is a volcanoplot with the ability to visualize the eect of
# dierent post-test lter cut-o values.
par(mar=c(5,4,0.5,2)+0.1)
res.volcanoplot(eR.tbl,max.pval=0.05,min.LFC=1,maxx=3,maxy=NULL,
		ylbls=3)



dev.off()

