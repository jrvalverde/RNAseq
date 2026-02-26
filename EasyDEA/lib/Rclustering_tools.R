
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
    cat("
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

    cat("
    
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

