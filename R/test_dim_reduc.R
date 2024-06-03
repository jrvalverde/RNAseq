# NOTE: requires starting R with
#	R --max-ppsize=500000
# to be able to make linear models.
#
# libraries to manipulate formulas
library(tools)
library(rlang)
# libraries for variable reduction
library(caret)
library(randomForest)
library(xgboost)
library(relaimpo)
library(earth)
library(Boruta)
library(rFerns)
library(randomForest)
library(DALEX)
#library(DALEXtra)
library(vita)
library(lasso)
# libraries for factor analysis
library(corrplot)
library(factoextra)
library(ctv)
library(psych)
library(psychTools)
library(phyloseq)
library(nFactors)
library("FactoMineR")
library("factoextra")
ibrary("gplots")
library("corrplot")
library("ggpubr")
library("ade4")
library("ExPosition")
library(vegan)
library(dplyr)

# the next ones are no longer maintained and do not work, do not use
# library(devtools)
#    install_github("tomasgreif/woe")
#    install_github("tomasgreif/riv")
# library(woe)
# library(riv)

#
# CONFIGURATION
#

SEED <- 20240429

# whether we want to test the functions in this file
test_coturnix <- FALSE
test_gallus2 <- TRUE

# IMPORTANT NOTE: 
#	The next two variables should be defined in the conditional
# block corresponding to the experiment we want to analyze:
# if (test_coturnix) or if (test_gallus2) which are below, after
# the section that defines the functions that we use. We should
# scroll down to those sections and define the correct values
# there.
#
#ORG <- 'Cjaponica'
#formula <- pgu ~ .
#
#ORG <- 'Ggallus2'
#formula <- sample ~ .


##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################

# These are some convenience functions to help with printing messages, loading libraries and plotting

#' my.name()
#'
#' my.name() obtains the name of the function that called it
#'
#' This function will return the name of the calling function.
#' The CALLING function !!!
#' THE CALLING FUNCTION !!!!!!
#' It is intended as a way to obtain a function's own name for error
#' reporting. 
#  We could do it directly, within a function using 
#  my.name <- deparse(sys.call())
#  but it would be difficult to understand. Isolating this in a function
#  call makes it more readable.
#'
#' 
#' @return	the name of the calling function (one up in the call stack)
#'
#' @usage	me <- my.name()
#'
#' @examples
#'	f <- function() { print(my.name()) }
#'	f()
#'	##[1] "print(my.name())"
#'	## my.name() is being passed to print() as an argument and is thus
#'	## called by 'print', hence the output
#'	##
#'	f <- function() { me <- my.name() ; print(me) }
#'	f()
#'	##[1] "f()"
#'	## my.name() is called by f() to obtain its own name and then the name
#'	## is handled to print().
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
my.name <- function(verbose=F) { 
	my.call.number <- sys.nframe()
	if (my.call.number == 1) {
	    # we are at the script/command line top level
		top.file <- sys.frame(1)$ofile
	    if (is.null(top.file))
		    caller <- "R"
		else {
		    caller <- top.file
            if (verbose == FALSE) caller <- basename(caller)
		}
	} else {
    	caller <- deparse(sys.call(-1)) 
		if (! verbose) caller <- gsub('\\(.*', '', as.character(caller))
		#who.called.me <- deparse( sys.calls()[[ my.parent.number ]] )
		#return( who.called.me )
    }
    return(caller)
}

# alternative (and simpler) implementation
me <- function(..., verbose=F) {
    if (sys.nframe() == 1) callingFun <- "R"
    else callingFun = as.list(sys.call(-1))[[1]]
	
    calledFun = as.list(sys.call())[[1]]
    if (verbose) message(paste(callingFun, " is calling ", calledFun, sep=""))
	return(callingFun)
}


my.caller <- function(verbose=F) {
	my.call.number <- sys.nframe()
	if (my.call.number <= 2) {
	    # we are at the script/command line top level
		# or colled by a function whose parent is the
		# top level
		top.file <- sys.frame(1)$ofile
	    if (is.null(top.file))
		    caller <- "R"
		else {
		    caller <- top.file
            if (verbose == FALSE) caller <- basename(caller)
		}
	} else {
    	caller <- deparse(sys.call(-2)) 
		if (! verbose) caller <- gsub('\\(.*', '', as.character(caller))
		#who.called.me <- deparse( sys.calls()[[ my.parent.number ]] )
		#return( who.called.me )
    }
    return(caller)
}


#' cat.err()
#'
#' err() prints a generic error message using cat()
#'
#' @param	abort	(boolean) whether the program should be stopped
#'			(defaults to FALSE)
#' @param	...	The message to print and cat options (see cat())
#'
#' @return	whatever cat() returns
#'
#' @usage	cat.err('some message\n', sep='')
#'
#' @examples
#'		f <- function(x) { if (x < 10) {print(x)} else {cat.err(x, "is too big\n")}}
#'		f(13)
#'		##ERROR in  f(13) : 13 is too big
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
cat.err <- function(..., abort=FALSE, verbose=T) { 
    # error messages should be always printed
	# so we will ignore the verbose argument
	verbose <- TRUE

    # caller <- mycaller()
    #caller <- deparse( sys.calls()[[sys.nframe()-1]] )
    caller <- my.caller(verbose=verbose)
    cat('ERROR in ', caller, ":", ...)
    if (abort) {
        quit(save="no", status=1, runLast=FALSE)
    }
}


#' cat.warn()
#'
#' cat.warn() prints a generic warning message using cat()
#'
#' @param	...	The message to print and cat options (see cat())
#'
#' @return	whatever cat() returns
#'
#' @usage	cat.warn('some message\n', sep='')
#'
#' @examples
#'		f <- function(x) { if (x < 10) {print(x)} else {cat.warn(x, "is too big\n")}}
#'		f(13)
#'		##WARNING in  f(13) : 13 is too big
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
cat.warn <- function(..., verbose=F) { 
    #caller <- deparse( sys.calls()[[sys.nframe()-1]] )
    #me <- as.list(sys.call())[[1]]
    #parent <- as.list(sys.call(-1))[[1]]
	caller <- my.caller(verbose=verbose)
	cat('WARNING in ', caller, ":", ...)
}


#' cat.info()
#'
#' cat.info() prints a generic information message using cat()
#'
#' @param	...	The message to print and cat options (see cat())
#'
#' @return	whatever cat() returns
#'
#' @usage	cat.warn('some message\n', sep='')
#'
#' @examples
#'		f <- function(x) { if (x < 10) {print(x)} else {cat.info(x, "is too big\n")}}
#'		f(13)
#'		##INFO  f(13) : 13 is too big
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#

cat.info <- function(..., verbose=F) { 
    #caller <- deparse( sys.calls()[[sys.nframe()-1]] )
    #me <- as.list(sys.call())[[1]]
    #parent <- as.list(sys.call(-1))[[1]]
    caller <- my.caller(verbose=verbose)
    cat('INFO ', caller, ":", ...)
}


# to be used instead of library(): this function ensures
# that the package is installed if not present in the system
#' use.package
#'
#' checks if a package is installed and loads it, if it is not then it 
#' will install it and then load it
#' 
#' @param	p	a single package name to load
#'
#' @param	silent	whether we want messages produced during the load displayed
#' 
#' @return	the package is guaranteed to be installed and loaded
#'
#' @usage	use.package(p, silent=T)
#' 
#' @examples	use.package('stringr')
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
use.package <- function(p, silent = F, verbose = T) {
    ### NOTE: we should likely try CRAN first, Bioconductor afterwards,
    ### install.packages offline, devtools::install and install_github
    ### for now we will kepp it simple
    if (!is.element(p, installed.packages()[,1])) {
        # we prefer BiocManager because it can install both Bioconductor
        # and CRAN packages
        #
        #install.packages(p, dep = TRUE)
        if (!is.element('BiocManager', installed.packages()[,1])) {
            install.packages('BiocManager', dependencies=TRUE)
        }
        BiocManager::install(p, dep = TRUE)

    }
    # this will fail for 'banneR' which is neither on CRAN nor on
    # Bioconductor, so we'll check for it as a special case for now
    if (p == 'banneR')
       library(banneR)		# if not available we will do without it
    if (verbose)
        require(p, character.only = TRUE)
	else
        suppressPackageStartupMessages(require(p, character.only = TRUE))
}



#' use.packages
#'
#' given a vector with one or more package names, it first verifies if all
#' are installed, if not, then missing packages are installed, and finally
#' all the packages are loaded
#' 
#' @param	p	vector of package names to load
#' 
#' @return	the packages are guaranteed to be installed and loaded
#'
#' @usage	use.packages(pkgs)
#' 
#' @examples	use.packages(c('stringr', 'tidyr'))
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
use.packages <- function(pkgs){
    new.pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
    if (length(new.pkgs))
        # if there are uninstalled packages
        install.packages(new.pkgs, dependencies = TRUE)

    sapply(pkgs, require, character.only = TRUE)
}


#' as.png
#'
#' run specified code to generate a graphic and save as PNG file
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) José R. Valverde, CNB-CSIC, 2019
#'
#' @license	EU-GPL
#'
#' @export
#
as.png <- function(PLOT=NULL, 
               file='out.png', width=1024, height=1024, 
               overwrite=TRUE, verbose=T) {

    if (is.null(file)) {
        tryCatch(print(PLOT))
    } else if ( overwrite || ! file.exists(file) ) {
        if (verbose){
            cat("as.png(): creating", file, "\n")
        }
	    tryCatch( {
                png(file, width=width, height=height)
                print(PLOT)
            },
            finally = dev.off()
        )
    }
    return()
}


#' get_annotation
#'
get_annotation <- function(ensembl.db=NULL, biomart.db=NULL, verbose=F) {
    if (verbose) cat.info('loading annotation\n')
	
    ensembl.annot <- NULL
    bm.annot.1 <- NULL
    if (! is.null(biomart.db) && file.exists(biomart.db)) {
	    if (verbose) cat.info('loading', biomart.db, '\n')
        # get file extension and convert to lower case to simplify checks
		ext <- tolower(file_ext(biomart.db))
        if (ext == 'tab' || ext == 'tsv' || ext == 'txt') {
            # from .../signif: '../../biomaRt.annotation.1st.txt'
            bm.annot.1 <- read.table(file=biomart.db, header=T)
        } else if (ext == 'rds') {
            bm.annot.1 <- readRDS(biomart.db)
        } else {
            cat.warn('extension of biomart database must be txt, tsv, tab or rds\n')
        }
    }
    if (! is.null(ensembl.db) && file.exists(ensembl.db)) {
	    if (verbose) cat.info("Loading ENSMBL annotation from", ensembl.db, '\n')
        # get file extension and convert to lower case to simplify checks
        ext <- tolower(file_ext(ensembl.db))
        if (ext == 'tab' || ext == 'tsv' || ext == 'txt') {
            # from .../signif: '../../../../net.EnsDb/net.ens.db.rds'
            ensembl.annot <- readRDS(file=ensembl.db)
        } else if (ext == 'rds') {
            ensembl.annot <- readRDS(ensembl.db)
        } else {
            cat.warn('extension of ensembl database must be txt, tsv, tab or rds\n')
        }
    }
    # we could also load an endsb.org.db package
    # and then (but this currently doesn't work for some reason)
    #ens.annot <- ensembldb::select(ens.db, column='GENEID', keytype='GENEID', keys=rownames(lfc), columns='SEQNAME')
    # ensembldb::supportedFilters(ens.db)
    # columns(edb)
    # listColumns(edb)
    # keytypes(edb)
    # gids <- keys(edb, keytype = "GENEID")
    # length(gids)
    #ens.ann <- AnnotationDbi::select(ens.db, 
    #              column='GENID', keytype= 'GENEID', keys=rownames(lfc), 
    #              columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION', # no longer available
    #                         'GENENAME', 'GENEID', 'ENTREZID', # empty
    #                         'TXID', 'TXBIOTYPE', # these make the call fail
    #                         'PROTEINID', 'UNIPROTID' # no longer available
    #                          )
    #              )

    return( list(ensembl=ensembl.annot, biomart=bm.annot.1) )
}

#' enrich_genes
#'
enrich_genes <- function(data, ann, data.genes='ensembl_gene_id', ann.genes=data.genes) {
    # add the annotation as additional columns to data, using as the
    # column data.genes from data as index into ann.genes to find the
    # corresponding annotation.
    #
    # use as
    #	annotated.data <- enrich_genes(data, "ensembl.gene.id", ensembl.ann, "GENEID")
    #	annotated.data <- enrich_genes(data, "ensembl.gene.id", biomart.ann, "ensembl_gene_id")
    #   annotated.data <- enrich
	# for ENTREZGENE data
	# annotated.data <- enrich_genes(data, "gene.id", biomart.ann, "ensembl_gene_id")
	
    enriched_genes <- cbind(data, ann[match(data[ , data.genes], ann[ ,ann.genes]), ]) 

    return(enriched_genes)
}


load_coturnix_l2fc_data <- function(file="all_log2FC_transpose.txt") {
	# return a list with the datasets that we will use in the subsequent analysis


    # Load log2FC of significantly altered genes:
    #df <- read.table("all_log2FC.txt",header=T)
    # unused    df <- read.table("all_log2FC.txt",header=T)
    fd <- read.table(file, header=T)

    # this ensures that rows get the proper PFU irrespective
    # of order
    #	1. remove prefix
    pfu <- sub("wt.*pfu", "", fd$gene)
    #	2. remove suffix
    pfu <- sub("_l2fc", "", pfu)


    # set name of first column
    #names(fd)[1] <- 'wt.vs.PFU'
    #fd[,1] <- c(0.1, 1.0, 10.0, 100.0)

    # substitute column1 (file names) by column "$pfu"
    fd <- cbind(pfu=as.numeric(pfu), fd[ , -1])

    # define working data sets
    # 1. add a row of all zeros for the wild.type
    #    (inf.pfu/wt.vs.pfu = 0, l2fc=0 for all genes)
    lfc  <- as.data.frame(rbind(rep(0,dim(fd)[2]), fd))
    # 2. remove columns with NAs
    lfc.noNA <- lfc[ , colSums(is.na(lfc))==0]			# remove NAs
    clfc <- lfc.noNA
    # we now have
    #	fd 	-> the log2FC transposed data for infected cells (4x8594)
    #	lfc -> the log2FC transposed data for wt and infected cells (5x8594)
    #   clfc -> the log2FC transposed data for wt and infected cells without NA-containing columns (5x1182)
    return( list( signif.l2fc.all.data=lfc, signif.l2fc.common.data=clfc ) )
}


load_coturnix_norm_count_data <- function(file="normalized_counts.tab", select=NULL) {
    # load normalized counts of all genes
    nc <- read.table(file=file, row.names=1)
    tnc <- t(nc)
    # prepend a column with infection levels
    # upon read columns are not preserved
    #tnc <- cbind(c(rep(0, 3), rep(0.1, 3), rep(1.0,3), rep(10.0,3), rep(100.0,3)), tnc)
    # so we will fix the labels
    pfu <- rownames(tnc)
    pfu <- sub('^X', '', pfu)
    pfu <- sub('pfu.*', '', pfu)
    pfu <- sub('wt.*', '0.0', pfu)
    pfu <- as.numeric(pfu)
    
    # amend nc and tnc to include a 'pfu' row/column
    nc <- rbind(pfu=pfu, nc)
    tnc <- cbind(pfu=pfu, tnc)
    tnc <- as.data.frame( tnc )		# just in case

    if ( ! is.null(select)) {
        # select only counts of sigfigicantly altered genes:
        sig.nc <- rbind(pfu=nc[1, ], nc[rownames(nc) %in% select])
        sig.tc <- cbind(pfu=tnc[ ,1], tnc[ , colnames(tnc) %in% select ] )
    #tc <- sig.tc
    } else {
        sig.nc <- nc
        sig.tc <- tnc
    }

    # tc now contains the significant transposed normalized counts	(15x8954)
    # nc contains the significant normalized counts	(8954x15)

    return(list( all_norm_counts=nc, transp_norm_counts=tnc,
                 signif_norm_counts=sig.nc,
                 signif_transp_norm_counts=sig.tc))
}

get_gallus_significant_genes <- function(pattern='signif.*.tab', verbose=T) {
    # one way is to load all saved significant comparisons and find
	# the union of the gene names
	files <- list.files(pattern='signif_sorted_sample.*_annotated.tab')
	genes <- c()
	data <- list()
	for (i in files) {
	    df <- read.table(i, header=T)
		# gene annotation in gallus was defined for all genes only in 'gene.name'
		genes <- union(genes, df$gene.name)
		data[[i]] <- df
		# safety check
		if (verbose) cat(i,sum(!is.na(df$gene)), 'OK', sum(is.na(df$gene)), 'KO', length(genes), 'tot\n')
	}
	data$sig.genes <- genes
	return(data)
}

# for gallus we will only use normalized count data
# but this function will load count data for all genes, irrespective
# of whether they are significant or not, that is why we need a select
# argument
load_gallus_norm_count_data <- function(file="normalized_counts.tab", select=NULL) { 
    # for gallus we have also annotated normalized counts that we can use
	nc <- read.table(file, header=T, row.names=1, sep='\t')
	
	# select significant genes
    if ( ! is.null(select)) {
        # select only counts of sigfigicantly altered genes:
        sig.nc <- nc[rownames(nc) %in% select, ]
    #tc <- sig.tc
    } else {
        sig.nc <- nc
    }
	return(list(total.counts=nc, significant.counts=sig.nc))
}

#
# -- Correlation Method
#
# Pearson Correlation
imp_cor <- function (formula, data, verbose=T) {
    # get the left hand side of the formula (the dependent variable) as a sub-formula
    fy <- f_lhs(formula)
    # convert to column names
    ny <- all.vars(fy)
    # get the right hand side (the predictor variables) as a sub-formula
    fx <- f_rhs(formula)
    # convert to column names
    nx <- all.vars(fx)
    if ( nx == '.' ) {
       nx <- colnames(data)[ colnames(data) != ny ]
    } 
    # get x and y
    y <- data[ , ny]
    x <- data[ , nx]

    if(verbose) cat.info("calculating correlation between variables\n")

    data_cor <- cor(x, y)	

    if (verbose) cat.info("selecting variables with cor > 0.8\n")
    highlyCorrelated <- rownames(data_cor[data_cor[, 1] > 0.8, , drop = FALSE])

    return(highlyCorrelated)
}



# -- relimp
#library(relaimpo)
#
# Apply a linear regression model and get the relative importance of 
# each variable using R² partitioned by averaging over orders. Observations
# with missing values are automatically excluded.
#
imp_relimp <- function(formula, data, tune=F, verbose=T) {
    if (tune == F) {
        if (verbose) cat.info("building linear model\n")
        # directly build a linear model
        regressor <- lm(formula, data=data) # fit a lm() model

    } else { 
        if (verbose) cat.info("building linear model with CARET\n")
         # tune with caret
        ### XXX WARNING!!!
        ### MAGIC NUMBERS!!!
        ### WE need to think if this is general or not
        control <- trainControl(method="repeatedcv", number=10, repeats=3)
        tuneGrid <- expand.grid(k=1:3, size=seq(5,20,by=1))	
        regressor <- train(formula, data=data, 
                           method="lm",  # Use "lm" for linear regression
                           preProcess="scale", 
                           trControl=control#,
                           #tuneGrid=tuneGrid
                           )
        # if verbose, plot the regressor
        # use try() to ignore errors at this step.
        if (verbose) try(plot(regressor))		
        if (verbose) cat.info("calculating importance\n")
        if (FALSE) {
            # this is variable importance using CARET
            #result <- calc.relimp(regressor$finalModel, type = "lmg", rela = TRUE) # check how to retrieve the lm from caret result
            result <- caret::varImp(regressor, conditional=TRUE) 	# conditional=True, adjusts for correlations 
    													 	        # between predictors
            return(result)
       } else {
            #to use relimp, we need to access directly the  model
            regressor <- regressor$finalModel
        }
    }
    if (verbose) try(summary(regressor))
    if (verbose) cat.info("calculating importance\n")
    relImportance <- calc.relimp(regressor, type = "lmg", rela = TRUE) # calculate relative importance scaled to 100

    return(sort(relImportance$lmg, decreasing=TRUE)) # relative importance
}



# -- Step-wise Regression Method
# 
# apply a linear model and reduce automatically the number of variable using
# a stepwise approach: variables are sequentially entered or removed from
# the model, and their added informative contribution is used to decide
# whether to keep them or not. We use bidirectional elimination testing at
# each step for variables to include or exclude.
#
# Prone to overfitting (i.e. results may be too specific).
#
# If you have large number of predictors, consider splitting the Data in 
# chunks of 10 predictors with each chunk holding the responseVar.
#
imp_stepwise <- function(formula, data, direction="both", tune=F, verbose=T) {

    if (verbose) cat.info("building linear model\n")
    # for the base model we need a formula that looks as "dep.var ~ 1"
    # take LHS of the formula (the dependent variable) and convert (deparse) to string
    y=deparse(f_lhs(formula))
    base.formula <- as.formula(paste(y, "~ 1"))		# create the needed base formula
    base.mod <- lm(base.formula , data=data) 
    if (tune == F) 
        full.mod <- lm(formula, data=data) # full model with all predictors
    else  {
        control <- trainControl(method="repeatedcv", number=10, repeats=3)
        full.mod <- train(formula, data=tc, method='lm', trControl=control)$finalModel
    }

    if (verbose) cat.info("calculating importance\n")
    stepMod <- step(base.mod, 
		            scope = list(lower = base.mod, upper = full.mod), 
                    direction = direction, 
                    trace = 1, 
                    steps = 1000) 

    # according to str() we want to save $Coefficients
    return( stepMod$coefficients )
}
# The output might include levels within categorical variables, since
# 'stepwise' is a linear regression based technique.
# 
# If you have a large number of predictor variables, the above code may need to
# be placed in a loop that will run stepwise on sequential chunks of
# predictors. The shortlisted variables can be accumulated for further analysis
# towards the end of each iteration. This can be very effective method, if you
# want to
# 
# * Be highly selective about discarding valuable predictor variables. 
# 
# * Build multiple models on the response variable.




#
# -- Recursive Feature Elimination RFE Method
#
# Search for a subset of variables using a backward selection (removing
# variables) to find the optimal combination removing the ones with less
# importance according to some metric.
#
#	We tune RFE using caret, controlling the output sizes from 4 to 2^8
# an letting the algorithm select the optimal number. We use a custom
# random forests model.
#
#	Alternatives are linear models (lmFuncs), naive Bayes (nbFuncs), 
# bagged trees (treebagFuncs) and caret-available functions (caretFuncs)
#
#library(caret)

imp_rfe <- function (formula, data, verbose=T) {
    # WE HAVE A FORMULA    y ~ a + b + c ...  or y ~ .
    # get the left hand side of the formula (the dependent variable) as a sub-formula
    fy <- f_lhs(formula)
    # convert to column names
    ny <- all.vars(fy)
    # get the right hand side (the predictor variables) as a sub-formula
    fx <- f_rhs(formula)
    # convert to column names
    nx <- all.vars(fx)
    if ( nx == '.' ) {
       nx <- colnames(data)[ colnames(data) != ny ]
    } 
    # get x and y
    y <- data[ , ny]
    x <- data[ , nx]

    if (verbose) cat.info("building RFE model\n")
    control <- rfeControl(functions=rfFuncs, method="cv", number=10)
    # caution MAGIC NUMBERS HERE
    results <- rfe(x, y, sizes=2^(2:8), rfeControl=control)

    if (verbose) cat.info("calculating importance\n")
    predictors <- predictors(results)
    if (verbose) {
        # summarize the results
        # list the chosen features
        # plot the results
        plot(results, type=c("g", "o"))
    }


    return(predictors)
}



# -- MARS (earth package)
#
# apply a non-linear model using spline curves estimated automatically to
# model automatically non-linear dependencies and interactions between
# variables.
#
# The earth package implements variable importance based on Generalized cross
# validation (GCV), number of subset models the variable occurs (nsubsets) and
# residual sum of squares (RSS).
library(earth)

imp_mars <- function (formula, data, verbose=T) {
    if (verbose) cat.info("building MARS model\n")
    regressor <- earth(formula, data=data) 

    if (verbose) cat.info("calculating importance\n")
    ev <- evimp (regressor) # estimate variable importance

    return(ev)
}



# 
# -- Learning Vector Quantization (LVQ) Method for CATEGORICAL variables
#
# can be considered a special kind of neural network, precursor of SOM and
# related o kNN, that applies a winner takes it all learning approach. During
# learning it tends to select the most meaningful variables.

#library(caret)

imp_lvq <- function(formula, data, verbose=T) {
    if (verbose) cat.info("building LVQ model\n")
    control <- trainControl(method="repeatedcv", number=10, repeats=3)
    # to avoid error 'subscript out of bounds'
    tuneGrid = expand.grid(k=1:3,size=seq(5,20,by=1))
    # train the model
    regressor <- train(formula, data=data, 
                      method="lvq", 			# error: wrong model for regression
                      preProcess="scale", 
                      trControl=control,
                      tuneGrid=tuneGrid)
    if (verbose) cat.info("calculating importance\n")
    # estimate variable importance
    importance <- caret::varImp(regressor, scale=FALSE)$importance

    return(importance)
}



#
# -- rpart (Recursive partitioning and regression trees)
#
# is an extension of decision trees from tables, and can be used for both
# classification, regression and survival like the popular and successful
# CART (Classification and Regression Trees, which builds a binary tree
# of binary or numeric data using the GINI index and cost-complexity pruning) 
# or C5.0 (which gives a binary or multibranch classification tree using
# information gain or entropy pruned by the binomial confidence limit)
#
#library(rpart)
#library(caret)

imp_rpart <- function(formula, data, tune=F, verbose=T) {
    if (tune == F) {
        if (verbose) cat.info("building rpart model\n")
        # this rcontrol was taken from iv_num in package 'riv'
        #	CAUTION: MAGIC NUMBERS
        minbucket <- nrow(data)/10
        if (minbucket < 1) minbucket <- 1
        rcontrol <- rpart.control(cp=0.001, minbucket=minbucket)
        model <- rpart(data=tc, formula=formula, control=rcontrol)
        # we need now a way to calculate importance. woe/riv might
        # be a good one if it did work, so we will use the ones
        # rpart itself calculates
        # model$variable.importance is a named vector of weights with variables as names
        return(model$variable.importance)
    } else {	# (tune == TRUE)
        if (verbose) cat.info("building rpart model with CART\n")
        set.seed(SEED)
        # CAUTION: MAGIC NUMBERS
        trainControl= trainControl(method="repeatedcv", repeats=5)
        rPartMod <- train(formula, 
                          data=data, 
                          method="rpart",
                          importance=T,
                          trControl=trainControl)

        if (verbose) cat.info('calculating importance\n')
        rpartImp <- caret::varImp(rPartMod, scale=F)

        return(rpartImp)
    }
}





# -- random forest
#
# a random forest uses an ensemble learning method to combines the output 
# of multiple decision trees to reach a single result and works for both, 
# classification and regression (using a najority rule for classification
# or the average for regression).
#
# Random forest can be very effective to find a set of predictors that best
# explains the variance in the response variable.
# library(caret)
# library(randomForest)
# This function builds a Random Forest to find the most relevant variables
#
# A Random Forest is like 'rpart' but with many decision trees instead of only one
#
imp_random_forest <- function(formula=NULL, data, dep.var='', verbose=T, method='rangerRF', ...) {
    # create a RF model
    if (verbose) cat.info("building RF model\n")
    if (method == 'rangerRF') {
        # Ranger is a (not so) fast implementation of random forest particularly suited for high
        # dimensional data... 
        # ..and seems to work for <= 5 samples
        if (verbose) cat.info("calculating importance\n")	# we can do both together with ranger
        regressor <- ranger::ranger(formula, data=data, 
                                    importance='permutation', 
                                    local.importance=T,
                                    scale.permutation.importance=T)
        return(regressor$variable.importance)			# importance for each independent variable
        #return(regressor$variable.importance.local)	# importance for each variable and for each sample
 	} else if (method == 'rangerRF.dep.var' ) {
	    # formula is not a formula, but the name of the dependent variable
		# we need to use this when dimensionality is too high and the
		# formula interface fails
		regressor <- ranger::ranger(dependent.variable.name=dep.var, 
		                            data=tc,
									importance='permutation',
									local.importance=T,
									scale.permutation.importance=T)
        return(regressor$variable.importance)			# importance for each independent variable
   }
    # else (if method == 'RF')
    regressor <- randomForest(formula , data=data, importance=TRUE, ...)  # fit the random forest with 
                                                                          # provided parameters (if any)
    if (verbose) cat.info("calculating importance\n")
    # get variable importance, based on mean decrease in accuracy using
    # varImp in the caret package (we indicate it explicitly to avoid name clashes due to masking)
    if ( TRUE ) {
        result <- caret::varImp(regressor, conditional=TRUE) 	# conditional=True, adjusts for correlations 
    													 	    # between predictors
        return(result)
    } else {
        # unused for now, we will activate it eventually and add a function-call argument
        # to allow choosing
        varImp::varimpAUC(regressor)    # more robust towards class imbalance, returns a list
        return(result$importance)
    }
}



# -- Boruta Method
# 
# The 'Boruta' method can be used to decide if a variable is
# important or not based on RF. It is an all relevant method, not a 
# minimal optimal one (may include redundant variables).
#
# Random Ferns, originally proposed for vision can be considered a
# constrained decision tree ensemble, as reliable as randomized trees
# but faster, and compete with SVM in some tasks.
#
# We allow selecting each of the two methods as the base for Boruta.
#

# install.packages('Boruta')
#library(Boruta)

imp_boruta <- function (formula, data, getImp='getImpRfZ', out.png=NULL, verbose=2) {
    if (verbose) cat.info("building boruta model\n")

    boruta_output <- Boruta(formula, data, doTrace=verbose) # perform Boruta search

    # plot results
    if (verbose) 
	    as.png(plot(boruta_output,
		       cex.axis=.7, las=2, xlab="", main="Variable Importance"),
		       file=out.png)
    #as.png(
    #    plot(boruta_output, cex.axis=.7, las=2, xlab="", main="Variable Importance"),
    #    file="impgenes/imp_boruta_tc.png") 


    if (verbose) cat.info("calculating importance\n")
    # get confirmed and tentative factors
    # boruta_signif <- names(boruta_output$finalDecision)[boruta_output$finalDecision 
    # 									      %in% c("Confirmed", "Tentative")]
    # alternate (and more elegant) way
    boruta_signif <- getSelectedAttributes(boruta_output, withTentative=T)
    if (verbose) attStats(boruta_output)

    return(boruta_signif)
}





#
# -- Regularized Random Forest (RRF)
#
# regularized regression estimates a penalization function to produce
# more parsimonious models that have less variables with predictive power 
# similar to full models and robust to correlated variables, whereas RF 
# rely only on tree randomization. RRF can lead to substantially more
# accurate models.
#
imp_rrf <- function ( formula, data, verbose=T ) {
    if (verbose) cat.info("building RRF model\n")

    set.seed(SEED)

    if (verbose) cat.info("building rrf model\n")
    rrfMod <- train(formula, 
                data = data, 
                method = "RRF",
                importance = T) #this argument is required for varImp

    if (verbose) cat.info('calculating importance\n')
    rrfImp <- caret::varImp(rrfMod, scale=F)

    if (verbose) as.png(
				    plot(rrfImp, top = 20, main='Variable Importance'),
                    file="impgenes/imp_rrf_tc.png", verbose)

    if (verbose) plot(rrfImp$importance)

    return(rrfImp$importance)
}






#
#
# -- DALEX (Model Agnostic Language for Exploration and Explanation)
#
# DALEX is a model-agnostic framework to explore any model and 
# explain its behavior calculating the contribution of each 
# variable as changes in the expected model response given other 
# variables. 
#
# DALEX does not seem to work with ranger and categorical variables.
# RandomForest is limited by gene names, and takes much, much longer
# than ranger, but should work.
#
#library(randomForest)
#library(DALEX)
#
imp_dalex <- function(formula=NULL, data, dep.var='', method='rf', verbose=T) {
    if (verbose) cat.info("building DALEX model with", method, "\n")
    if (method == "rf") {
        regressor <- randomForest(formula , data=data, importance=TRUE) 
        # fit the random forest with default parameters
	} else if (method == "rangerRF") {
	    # DO NOT USE
        # Ranger is a (not so) fast implementation of random forest particularly suited for high
        # dimensional data... 
        # NOTE THAT DALEX DOES NOT WORK WITH RANGER
		regressor <- ranger::ranger(formula, data=data, 
                                    importance='permutation', 
                                               local.importance=T,
                                               scale.permutation.importance=T)
	} else if (method == 'rangerRF.dep.var' ) {
	    # DO NOT USE
	    # formula is not a formula, but the name of the dependent variable
		# we need to use this when dimensionality is too high and the
		# formula interface fails
        # NOTE THAT DALEX DOES NOT WORK WITH RANGER
		regressor <- ranger::ranger(dependent.variable.name=dep.var, 
		                            data=data,
									importance='permutation',
									local.importance=T,
									scale.permutation.importance=T)
    } else if (method == "bGLM") {
        regressor <- glm(formula , data=data, family='binomial')
    } else if (method == "mGLM") {
        cat.warn('multinomial models are unsupported\n')
        regressor <- glm(formula , data=data, family='multinomial')
    } else if (method == "gGLM") {
        regressor <- glm(formula , data=data, family='gaussian')
    } else {
        cat.err("unknown method", method, '\n')
    }

    if (verbose) cat.info("calculating importance\n")
    # Variable importance with DALEX
    # we need 'y', and to get it we need to parse the formula:
    #    all.vars returns all the variables that have been used in the formula
    #    so, the first element is the 'y' component (y ~ a ? b ? ...)
    if (dep.var == '') {
    	y <- all.vars(formula)[1]
    } else {
	    y = dep.var
	}
	# and now we use as y the column named y
    explained <- DALEX::explain(regressor, data=data, y=data[ , y], label=method)
    #x <- data[ , ! names(data) %in% c(y) ]
	#explained <- DALEX::explain(regressor, data=x, y=data[ , y], label=method)
    #
    # residuals <- model_performance(explained)
    # plot(residuals)
    # vip <- DALEX::variable_importance(explained) 
    # pdp <- DALEX::variable_response(explained, variable=*select one var*, type='pdp')

    DALEX=T
    if (DALEX) {
        if (verbose) cat.info('extracting variable importance (slow)\n')
        vip <- DALEX::variable_importance(explained)
		return( vip )
    } else {
        if (verbose) cat.info('extracting feature importance (slow)\n')
        # Get the variable importances
        #varimps <- variable_dropout(explained, type='raw')
        #print(varimps)
        #plot(varimps)
        vip <- DALEX::feature_importance(explained, type='variable_importance') # very slow
        
		# plot 20 most important
        if (verbose) {
		    cat.info('plotting importance')
            nv<-dim(vip)[1] / length(levels(as.factor(vip$permutation)))
            as.png(plot(vip[c(1,(nv-20):nv),]), 'impgenes/top20_dalex.png')
        }
        # get mean loss and sort to get variable importance
        vip.mean.loss <- aggregate(vip$dropout_loss, list(vip$variable), mean)
        imp <- (vip.mean.loss[ order(vip.mean.loss$x, decreasing=F),])
        colnames(imp) <- c('variable', 'mean_dropout_loss')
        #print(imp)
        #head(imp)
        #head(sort(imp))
        return ( imp )
    } 
}




#
# -- VITA
#
#		PIPM implements the test approach of Altmann et 
# al. (2010) for the permutation variable importance measure 
# in a random forest for classification and regression by randomly
# shuffling the values of a single variable and observing the
# resulting degradation of the model score.
#
#library(vita)
#
imp_vita <- function (formula=NULL, data, dep.var="", method='RF', verbose=T) {
    if (verbose) cat.info("building RF model\n")

    if (method == 'rangerRF') {
        # Ranger is a (not so) fast implementation of random forest particularly suited for high
        # dimensional data... 
        # ..and seems to work for <= 5 samples
        regressor <- ranger::ranger(formula, data=data, 
                                    importance='permutation', 
                                               local.importance=T,
                                               scale.permutation.importance=T)
		# NOTE that PIPM implements the test approach of Altmann et 
		# al. (2010) for the permutation variable importance measure 
		# ‘VarImp’ in a random forest for classification and regression.
		# When we use ranger, we are already asking for the permutation
		# approach, so this is the same as using 'ramger' in imp_RF
		return(regressor$variable.importance)
    } else if (method == 'rangerRF.dep.var' ) {
	    # formula is not a formula, but the name of the dependent variable
		# we need to use this when dimensionality is too high and the
		# formula interface fails
        regressor <- ranger::ranger(dependent.variable.name=dep.var, 
		                            data=tc,
									importance='permutation',
									local.importance=T,
									scale.permutation.importance=T)
		# NOTE that PIPM implements the test approach of Altmann et 
		# al. (2010) for the permutation variable importance measure 
		# ‘VarImp’ in a random forest for classification and regression.
		# When we use ranger, we are already asking for the permutation
		# approach, so this is the same as using 'ramger' in imp_RF
		return(regressor$variable.importance)
	} else if (method == 'RF') {
        regressor <- randomForest(formula , data=data, importance=TRUE)  # fit the random forest with 
    }

    if (verbose) cat.info("calculating importance\n")
    y <- all.vars(formula)[1]
    pimp.varImp.reg <- PIMP(data, data[ , y], regressor,S =10, parallel=TRUE)
    return(pimp.varImp.reg$VarImp[order(pimp.varImp.reg$VarImp, decreasing=T),])
}




# -- LASSO regression (least absolute shrinkage and selection operator)
#
# selects variables and regularizes models to improve their precision and
# interpretability. Originally formulated for linear regression.
#
# Least Absolute Shrinkage and Selection Operator (LASSO) penalizes 
# with L1-norm.
#     It imposes a cost to having large weights (value of coefficients). 
#	  It is called L1 regularization, because the cost added, is
# proportional to the absolute value of weight coefficients.
# As a result, it reduces the coefficients of unwanted variables to zero leaving
# only the important ones.
#
#library(lasso)
#
# this works for both continuous and categorical variables
imp_lasso <- function (formula, data, verbose=T) {
    # get the left hand side of the formula (the dependent variable) as a sub-formula
    fy <- f_lhs(formula)
    # convert to column names
    ny <- all.vars(fy)
    # get the right hand side (the predictor variables) as a sub-formula
    fx <- f_rhs(formula)
    # convert to column names
    nx <- all.vars(fx)
    if ( nx == '.' ) {
       nx <- colnames(data)[ colnames(data) != ny ]
    } 
    # get x and y
    y <- data[ , ny]
    x <- data[ , nx]


    # make y binomial
    #y <- as.double(as.matrix(ifelse(y == 0.0, 0, 1))) 
    #y <- as.double(as.matrix(y))

    if(verbose) cat.info("building lasso model\n")
    set.seed(SEED)
    if (is.factor(y)) {
        if (length(levels(y)) == 2) {
            family <-'binomial'
        } else {
            family='multinomial'
        }
        if (verbose) cat.info("y is a factor, using a", family, "lasso model\n")
    } else {
        if (verbose) cat.info("assuming gaussian y\n")
        family='gaussian'
    }        
    cv.lasso <- cv.lasso(as.matrix(x), y, family=family, 
    		 		      alpha=1, 			# 0 -> ridge, 1-> LASSO
                          parallel=TRUE, 
                          standardize=TRUE #, 
                          #type.measure='auc'
                         )

    if (verbose) plot(cv.lasso)
    # plot with two vertical lines: first is the lambda with lowers mean squared error
    # and second the number of variables with highest deviance within 1 SD

    if (verbose) cat.info("calculating importance\n")
    # the best lambda value is stored in $lambda,min
    if (family != 'multinomial') {
        coef<- round(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min)), 2)
        coef <- coef[coef[ ,1] !=0, ]
    } else {
        coef <- list()
        for (i in levels(y)) {   # we get a set of coefficients for each level
            level.coef <- round(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min)[[ i ]]), 2)
            coef[[i]] <- level.coef[ level.coef[ ,1] != 0, ]
        }
    }

    # another way
    #coefList <- coef(cv.lasso, s='lambda.1se')
    #coefList <- data.frame(coefList@Dimnames[[1]][coefList@i+1],coefList@x)
    #names(coefList) <- c('var','val')


    # return coefficients with lambda value !=0
    return( coef )

}




# -- xgboost (Extreme Gradient Boosting)
#
# Gradient boosting tries to find the function that best approximates
# the data by taking lots of simple functions and adding them together
# and training on the gradient of the error with respect to loss 
# predictions. In gradient boosting trees we consider as the simple model
# a decision tree. XGBoost looks at the distribution of variables across
# all data points in a leaf to reduce the search space, adds some
# regularizations and speed up the calculations.
#
# Gradient boosting trees can be more accurate than RF and RRF as they 
# can correct each other's errors.
#
#library(caret)
#library(xgboost)
#
imp_xgboost <- function(formula, data, verbose=T) {
    if (verbose) cat.info("building xgboost model\n")

    regressor <- train(formula, data, 
                       method = "xgbTree", 
                       trControl = trainControl("cv", number = 10))
                       # scale=TRUE)
    if (verbose) cat.info("calculating importance\n")
    result <- caret::varImp(regressor, conditional=T)

    # according to str() we want to save $importance from caret::varImp() results
    return(result$importance)
}




# -- Genetic Algorithm
#
# repeatedly uses a huge number of random selections of the variables and 
# calculates the RMSECV (RMS error of cross-validation) and select the best
# candidates discrading lower half combination on each iteration and "breeding"
# and mutating the best half until the ending criteria is met..
#
#library(caret)

imp_ga <- function (formula, data, verbose=T) {
    # get the left hand side of the formula (the dependent variable) as a sub-formula
    fy <- f_lhs(formula)
    # convert to column names
    ny <- all.vars(fy)
    # get the right hand side (the predictor variables) as a sub-formula
    fx <- f_rhs(formula)
    # convert to column names
    nx <- all.vars(fx)
    if ( nx == '.' ) {
       nx <- colnames(data)[ colnames(data) != ny ]
    } 
    # get x and y
    y <- data[ , ny]
    x <- data[ , nx]

    if (verbose) cat.info("building GA model\n")
    # Define control function
    ga_ctrl <- gafsControl(functions = rfGA,  method = "cv", repeats = 3)
    # another option is `caretGA`
    # Genetic Algorithm feature selection
    ga_obj <- gafs(x,
 	              y,
 	              iters = 100,  # normally much higher (100+)        
  	              gafsControl = ga_ctrl)

    return(ga_obj$optVariables)
}




# -- Simulated Annealing
#
# is a global search method that makes small perturbations (removing a
# small subset of model variable members) to an initial candidate solution 
# trying to improve it with an acceptance probability.
# A subotimal solution should eventually improve over interations.
#
#library(caret)
#
imp_sa <- function (formula, data, verbose=T) {
    # get the left hand side of the formula (the dependent variable) as a sub-formula
    fy <- f_lhs(formula)
    # convert to column names
    ny <- all.vars(fy)
    # get the right hand side (the predictor variables) as a sub-formula
    fx <- f_rhs(formula)
    # convert to column names
    nx <- all.vars(fx)
    if ( nx == '.' ) {
       nx <- colnames(data)[ colnames(data) != ny ]
    } 
    # get x and y
    y <- data[ , ny]
    x <- data[ , nx]

    if (verbose) cat.info("building SA model\n")
    # Define control function  
    sa_ctrl <- safsControl(functions = rfSA,
                           method = "repeatedcv",            
                           repeats = 3,            
                           improve = 5)

    set.seed(SEED)
    # Simulated Annealing feature selection
    sa_obj <- safs(x, 
    		       y,
                   safsControl = sa_ctrl)

    return(sa_obj$optVariables)
}


# EFA
#
#	Try to uncover the hidden, underlying structure of a large set of
# variables identifying subsets of related variables.
#
get_efa_factors <- function(data, n.factors, annotation, force.recalculation=F, out.prefix=quote(data), verbose=T) {
    
	efa.save.file <- paste(out.prefix, '.efa.', n.factors, '.Rds', sep='')
    
    # check if already saved and, if so, load precomputed results
    if ( file.exists(efa.save.file) && ( ! force.recalculation ) ) {
        if (verbose) cat.info('loading precomputed EFA from', efa.save.file, '\n')
        data.fa <- readRDS(efa.save.file)
    } else {
        # efa.save.file does not exist or force.recalculation is TRUE
        if (verbose) cat.info('computing EFA.\n')
        # compute correlation matrix
        if (verbose) cat.info('calculating correlation matrix\n')
        data.cor <- cor(data, use="pairwise.complete.obs")
        # calculate EFA
        #    (we should consider possibly passing the extra arguments as '...')
        # 	this can take a significant amount of time for many variables
        if (verbose) cat.info('calculating EFA for', n.factors, 'factors\n')
        # since we are only interested in the hidden factors, we do not need
        # the confidence intervals, and can use fac() instead of fa() to try
        # to speed up the calculations.
        # fac() already uses all the processors available when possible.
        data.fa <- fac(r=data.cor, n.factors=n.factors, 
                     fm="minres", rotate='varimax', scores='tenBerge') #works

        # save results
        #	this can take long and a large amount of space (on the order of GB)
        # on disk
        if (verbose) {
            cat.info('saving results in', efa.save.file, '\n')
            cat.info('this may take a large amount of space.\n')
        }
        saveRDS(data.fa, efa.save.file)

        out <- paste(out.prefix, '.efa.', n.factors, '.txt', sep='')
        sink(out)
        print(data.fa)
        sink()

        # and these two may take also long to create
        if (verbose) cat.info('plotting the results.\n')
        out <- paste(out.prefix, '.efa.', n.factors, '.plot.png', sep='')
        as.png(fa.plot(data.fa), file=out)
        out <- paste(out.prefix, '.efa.', n.factors, '.diagram.png', sep='')
        as.png(fa.diagram(data.fa), file=out)
    }
    
    # the following is quick and, for now, we do not mind repeating each time
    data.fa.clusters <- factor2cluster(data.fa, aslist=F)
    # write a table with values -1, 0, 1 for each variable and factor
    out <-  paste(out.prefix, '.efa.', n.factors, '.clusters.png', sep='')
    write.table(data.fa.clusters, file=out, sep='\t', col.names=T)

    # it may be better to get them as a list and then annotate and save each
	# separately
    data.fa.clusters <- factor2cluster(data.fa, aslist=T)
    for (i in 1:length(data.fa.clusters)) {
    #for (i in seq(along=data.fa.clusters)) {
        name <- names(data.fa.clusters)[i]
        cat('processing', i, name, '\n')
        genes <- data.fa.clusters[[ i ]]
        clus <- data.frame(ensembl_gene_id=genes)
        aclus <- enrich_genes(clus, annotation)
        out <- paste(out.prefix, '.efa.', n.factors, '.', name, '.tab', sep='')
        cat('saving as', out, '\n')
        write.table(aclus, out, sep='\t', row.names=F)
    }
    
    if (F) {
        # These take too long to be feasible for now.
        # Maybe later, we could activate them for sub-clasifying
        # smaller subsets.

        # use a hierarchical clustering solution and find omega
        data.iclust <- iclust(data)
        data.omega <- omega(data)
        data.ba <- bassAckward(data)

        # find alpha as an estimate of reliability
        data.alpha <- alpha(data)
    }
}




######################################################################
# load and prepare all required data
######################################################################

# create dir to save the results
dir.create('impgenes', showWarning=FALSE)
dir.create('efa', showWarning=FALSE)
dir.create('ca', showWarning=FALSE)
dir.create('dca', showWarning=FALSE)


###########################################################################
#
# load annotation
#
# The annotation is in the same location in both cases
#
#ann <- get_annotation('../../../../net.EnsDb/net.ens.db.rds', '../../biomaRt.annotation.1st.txt')
# ens.ann <- ann$ensembl, bio.ann <- ann$biomart
# for now, ensembl does not work during column selection
if (exists(quote(bm.annot.1)))
    bio.ann <- bm.annot.1
else
    bio.ann <- get_annotation(ensembl.db=NULL, biomart.db='../../biomaRt.annotation.1st.txt')$biomart

if (is.null(bio.ann)) cat.err("Could not read BIOMART database\n")


###########################################################################
#
# Analyze
#
# the analyses are different depending on the data we have
# although we might probably fuse both...
#
if (test_coturnix) { 
	ORG <- 'Cjaponica'
	formula <- pfu ~ .

	##########################################################################################
	#
	# Load experiment-specific data
	#
	##########################################################################################

	# for Coturnix japonica we use
	data <- load_coturnix_l2fc_data()
	lfc <- data$signif.l2fc.all.data			# log2FC data (for any pair of differences)
	clfc <- data$signif.l2fc.common.data		# common log2FC present in all samples

	all.signif.genes <- colnames(lfc)					# genes significant in any pairwise comparison
	norm_counts <- load_coturnix_norm_counts_data(
	                             file="../../../normalized_counts.tab", 
								 select=all.signif.genes )
	# get the dataset that we will use
	tc <- norm_counts$signif_transp_norm_counts

	# at this point, in all lfc, clfc, and tc, the 'pfu' column contains the 
	# PFUs used to infect the cells, 0 for wt

	# prepare multinomial factor versions of the data (pfu converted to factor) for use when
	# regression fails or with categorical methods
	fdat <- lfc ; fdat$pfu = as.factor(fdat$pfu)
	fdatz <- clfc ; fdatz$pfu <- as.factor(fdatz$pfu)
	ftc <- tc ; ftc$pfu <- as.factor(ftc$pfu)
	# binomial version (just in case) for binomial methods (this should identify wt-specific genes)
	btc <- ftc
	btc[ ,1] = as.factor(ifelse(btc[ , 1] == 0, 0, 1))



	##########################################################################################
	#
	# Now, we can proceed with the analysis of gene importance:
	# For Coturnix data (outcome is numeric)
	#	we expect to identify the most important genes to "explain" the PFU level (i.e., the
	# genes that are more relevantly associated with the response to different PFU infection
	# levels); they will likely represent constallations of related genes *whose association
	# with the outcome is similar*.
	#
	##########################################################################################

    # 1: correlation
	#
    lfc.i.cor <- imp_cor(formula, lfc)
    clfc.i.cor <- imp_cor(formula, clfc)
    tc.i.cor <- imp_cor(formula, tc)

    # we need to tune the following for each method because each
	# returns the results in a different way
	# (it would be better to ensure all functions return the same
	# data columns in the same format)
    if (exists(quote(lfc.i.cor)) && ! is.null(lfc.i.cor)) {
        edat.i.cor <- data.frame(ensembl_gene_id=(lfc.i.cor))
        edat.i.cor <- enrich_genes(edat.i.cor, bio.ann)
        write.table(edat.i.cor, 'impgenes/imp_cor_dat.txt', sep='\t', row.names=F)
    }
    if (exists(quote(clfc.i.cor)) && ! is.null(clfc.i.cor)) {
        eclfc.i.cor <-data.frame(ensembl_gene_id=clfc.i.cor)
        eclfc.i.cor <- enrich_genes(eclfc.i.cor, bio.ann)
        write.table(eclfc.i.cor, 'impgenes/imp_cor_datz.txt', sep='\t', row.names=F)
    }
    if (exists(quote(tc.i.cor)) && ! is.null(tc.i.cor)) {
        etc.i.cor <- data.frame(ensembl_gene_id=tc.i.cor)
        etc.i.cor <- enrich_genes(etc.i.cor, bio.ann)
        write.table(etc.i.cor, 'impgenes/imp_cor_tc.txt', sep='\t', row.names=F)
    }

    # 2: relimp
	#
    lfc.i.relimp <- imp_relimp(formula, lfc, tune=T)
    clfc.i.relimp <- imp_relimp(formula, clfc, tune=T)
    tc.i.relimp <- imp_relimp(formula, tc, tune=T)			# only one gene selected "ENSCJPG00005008112"

    if (! is.null(lfc.i.relimp)) {
        edat.i.relimp <- cbind(ensembl_gene_id=rownames(lfc.i.relimp), lfc.i.relimp)
        edat.i.relimp <- enrich_genes(edat.i.relimp, bio.ann)
        write.table(edat.i.relimp, 'impgenes/imp_relimp_dat.txt', sep='\t', row.names=F)
    }
    if (! is.null(clfc.i.relimp)) {
        eclfc.i.relimp <- cbind(ensembl_gene_id=rownames(eclfc.i.relimp), eclfc.i.relimp)
        eclfc.i.relimp <- enrich_genes(eclfc.i.relimp, bio.ann)
        write.table(eclfc.i.relimp, 'impgenes/imp_relimp_datz.txt', sep='\t', row.names=F)
    }
    if (! is.null(tc.i.relimp)) {
        etc.i.relimp <- cbind(ensembl_gene_id=rownames(tc.i.relimp), tc.i.relimp)
        etc.i.relimp <- enrich_genes(etc.i.relimp, bio.ann)
        write.table(etc.i.relimp, 'impgenes/imp_relimp_tc.txt', sep='\t', row.names=F)
    }

    # 3: stepwise
	#
    lfc.i.step <- imp_stepwise(formula, lfc)		# NA
    clfc.i.step <- imp_stepwise(formula, clfc)		# NA
    tc.i.step <- imp_stepwise(formula, tc)			# NA

    if (! is.null(lfc.i.step)) {
        edat.i.step <- data.frame(ensembl_gene_id=names(lfc.i.step), cpeffocient=lfc.i.step)
        edat.i.step <- enrich_genes(edat.i.step, bio.ann)
        write.table(edat.i.step, 'impgenes/imp_step_dat.txt', sep='\t', row.names=F)
    }
    if (! is.null(clfc.i.step)) {
        eclfc.i.step <-data.frame(ensembl_gene_id=names(clfc.i.step), cpeffocient=clfc.i.step)
        eclfc.i.step <- enrich_genes(eclfc.i.step, bio.ann)
        write.table(eclfc.i.step, 'impgenes/imp_step_datz.txt', sep='\t', row.names=F)
    }
    if (! is.null(tc.i.step)) {
        etc.i.step <- data.frame(ensembl_gene_id=names(tc.i.step), cpeffocient=tc.i.step)
        etc.i.step <- enrich_genes(etc.i.step, bio.ann)
        write.table(etc.i.step, 'impgenes/imp_step_tc.txt', sep='\t', row.names=F)
    }

    # 4: recursive feature elimination
	#
    fdat <- lfc ; fdat$pfu <- as.factor(fdat$pfu)
    fdatz <- clfc ; fdatz$pfu <- as.factor(fdatz$pfu)
    ftc <- tc ; ftc$pfu <- as.factor(ftc$pfu)

    lfc.i.rfe <- imp_rfe(formula, fdat)
    clfc.i.rfe <- imp_rfe(formula, fdatz)
    tc.i.rfe <- imp_rfe(formula, ftc)

    if (! is.null(lfc.i.rfe)) {
        edat.i.rfe <- data.frame(ensembl_gene_id=lfc.i.rfe)
        edat.i.rfe <- enrich_genes(edat.i.rfe, bio.ann)
        write.table(edat.i.rfe, 'impgenes/imp_rfe_dat.txt', sep='\t', row.names=F)
    }
    if (! is.null(clfc.i.rfe)) {
        eclfc.i.rfe <- data.frame(ensembl_gene_id=clfc.i.rfe)
        eclfc.i.rfe <- enrich_genes(eclfc.i.rfe, bio.ann)
        write.table(eclfc.i.rfe, 'impgenes/imp_rfe_datz.txt', sep='\t', row.names=F)
    }
    if (! is.null(tc.i.rfe)) {
        etc.i.rfe <- data.frame(ensembl_gene_id=tc.i.rfe)
        etc.i.rfe <- enrich_genes(etc.i.rfe, bio.ann)
        write.table(etc.i.rfe, 'impgenes/imp_rfe_tc.txt', sep='\t', row.names=F)
    }

    # 5: mars
	#
    lfc.i.mars <- imp_mars(formula, data=lfc)		# fails because of NAs 	ENSCJPG00005000039
    clfc.i.mars <- imp_mars(formula, data=clfc)		# same					ENSCJPG00005000039
    tc.i.mars <- imp_mars(formula, data=tc)			# only one important gene: ENSCJPG00005008112

    if (! is.null(lfc.i.mars)) {
        edat.i.mars <- cbind(ensembl_gene_id=rownames(lfc.i.mars), lfc.i.mars)
        edat.i.mars <- enrich_genes(edat.i.mars, bio.ann)
        write.table(edat.i.mars, 'impgenes/imp_mars_dat.txt', sep='\t', row.names=F)
    }
    if (! is.null(clfc.i.mars)) {
        eclfc.i.mars <- cbind(ensembl_gene_id=rownames(eclfc.i.mars), eclfc.i.mars)
        eclfc.i.mars <- enrich_genes(eclfc.i.mars, bio.ann)
        write.table(eclfc.i.mars, 'impgenes/imp_mars_datz.txt', sep='\t', row.names=F)
    }
    if (! is.null(tc.i.mars)) {
        etc.i.mars <- cbind(ensembl_gene_id=rownames(tc.i.mars), tc.i.mars)
        etc.i.mars <- enrich_genes(etc.i.mars, bio.ann)
        write.table(etc.i.mars, 'impgenes/imp_mars_tc.txt', sep='\t', row.names=F)
    }

    # 6: lvq
	#
	# for tc only one gene was selected
	#	ENSCJPG00005008112
    fdat <- lfc ; fdat$pfu <- as.factor(fdat$pfu)
    fdatz <- clfc ; fdatz$pfu <- as.factor(fdatz$pfu)
    ftc <- tc ; ftc$pfu <- as.factor(ftc$pfu)

    lfc.i.lvq <- imp_lvq(formula, fdat)
    clfc.i.lvq <- imp_lvq(formula, fdatz)
    tc.i.lvq <- imp_lvq(formula, ftc)

    if (! is.null(lfc.i.lvq)) {
        edat.i.lvq <- cbind(ensembl_gene_id=rownames(lfc.i.lvq), lfc.i.lvq)
        edat.i.lvq <- enrich_genes(edat.i.lvq, bio.ann)
        write.table(edat.i.lvq, 'impgenes/imp_lvq_dat.txt', sep='\t', row.names=F)
    }
    if (! is.null(clfc.i.lvq)) {
        eclfc.i.lvq <- cbind(ensembl_gene_id=rownames(eclfc.i.lvq), eclfc.i.lvq)
        eclfc.i.lvq <- enrich_genes(eclfc.i.lvq, bio.ann)
        write.table(eclfc.i.lvq, 'impgenes/imp_lvq_datz.txt', sep='\t', row.names=F)
    }
    if (! is.null(tc.i.lvq)) {
        etc.i.lvq <- cbind(ensembl_gene_id=rownames(tc.i.lvq), tc.i.lvq)
        etc.i.lvq <- enrich_genes(etc.i.lvq, bio.ann)
        write.table(etc.i.lvq, 'impgenes/imp_lvq_tc.txt', sep='\t', row.names=F)
    }


    # 7: Boruta
	#
    lfc.i.boruta <- imp_boruta(formula, lfc)
    clfc.i.boruta <- imp_boruta(formula, clfc)
    tc.i.boruta <- imp_boruta(formula, tc)

    # repeat again to clean up the results for plotting
	# otherwise we get too much data and cannot see anything
    selected_columns_dat <- cbind(pfu=lfc[ ,1], lfc[, colnames(lfc) %in% boruta_signif])
    selected_columns_datz <- cbind(pfu=clfc[ ,1], clfc[, colnames(clfc) %in% boruta_signif])
    selected_columns_tc <- cbind(pfu=tc[ ,1], tc[, colnames(tc) %in% tc.i.boruta])

    selected.lfc.i.boruta <- imp_boruta(formula, selected_columns_dat)
    selected.clfc.i.boruta <- imp_boruta(formula, selected_columns_datz)
    selected.tc.i.boruta <- imp_boruta(formula, selected_columns_tc)

    as.png(
        plot(boruta_signif, cex.axis=.7, las=2, xlab="", main="Variable Importance")
        "impgenes/imp_selected_tc.png")

    if (! is.null(lfc.i.boruta)) {
        edat.i.boruta <- data.frame(ensembl_gene_id=lfc.i.boruta)
        edat.i.boruta <- enrich_genes(edat.i.boruta, bio.ann)
        write.table(edat.i.boruta, 'impgenes/imp_boruta_dat.txt', sep='\t', row.names=F)
    }
    if (! is.null(clfc.i.boruta)) {
        eclfc.i.boruta <- data.frame(ensembl_gene_id=clfc.i.boruta)
        eclfc.i.boruta <- enrich_genes(eclfc.i.boruta, bio.ann)
        write.table(eclfc.i.boruta, 'impgenes/imp_boruta_datz.txt', sep='\t', row.names=F)
    }
    if (! is.null(tc.i.boruta)) {
        etc.i.boruta <- data.frame(ensembl_gene_id=tc.i.boruta)
        etc.i.boruta <- enrich_genes(etc.i.boruta, bio.ann)
        write.table(etc.i.boruta, 'impgenes/imp_boruta_tc.txt', sep='\t', row.names=F)
    }
	

    # 8: ferns
	#
    lfc.i.boruta.ferns <- imp_boruta(formula, lfc, getImp='getImpFerns')
    clfc.i.boruta.ferns <- imp_boruta(formula, clfc, getImp='getImpFerns')
    tc.i.boruta.ferns <- imp_boruta(formula, tc, getImp='getImpFerns')


    if (! is.null(lfc.i.boruta.ferns)) {
        edat.i.boruta.ferns <- data.frame(ensembl_gene_id=lfc.i.boruta.ferns)
        edat.i.boruta.ferns <- enrich_genes(edat.i.boruta.ferns, bio.ann)
        write.table(edat.i.boruta.ferns, 'impgenes/imp_boruta.ferns_dat.txt', sep='\t', row.names=F)
    }
    if (! is.null(clfc.i.boruta.ferns)) {
        eclfc.i.boruta.ferns <- data.frame(ensembl_gene_id=clfc.i.boruta.ferns)
        eclfc.i.boruta.ferns <- enrich_genes(eclfc.i.boruta.ferns, bio.ann)
        write.table(eclfc.i.boruta.ferns, 'impgenes/imp_boruta.ferns_datz.txt', sep='\t', row.names=F)
    }
    if (! is.null(tc.i.boruta.ferns)) {
        etc.i.boruta.ferns <- data.frame(ensembl_gene_id=tc.i.boruta.ferns)
        etc.i.boruta.ferns <- enrich_genes(etc.i.boruta.ferns, bio.ann)
        write.table(etc.i.boruta.ferns, 'impgenes/imp_boruta.ferns_tc.txt', sep='\t', row.names=F)
    }

    
	# 9: rpart
	#
    # with lfc: Something is wrong; all the RMSE metric values are missing
    # with clfc: Something is wrong; all the RMSE metric values are missing Model fit failed Argument importance not matched
    # with tc: Something is wrong; all the RMSE metric values are missing Model fit failed Argument importance not matched
    # with ftc: Something is wrong; all the Accuracy metric values are missing: etc.
    lfc.i.rpart <- imp_rpart(formula , data=lfc)		# with tune=TRUE gives ERROR: RMSD values missing
    clfc.i.rpart <- imp_rpart(formula , data=clfc)		# with tune=TRUE gives ERROR: RMSD values missing
    tc.i.rpart <- imp_rpart(formula , data=tc)
    # tc.i.rpart <- imp_rpart(formula , data=ftc)

    if (! is.null(lfc.i.rpart)) {
        edat.i.rpart <- data.frame(ensembl_gene_id=names(lfc.i.rpart), lfc.i.rpart)
        edat.i.rpart <- enrich_genes(edat.i.rpart, bio.ann)
        write.table(edat.i.rpart, 'impgenes/imp_rpart_dat.txt', sep='\t', row.names=F)
    }
    if (! is.null(clfc.i.rpart)) {
        eclfc.i.rpart <- data.frame(ensembl_gene_id=names(clfc.i.rpart), clfc.i.rpart)
        eclfc.i.rpart <- enrich_genes(eclfc.i.rpart, bio.ann)
        write.table(eclfc.i.rpart, 'impgenes/imp_rpart_datz.txt', sep='\t', row.names=F)
    }
    if (! is.null(tc.i.rpart)) {
        etc.i.rpart <- data.frame(ensembl_gene_id=names(tc.i.rpart), tc.i.rpart)
        etc.i.rpart <- enrich_genes(etc.i.rpart, bio.ann)
        write.table(etc.i.rpart, 'impgenes/imp_rpart_tc.txt', sep='\t', row.names=F)
    }
	

    # 10: random forest
	#
	# the 'RandomForest' package implementation fails to process our data,
	# but the 'ranger' package one (designed for high dimensionality) does
	# work: we have made it the default now
	#
    # lfc can have NAs, the others don't
    #lfc.i.rf <- imp_random_forest(formula, lfc, na.action=na.omit)
    clfc.i.rf <- imp_random_forest(formula, clfc)
    tc.i.rf <- imp_random_forest(formula, tc)

    # all of these fail with method='RF' because we only have five infection levels
    #	==> according to Google searches we cannot use RF for regression with <= 5 levels
    # this means we cannot use regression, and we need to resort to 
    # classification, treating the infection levels as a factor
    # (this is, we treat PFU levels as categorical variables, which
    # implies we will ignore the magnitude of $pfu):
    #	convert $pfu to factor and try again
    # now do the (categorical) importance calculations
    fdat.i.rf <- imp_random_forest(formula, fdat, na.action=na.omit)		# need at least two classes
    fdatz.i.rf <- imp_random_forest(formula, fdatz)						# 
    ftc.i.rf <- imp_random_forest(formula, ftc)							# 

    # for rangerRF (for RF change names() by rownames()
    if (exists(quote(lfc.i.rf)) && ! is.null(lfc.i.rf)) {
        edat.i.rf <- cbind(ensembl_gene_id=names(lfc.i.rf), lfc.i.rf)
        edat.i.rf <- enrich_genes(edat.i.rf, bio.ann)
        write.table(edat.i.rf, 'impgenes/imp_rf_dat.txt', sep='\t', row.names=F)
    }
    if (exists(quote(clfc.i.rf)) && ! is.null(clfc.i.rf)) {
        eclfc.i.rf <- cbind(ensembl_gene_id=names(clfc.i.rf), clfc.i.rf)
        eclfc.i.rf <- enrich_genes(eclfc.i.rf, bio.ann)
        write.table(eclfc.i.rf, 'impgenes/imp_rf_datz.txt', sep='\t', row.names=F)
    }
    if (exists(quote(tc.i.rf)) && ! is.null(tc.i.rf)) {
        etc.i.rf <- cbind(ensembl_gene_id=names(tc.i.rf), tc.i.rf)
        etc.i.rf <- enrich_genes(etc.i.rf, bio.ann)
        write.table(etc.i.rf, 'impgenes/imp_rf_tc.txt', sep='\t', row.names=F)
    }

    # 11: regularized random forest
	#
    #lfc.i.rrf <- imp_rrf(formula , data=lfc)		# ERROR: RMSD values missing
    #clfc.i.rrf <- imp_rrf(formula , data=clfc)		# same
    tc.i.rrf <- imp_rrf(formula , data=tc)

    if (exists(quote(lfc.i.rrf)) && ! is.null(lfc.i.rrf)) {
        edat.i.rrf <- cbind(ensembl_gene_id=rownames(lfc.i.rrf), lfc.i.rrf)
        edat.i.rrf <- enrich_genes(edat.i.rrf, bio.ann)
        write.table(edat.i.rrf, 'impgenes/imp_rrf_dat.txt', sep='\t', row.names=F)
    }
    if (exists(quote(clfc.i.rrf)) && ! is.null(clfc.i.rrf)) {
        eclfc.i.rrf <- cbind(ensembl_gene_id=rownames(clfc.i.rrf), clfc.i.rrf)
        eclfc.i.rrf <- enrich_genes(eclfc.i.rrf, bio.ann)
        write.table(eclfc.i.rrf, 'impgenes/imp_rrf_datz.txt', sep='\t', row.names=F)
    }
    if (exists(quote(tc.i.rrf)) && ! is.null(tc.i.rrf)) {
        etc.i.rrf <- cbind(ensembl_gene_id=rownames(tc.i.rrf), tc.i.rrf)
        etc.i.rrf <- enrich_genes(etc.i.rrf, bio.ann)
        write.table(etc.i.rrf, 'impgenes/imp_rrf_tc.txt', sep='\t', row.names=F)
    }


    # 12: DALEX
	#
    # RF does not work for our dataset, so this selection method
    # will neither, unless we use rangerRF which does work.

    # If we use an RF model, since we have less than six values, we need
    # to use a categorical RF
    fdat <- lfc ; fdat$pfu <- as.factor(fdat$pfu)
    fdatz <- clfc ; fdatz$pfu <- as.factor(fdatz$pfu)
    ftc <- tc ; ftc$pfu <- as.factor(ftc$pfu)

    # but if we use rangerRF, it works
    # and so it does with gGLM for numerical or binary data
    fdat.i.dalex <- imp_dalex(formula, fdat)
    fdatz.i.dalex <- imp_dalex(formula, fdatz)
    ftc.i.dalex <- imp_dalex(formula, ftc)

    lfc.i.dalex <- imp_dalex(formula, lfc, method='rangerRF')
    clfc.i.dalex <- imp_dalex(formula, clfc, method='rangerRF')
    tc.i.dalex <- imp_dalex(formula, tc, method='rangerRF')


    if (exists(quote(lfc.i.dalex)) && ! is.null(lfc.i.dalex)) {
        edat.i.dalex <- lfc.i.dalex ; colnames(edat.i.dalex) <- c('ensembl_gene_id', 'mean_dropout_loss')
        edat.i.dalex <- enrich_genes(edat.i.dalex, bio.ann)
        write.table(edat.i.dalex, 'impgenes/imp_dalex_dat.txt', sep='\t', row.names=F)
    }
    if (exists(quote(clfc.i.dalex)) && ! is.null(clfc.i.dalex)) {
        eclfc.i.dalex <- clfc.i.dalex ; colnames(eclfc.i.dalex) <- c('ensembl_gene_id', 'mean_dropout_loss')
        eclfc.i.dalex <- enrich_genes(eclfc.i.dalex, bio.ann)
        write.table(eclfc.i.dalex, 'impgenes/imp_dalex_datz.txt', sep='\t', row.names=F)
    }
    if (exists(quote(tc.i.dalex)) && ! is.null(tc.i.dalex)) {
        etc.i.dalex <- tc.i.dalex ; colnames(etc.i.dalex) <- c('ensembl_gene_id', 'mean_dropout_loss')
        etc.i.dalex <- enrich_genes(etc.i.dalex, bio.ann)
        write.table(etc.i.dalex, 'impgenes/imp_dalex_tc.txt', sep='\t', row.names=F)
    }


    # 13: vita
	#
    lfc.i.vita <- imp_vita(formula , data=fdat)
    clfc.i.vita <- imp_vita(formula , data=fdatz)
    tc.i.vita <- imp_vita(formula , data=ftc)

    if (exists(quote(lfc.i.vita)) && ! is.null(lfc.i.vita)) {
        edat.i.vita <- cbind(ensembl_gene_id=rownames(t(lfc.i.vita)), t(lfc.i.vita))
        edat.i.vita <- enrich_genes(edat.i.vita, bio.ann)
        write.table(edat.i.vita, 'impgenes/imp_vita_dat.txt', sep='\t', row.names=F)
    }
    if (exists(quote(clfc.i.vita)) && ! is.null(clfc.i.vita)) {
        eclfc.i.vita <- cbind(ensembl_gene_id=names(t(clfc.i.vita)), clfc.i.vita)
        eclfc.i.vita <- enrich_genes(eclfc.i.vita, bio.ann)
        write.table(eclfc.i.vita, 'impgenes/imp_vita_datz.txt', sep='\t', row.names=F)
    }
    if (exists(quote(tc.i.vita)) && ! is.null(tc.i.vita)) {
        etc.i.vita <- cbind(ensembl_gene_id=names(tc.i.vita), tc.i.vita)
        etc.i.vita <- enrich_genes(etc.i.vita, bio.ann)
        write.table(etc.i.vita, 'impgenes/imp_vita_tc.txt', sep='\t', row.names=F)
    }

    # 14: LASSO
	#
    lfc.i.lasso <- imp_lasso(formula , data=lfc)
    clfc.i.lasso <- imp_lasso(formula , data=clfc)
    tc.i.lasso <- imp_lasso(formula , data=tc)
    ftc <- tc ; ftc$pfu <- as.factor(ftc$pfu)
    ftc.i.lasso <- imp_lasso(formula , data=ftc)
    btc <- ftc
    btc[ ,1] = as.factor(ifelse(btc[ , 1] == 0, 0, 1))
    btc.i.lasso <- imp_lasso(formula , data=btc)

    if (exists(quote(lfc.i.lasso)) && ! is.null(lfc.i.lasso)) {
        edat.i.lasso <- cbind(ensembl_gene_id=names(lfc.i.lasso), lfc.i.lasso)
        edat.i.lasso <- enrich_genes(edat.i.lasso, bio.ann)
        write.table(edat.i.lasso, 'impgenes/imp_lasso_dat.txt', sep='\t', row.names=F)
    }
    if (exists(quote(clfc.i.lasso)) && ! is.null(clfc.i.lasso)) {
        eclfc.i.lasso <- cbind(ensembl_gene_id=names(clfc.i.lasso), clfc.i.lasso)
        eclfc.i.lasso <- enrich_genes(eclfc.i.lasso, bio.ann)
        write.table(eclfc.i.lasso, 'impgenes/imp_lasso_datz.txt', sep='\t', row.names=F)
    }
    if (exists(quote(tc.i.lasso)) && ! is.null(tc.i.lasso)) {
        etc.i.lasso <- cbind(ensembl_gene_id=names(tc.i.lasso), tc.i.lasso)
        etc.i.lasso <- enrich_genes(etc.i.lasso, bio.ann)
        write.table(etc.i.lasso, 'impgenes/imp_lasso_tc.txt', sep='\t', row.names=F)
    }


    # 15: xgboost
	#
	lfc.i.xgboost <- imp_xgboost(formula, lfc)
    clfc.i.xgboost <- imp_xgboost(formula, clfc)
    tc.i.xgboost <- imp_xgboost(formula, tc)

    if (! is.null(lfc.i.xgboost)) {
        edat.i.xgboost <- cbind(ensembl_gene_id=rownames(lfc.i.xgboost), lfc.i.xgboost)
        edat.i.xgboost <- enrich_genes(edat.i.xgboost, bio.ann)
        write.table(edat.i.xgboost, 'impgenes/imp_xgboost_dat.txt', sep='\t', row.names=F)
    }
    if (! is.null(clfc.i.xgboost)) {
        eclfc.i.xgboost <- cbind(ensembl_gene_id=rownames(eclfc.i.xgboost), eclfc.i.xgboost)
        eclfc.i.xgboost <- enrich_genes(eclfc.i.xgboost, bio.ann)
        write.table(eclfc.i.xgboost, 'impgenes/imp_xgboost_datz.txt', sep='\t', row.names=F)
    }
    if (! is.null(tc.i.xgboost)) {
        etc.i.xgboost <- cbind(ensembl_gene_id=rownames(tc.i.xgboost), tc.i.xgboost)
        etc.i.xgboost <- enrich_genes(etc.i.xgboost, bio.ann)
        write.table(etc.i.xgboost, 'impgenes/imp_xgboost_tc.txt', sep='\t', row.names=F)
    }


    # 16: genetic algorithm
	#
    lfc.i.ga <- imp_ga(formula, lfc)
    clfc.i.ga <- imp_ga(formula, clfc)
    tc.i.ga <- imp_ga(formula, tc)

    if (exists(quote(lfc.i.ga)) && ! is.null(lfc.i.ga)) {
        edat.i.ga <- data.frame(ensembl_gene_id=(lfc.i.ga))
        edat.i.ga <- enrich_genes(edat.i.ga, bio.ann)
        write.table(edat.i.ga, 'impgenes/imp_ga_dat.txt', sep='\t', row.names=F)
    }
    if (exists(quote(clfc.i.ga)) && ! is.null(clfc.i.ga)) {
        eclfc.i.ga <-data.frame(ensembl_gene_id=clfc.i.ga)
        eclfc.i.ga <- enrich_genes(eclfc.i.ga, bio.ann)
        write.table(eclfc.i.ga, 'impgenes/imp_ga_datz.txt', sep='\t', row.names=F)
    }
    if (exists(quote(tc.i.ga)) && ! is.null(tc.i.ga)) {
        etc.i.ga <- data.frame(ensembl_gene_id=tc.i.ga)
        etc.i.ga <- enrich_genes(etc.i.ga, bio.ann)
        write.table(etc.i.ga, 'impgenes/imp_ga_tc.txt', sep='\t', row.names=F)
    }

    
	# 17: simulated annealing
	#
    lfc.i.sa <- imp_sa(formula, lfc)
    clfc.i.sa <- imp_sa(formula, clfc)
    tc.i.sa <- imp_sa(formula, tc)


    if (exists(quote(lfc.i.sa)) && ! is.null(lfc.i.sa)) {
        edat.i.sa <- data.frame(ensembl_gene_id=(lfc.i.sa))
        edat.i.sa <- enrich_genes(edat.i.sa, bio.ann)
        write.table(edat.i.sa, 'impgenes/imp_sa_dat.txt', sep='\t', row.names=F)
    }
    if (exists(quote(clfc.i.sa)) && ! is.null(clfc.i.sa)) {
        eclfc.i.sa <-data.frame(ensembl_gene_id=clfc.i.sa)
        eclfc.i.sa <- enrich_genes(eclfc.i.sa, bio.ann)
        write.table(eclfc.i.sa, 'impgenes/imp_sa_datz.txt', sep='\t', row.names=F)
    }
    if (exists(quote(tc.i.sa)) && ! is.null(tc.i.sa)) {
        etc.i.sa <- data.frame(ensembl_gene_id=tc.i.sa)
        etc.i.sa <- enrich_genes(etc.i.sa, bio.ann)
        write.table(etc.i.sa, 'impgenes/imp_sa_tc.txt', sep='\t', row.names=F)
    }



	##########################################################################################
	#
	# Now, we can proceed with the factor analysis:
	#	we expect to identify the most important groups of genes that determine the
	# outcome variable
	#
	##########################################################################################
	#
	# EFA
	#
	# with EFA, we calculate which genes behave coordinately across
	# the full dataset

	# NOTES:
	#
	#	After trying various packages unsuccessfully, I tried to repeat
	#  EFA using the normalized corrected counts. 
	#
	# With all data it takes too long. So, I tried using smaller subsets.
	# with 100 genes it worked for psych::fa (other implementations failed
	# at singular matrix steps, psych::fa detects the situation and takes
	# a detour (uses the pseudo-inverse) to complete the calculations.
	#
	# Trying with increasing numbers of genes (200, 400, 800) also completed
	# but took much longer
	#
	# This led to two observations:
	#	1) we do not need all genes, only ones with significant l2FC
	#	2) in the cases where fa.parallel was allowed to run to completion
	#	   the optimal number of hidden variables seems to be either 4 (by
	#      the intersection with simulated data) or 6 (by the elbow method)
	#
	# Applying EFA with psych::fa() to the significant l2FC dataset using
	# 6 as the number of factors works more or less... depending on the
	# method: some methods will result in unreliable factoring, some in
	# calculation problems, and some simply work. Increasing the number of
	# factors reaches a point where the calculation fails.
	#
	# This led to additional observations:
	#   1) choosing an appropriate number of factors will determine completion
	#	2) choosing an appropriate combination of methods will also determine
	#      completion.
	#
	# After this, I decided to retake the calculations with the log2FC data
	# directly again applying what I had learnt from the count data and now
	# we can make it work.
	#
	# The plots of hidden factors suggest that there is a major one that
	# includes the (0, 0.1, 1, 10, 100) column. These would be the genes
	# that correlate linearly with infection level. The other groups are
	# genes with coherent responses that do not follow the infection level
	# pattern. They may still be highly relevant and biologically significant
	# but we need to identify what they do.
	#
	#	Conclusion: we need now to separate the factor components and annotate
	# them.
	#
	# Another element to investigate is what happens if we isolate the group
	# of genes that correlate with infection level and run EFA again on them.
	# With a bit of luck, after removing extra "noise" we might be able to
	# further isolate additional subgroups.
	#
	# These are the tries of the analysis that worked (removing the partial analyses):
	#
	#

	WE_HAVE_LOTS_OF_TIME_AND_CPU <- FALSE
	# NOTE we have 100 CPUs, 1TB memory, and still set this to FALSE!!!
	if ( WE_HAVE_LOTS_OF_TIME_AND_CPU ) {
    	# try to find out optimal number of hidden factors using fa.parallel()

    	tc.cor <- cor(tc, use="pairwise.complete.obs")
    	as.png(
    		corrplot(tc.cor, order='hclust', tl.col="black"), #the plot is wrong
    		"efa/tc_corplot.png", width = 800, height = 600)

    	# this may take very long: it will run EFA many times, increasing the number
    	# of factors to find an optimal number.
    	# Maybe we should consider doing it with a subsample? Would that be representative?
    	tc.fa.paral <- fa.parallel(x=tc.cor, fm="minres", fa="fa")
    	print(tc.fa.paral)
	} else {
	    # 800 is on the verge of tolerability for just finding out
		# a prediction for the number of factors/components
		# doing it with different size subsets will give us an 
		# idea of a potential number and its dependence on size
		# (if we assume genes are at random, size should have little
		# to no impact or at least tend to stabilize at some point)
    	tc100.fa.paral <- fa.parallel(x=tc[ , 1:100], fm="minres", fa="fa")
    	print(tc100.fa.paral)
    	tc200.fa.paral <- fa.parallel(x=tc[ , 1:200], fm="minres", fa="fa")
    	print(tc200.fa.paral)
    	tc400.fa.paral <- fa.parallel(x=tc[ , 1:400], fm="minres", fa="fa")
    	print(tc400.fa.paral)
	    	tc800.fa.paral <- fa.parallel(x=tc[ , 1:800], fm="minres", fa="fa")
    	print(tc800.fa.paral)

	}
	# when run incrementally, there seems to be an agreement that
	# the optimal number of factors should be 4 according to fa.parallel()
	# and 6 by the elbow criterion:

	# nFactors = 6
	st.seed(SEED)
	get_efa_factors(tc, n.factors=6, annotation=bio.ann, 
	                force.recalculation=F, 
					out.prefix='efa/tc', 
					verbose=TRUE)
	# this works

	# nFactors = 4
	st.seed(SEED)
	get_efa_factors(data=tc, n.factors=4, annotation=bio.ann force.recalculation=F, out.prefix='efa/tc', verbose=T)
	# this one raises an error:
	#	Error in e$vectors %*% diag(inv.sqrt.ev) : non-conformable arguments

	### EGA
	#
	#	EGA failed.
}





if (test_gallus2) {
	ORG <- 'Ggallus2'
	formula <- sample ~ .

	##########################################################################################
	#
	# Load experiment-specific data
	#
	##########################################################################################

	# for Gallus gallus we'd use something else
	# but we better use a different script until we can identify commonalities
	signif <- get_gallus_significant_genes()
	nc <- load_gallus_norm_count_data(file="normalized_counts.tab",
                            	select=signif$sig.genes)
	snc <- nc$significant.counts
	samples <- colnames(snc)
	samples <- sub('_sample.*', '', samples)
	tc <- data.frame(sample=as.factor(samples))
	tc <- cbind(tc, t(snc))
	# and now 'tc' contains genes in columns, and its first column are sample
	# names, ready for categorical variable reduction
    #tc$samples <- as.factor(tc$samples)

    # specific row subsets from tc (in case we want specific comparisons)
    uninfected <- c(1:3, 7:9, 13:15)
	infected <- c(4:6, 10.12, 16:18)
	wt <- c(1:6)				# wt uninfected vs. infected
	persistent <- c(7:12)		# virus-resisting uninfected vs.  infected
	revertant <- c(13:18)		# revertant to virus-sensible uninf vs. inf

	##########################################################################################
	#
	# Now, we can proceed with the analysis of gene importance:
	# For Gallus data (outcome is categorical)
	#	we expect to identify the most important genes to "explain" the sample
	# (i.e., the genes that are more relevantly associated with the different 
	# cell types and infection status); they will likely represent 
	# constallations of related genes *whose association with the outcome is 
	# similar*.
	#
	##########################################################################################

    ### correlation
	#
	#tc.i.cor <- imp_cor(formula, tc)				 # 'y' must be numeric
	
	### relimp (relative importance)
	#
	#tc.i.relimp <- imp_relimp(formula, tc, tune=T)	 # wrong model type for classification
	
	### stepwise
	#
	#tc.i.step <- imp_stepwise(formula, tc)			 # ERROR: AIC is infinity for this model
	
	
	#### recursive feature elimination
	# 
	tc.i.rfe <- imp_rfe(formula, tc) 	# no results
	if (! is.null(tc.i.rfe)) {
        etc.i.rfe <- data.frame(gene.name=tc.i.rfe)
        etc.i.rfe <- enrich_genes(etc.i.rfe, bio.ann, 'gene.name', 'entrezgene_accession')
        write.table(etc.i.rfe, 'impgenes/imp_rfe_tc.txt', sep='\t', row.names=F)
    }
	
	
	#### mars
	#
	tc.i.mars <- imp_mars(formula, data=tc)		# y is a character variable
	if (! is.null(tc.i.mars)) {
        etc.i.mars <- cbind(gene.name=rownames(tc.i.mars), tc.i.mars)
        etc.i.mars <- enrich_genes(etc.i.mars, bio.ann, 
		                           'gene.name', 'entrezgene_accession')
        write.table(etc.i.mars, 'impgenes/imp_mars_tc.txt', sep='\t', row.names=F)
    }


	#### lvq (Learning Vector Quantization)
	#
	tc.i.lvq <- imp_lvq(formula, ftc) 
	if (! is.null(tc.i.lvq)) {
        etc.i.lvq <- cbind(gene.name=rownames(tc.i.lvq), tc.i.lvq)
        etc.i.lvq <- enrich_genes(etc.i.lvq, bio.ann, 
		                          'gene.name', 'entrezgene_accession')
        write.table(etc.i.lvq, 'impgenes/imp_lvq_tc.txt', sep='\t', row.names=F)
    }
	
	
	#### Boruta with RF
	#
	tc.i.boruta <- imp_boruta(formula, tc) 	# Unsupported type of dependent variable
	tc.i.boruta <- imp_boruta(formula, ftc)
	if (! is.null(tc.i.boruta)) {
        etc.i.boruta <- data.frame(gene.name=tc.i.boruta)
        etc.i.boruta <- enrich_genes(etc.i.boruta, bio.ann, 
		                             'gene.name', 'entrezgene_accession')
        write.table(etc.i.boruta, 'impgenes/imp_boruta_tc.txt', sep='\t', row.names=F)
    }
	selected_columns_tc <- cbind(sample=tc[ ,1], 
	                             tc[, colnames(tc) %in% tc.i.boruta])
	selected_columns_tc$sample <- as.factor(selected_columns_tc$sample)
    # generate plot using verbose=T
	selected.tc.i.boruta <- imp_boruta(formula, selected_columns_tc, 
	                                   out.png="impgenes/imp_selected_tc.png",
									   verbose=TRUE)
    
	# this will give us the cream of the cream (a refinement
	# of the genes initially selected)
	if (! is.null(selected.tc.i.boruta)) {
        sel.etc.i.boruta <- data.frame(gene.name=selected.tc.i.boruta)
        sel.etc.i.boruta <- enrich_genes(etc.i.boruta, bio.ann, 
		                             'gene.name', 'entrezgene_accession')
        write.table(sel.etc.i.boruta, 'impgenes/imp_boruta_selected_tc.txt', sep='\t', row.names=F)
    }
	
	
	# ### Boruta with ferns
	#
	tc.i.boruta.ferns <- imp_boruta(formula, tc, getImp='getImpFerns')
	
	if (! is.null(tc.i.boruta.ferns)) {
        etc.i.boruta.ferns <- data.frame(gene.name=tc.i.boruta.ferns)
        etc.i.boruta.ferns <- enrich_genes(etc.i.boruta.ferns, bio.ann, 
		                                   'gene.name', 'entrezgene_accession')
        write.table(etc.i.boruta.ferns, 'impgenes/imp_boruta.ferns_tc.txt', sep='\t', row.names=F)
    }
	
	# get a cleaner plot
	selected_columns_tc_f <- cbind(sample=tc[ ,1], tc[, 
	                               colnames(tc) %in% tc.i.boruta.ferns])
	selected_columns_tc_f$sample <- as.factor(selected_columns_tc_f$sample)
	# save plot using verbose=T
	selected.tc.i.boruta.ferns <- imp_boruta(formula, selected_columns_tc_f,
	               out.png= "impgenes/imp_selected_tc_ferns.png",
				   verbose=T)

	# this will give us the cream of the cream (a refinement
	# of the genes initially selected)
	if (! is.null(selected.tc.i.boruta.ferns)) {
        sel.etc.i.boruta.ferns <- data.frame(gene.name=selected.tc.i.boruta.ferns)
        sel.etc.i.boruta.ferns <- enrich_genes(etc.i.boruta.ferns, bio.ann, 
		                             'gene.name', 'entrezgene_accession')
        write.table(sel.etc.i.boruta.ferns, 'impgenes/imp_boruta_ferns_selected_tc.txt', sep='\t', row.names=F)
    }
	
	
	#### rpart
	#
	tc.i.rpart <- imp_rpart(formula , data=tc)
    
	if (! is.null(tc.i.rpart)) {
        etc.i.rpart <- data.frame(gene.name=names(tc.i.rpart), tc.i.rpart)
        etc.i.rpart <- enrich_genes(etc.i.rpart, bio.ann, 
		                            'gene.name', 'entrezgene_accession')
        write.table(etc.i.rpart, 'impgenes/imp_rpart_tc.txt', sep='\t', row.names=F)
    }
	
	#### random forest
	#
	# Note that when using 'rangerRF' we are already using a
	# permutation variable importance method, which is what
	# imp_vita (called later) is about.
	#tc.i.rf <- imp_random_forest(formula, tc, method='rf')	# Error in eval(predvars, data, env) : object 'TRNAD-GUC_1' not found
	tc.i.rfr <- imp_random_forest(dep.var='sample', tc, method='rangerRF.dep.var') 	
	
	if (exists(quote(tc.i.rf)) && ! is.null(tc.i.rf)) {
	    etc.i.rf <- tc.i.rf[rowSums(tc.i.rf) != 0,]
		etc.i.rf <- cbind(gene.name=rownames(etc.i.rf), etc.i.rf)
        etc.i.rf <- enrich_genes(etc.i.rf, bio.ann,
								  'gene.name', 'entrezgene_accession')
        write.table(etc.i.rf, 'impgenes/imp_rf_tc.txt', sep='\t', row.names=F)
    }
	if (exists(quote(tc.i.rfr)) && ! is.null(tc.i.rfr)) {
	    #etc.i.rfr <- tc.i.rfr[rowSums(tc.i.rfr) != 0,]
		#etc.i.rfr <- cbind(gene.name=rownames(etc.i.rfr), etc.i.rfr)
        # rangerRF
		etc.i.rfr <- data.frame(gene.name=names(tc.i.rfr), importance=tc.i.rfr)
        etc.i.rfr <- enrich_genes(etc.i.rfr, bio.ann,
								  'gene.name', 'entrezgene_accession')
        write.table(etc.i.rfr, 'impgenes/imp_rangerRF_tc.txt', sep='\t', row.names=F)
    }
	
	
	### regularized random forest
	#
	tc.i.rrf <- imp_rrf(formula , data=tc)
	
	if (exists(quote(tc.i.rrf)) && ! is.null(tc.i.rrf)) {
	    tc.i.rrf <- tc.i.rrf[rowSums(tc.i.rrf) != 0,]
        etc.i.rrf <- cbind(gene.name=rownames(tc.i.rrf), tc.i.rrf)
        etc.i.rrf <- enrich_genes(etc.i.rrf, bio.ann,
								  'gene.name', 'entrezgene_accession')
        write.table(etc.i.rrf, 'impgenes/imp_rrf_tc.txt', sep='\t', row.names=F)
    }

	
	#### DALEX
	# 
	# fails with random-forest
	# XXX with RF we get Illegal column names in formula interface 
	#tc.i.dalex <- imp_dalex(formula, tc, method='rangerRF')
	# this gives 'sum' is not meaningful for factors
	#tcr.i.dalex <- imp_dalex(data=tc, method='rangerRF.dep.var', dep.var='sample') 	
	# this gives 'sum' is not meaningful for factors
	#	(so, it seems DALEX does not support rangerRF for classification
	#
	# RF fails because of the gene names
	# we can try to "fix" them, do the calculation and "unfix" them
	tcx <- tc
	names(tcx) <- c('sample', paste("V", seq(2, length(names(tc))), sep=''))
	o.n <- data.frame(old=names(tc), new=names(tcx)
	# DALEX runs with RF, but as it is a lot slower than ranger
	# and DALEX makes many repeated calculations, it can take a
	# very long time: with 100CPUs this one took about one week
	# and a half to run!!!
	tcx.i.dalex <- imp_dalex(formula, data=tcx, method="RF")
	# we'd get a dataframe with first column being variable 
	# pseudo-names. We need to match them to o.n
	# and now recover old names
	
    if (exists(quote(tc.i.dalex)) && ! is.null(tc.i.dalex)) {
        etc.i.dalex <- tc.i.dalex ; colnames(etc.i.dalex) <- c('ensembl_gene_id', 'mean_dropout_loss')
        etc.i.dalex <- enrich_genes(etc.i.dalex, bio.ann)
        write.table(etc.i.dalex, 'impgenes/imp_dalex_tc.txt', sep='\t', row.names=F)
    }

	### vita
	#
	# XXX object 'TRNAD-GUC_1' not found
	#tc.i.vita <- imp_vita(formula , data=tc)
	tc.i.vita <- imp_vita(dep.var='sample', data=tc, method='rangerRF.dep.var')
	# SAVE ME (but no need)
	
	
	### LASSO
	#
	tc.i.lasso <- imp_lasso(formula , data=tc)

	if (exists(quote(tc.i.lasso)) && ! is.null(tc.i.lasso)) {
    	# convert list output to a table that we can annotate
		lasso.table <- data.frame(var=c('(Intercept)'))
		for (i in names(tc.i.lasso)) { 
    		if (length(names(tc.i.lasso[[i]])) == 0) names(tc.i.lasso[[i]]) <- '(Intercept)'

    		df <- data.frame(var=names(tc.i.lasso[[i]]), coef=tc.i.lasso[[i]])
    		lasso.table <- merge(x=lasso.table, y=df, by='var', all=T)
			colnames(lasso.table)[ dim(lasso.table)[2] ]  <- i

		}
		lasso.table[ is.na(lasso.table) ] <- 0
		# annotate
		elasso.table <- enrich_genes(lasso.table, bio.ann,
									 'var', 'entrezgene_accession')
        # save
		write.table(elasso.table, 'impgenes/imp_lasso_tc.txt', sep='\t', row.names=F)
    }
	
	
	#### xgboost
	#
	#WARNING: 'ntree_limit' is deprecated,use 'iteration_range' instead

    tc.i.xgboost <- imp_xgboost(formula, tc) #WARNING: use 'iteration_range' instead
	
    if (exists(quote(tc.i.xgboost)) && ! is.null(tc.i.xgboost)) {
        etc.i.xgboost <- data.frame(gene.name=rownames(tc.i.xgboost), tc.i.xgboost)
        etc.i.xgboost <- enrich_genes(etc.i.xgboost, bio.ann, 
								 'gene.name', 'entrezgene_accession')
        write.table(etc.i.xgboost, 'impgenes/imp_xgboost_tc.txt', sep='\t', row.names=F)
    }
	
	
	#### genetic algorithm
	#
	tc.i.ga <- imp_ga(formula, tc)
	
    if (exists(quote(tc.i.ga)) && ! is.null(tc.i.ga)) {
        etc.i.ga <- data.frame(gene.name=tc.i.ga)
        etc.i.ga <- enrich_genes(etc.i.ga, bio.ann, 
								 'gene.name', 'entrezgene_accession')
        write.table(etc.i.ga, 'impgenes/imp_ga_tc.txt', sep='\t', row.names=F)
    }
	

	#### simulated annealing
	#
	tc.i.sa <- imp_sa(formula, tc)
	
	if (exists(quote(tc.i.sa)) && ! is.null(tc.i.sa)) {
        etc.i.sa <- data.frame(gene.name=tc.i.sa)
        etc.i.sa <- enrich_genes(etc.i.sa, bio.ann, 
								 'gene.name', 'entrezgene_accession')
        write.table(etc.i.sa, 'impgenes/imp_sa_tc.txt', sep='\t', row.names=F)
    }
	

	##########################################################################################
	#
	# Now, we can proceed with the factor analysis:
	#	we expect to identify the most important groups of genes that determine the
	# outcome variable
	#
	##########################################################################################
	#
	# EFA
	#
	# with EFA, we calculate which genes behave coordinately across
	# the full dataset
	#
	# We can do EFA on categorical variables if we can provide a
	# correlation matrix.
    #
	# Unfortunately this does not work with counts because of excessive
	# categories (seems to take all variables as categorical, possibly
	# bcause categories are defined as a column).
	# 
	# To do EFA we need a correlation matrix. We can generate one for
	# categorical variables using polychoric():
	# library(polycor)
	# pc.cor <- polychoric(dat)$rho
	# or
	# pc <- hetcor(data, ML=TRUE)
	# and then run EFA
	#
	# we can determine the number of factors with fa.parallel() or
	# directly:
	# fap <- fa.parallel(data, cor='poly')   # parallel analysis for dichotomous data
	# vss(pc$correlations, n.obs=N, rotate="varimax")  # very simple structure
	#
	# Then, with the polychoric correlation matrix we can use fa():
	# faPC <- fa(r=pc$correlations, nfactors=NF, n.obs=N, rotate="varimax")
	# faPC$loadings
	#
	# it was once also possible to use fa.poly() from 'psych')
	# (but it is now deprecated in favor of fa() )
	# faPCdirect <- fa.poly(data, nfactors=NF, rotate="varimax")
	# faPCdirect$fa$loadings
	#
	# and do plots:
	# factor.plot(faPCdirect$fa, cut=0.5)
	# fa.diagram(faPCdirect)
	#
	# For factor scores, look at package ltm which has a factor.scores() 
	# function specifically for polytomous outcome data.
	# see https://github.com/drizopoulos/ltm/blob/master/Examples/Scoring.R
	# library(ltm)
	# fs <- factor.scores(faPC)
	# plot(fs, include.items = TRUE, main = "EFA.poly")
	#
	# We can try to remove the 'sample' column and make it into the
	# row names, but we cannot have duplicate row.names, and if we
	# leave the three replicas, they would be considered as different
	# outcomes, not collated as a single one.
	#
	# In conclusion: in the case of Coturnix, EFA worked best with
	# count data, but in the case of Gallus, we cannot rely on raw
	# (or normalized, multi-sample) count data for either EFA or CA.
	#
	# A possible solution might be to take the average of all three counts
	# for each category, move the categories to row-names and try again
	# but in that case, we might as well try to use log2FC instead, which
	# (in theory) better summarizes our data, although it would reflect
	# differences between categories and not the categories themselves.
	#
	# To start with, we will try to approximate "population counts" from
	# our sample counts by taking the average of the normalized counts
	# using aggregate(), and then try to do EFA and CA with them. We
	# can do this because we are using normalized counts.
	mtc <- aggregate( tc[ , 2:dim(tc)[2] ], by=list(tc$sample), FUN=mean )
	rownames(mtc) <- mtc$Group.1
	mtc$Group.1 <- NULL 

	# and now we do EFA on the 'estimated population counts'
	#mtc.fac.5 <- fac(mtc,nfactors=5, rotate="varimax")
	#as.png(fa.plot(mtc.fac.5),
	#	   'impgenes/mtc.fac.5.png')
    #as.png(fa.diagram(mtc.fac.5),
	#	   'impgenes/mtc.fac.5.diagram.png')
	#mtc.fac.6 <- fac(mtc,nfactors=6, rotate="varimax")
	#as.png(fa.plot(mtc.fac.6),
	#	   'impgenes/mtc.fac.6.png')
    #as.png(fa.diagram(mtc.fac.6),
	#	   'impgenes/mtc.fac.6.diagram.png')

	WE_HAVE_LOTS_OF_TIME_AND_CPU <- FALSE
	# NOTE we have 100 CPUs, 1TB memory, and still set this to FALSE!!!
	if ( WE_HAVE_LOTS_OF_TIME_AND_CPU ) {
    	# try to find out optimal number of hidden factors using fa.parallel()

    	mtc.cor <- cor(mtc, use="pairwise.complete.obs")
    	as.png(
    		corrplot(mtc.cor, order='hclust', tl.col="black"), #the plot is wrong
    		"efa/mtc_corplot.png", width = 800, height = 600)

    	# this may take very long: it will run EFA many times, increasing the number
    	# of factors to find an optimal number.
    	# Maybe we should consider doing it with a subsample? Would that be representative?
    	mtc.fa.paral <- fa.parallel(x=mtc.cor, fm="minres", fa="fa")
    	print(mtc.fa.paral)
	} else {
	    # 800 is on the verge of tolerability for just finding out
		# a prediction for the number of factors/components
		# doing it with different size subsets will give us an 
		# idea of a potential number and its dependence on size
		# (if we assume genes are at random, size should have little
		# to no impact or at least tend to stabilize at some point)
    	mtc100.fa.paral <- fa.parallel(x=mtc[ , 1:100], fm="minres", fa="fa")
    	print(mtc100.fa.paral)
    	mtc200.fa.paral <- fa.parallel(x=mtc[ , 1:200], fm="minres", fa="fa")
    	print(mtc200.fa.paral)
    	mtc400.fa.paral <- fa.parallel(x=mtc[ , 1:400], fm="minres", fa="fa")
    	print(mtc400.fa.paral)
	    	mtc800.fa.paral <- fa.parallel(x=mtc[ , 1:800], fm="minres", fa="fa")
    	print(mtc800.fa.paral)

	}
	# when run incrementally, there seems to be an agreement that
	# the optimal number of factors should be 2 according to fa.parallel()
	# using eigenvalue analysis with simulated data
	# and 6 by the elbow criterion as deduced from visual inspection of
	# the plots: 
	# we will try 6 and 5 (2 is too few).

	nFactors = 6
	st.seed(SEED)
	get_efa_factors(mtc, n.factors=6, annotation=bio.ann, 
	                force.recalculation=F, 
					out.prefix='efa/mtc', 
					verbose=TRUE)
	# this works

	nFactors = 5
	st.seed(SEED)
	get_efa_factors(data=mtc, n.factors=5, annotation=bio.ann,
				    force.recalculation=F, 
					out.prefix='efa/mtc', 
					verbose=TRUE)
	# this one raises an error:
	#	Error in e$vectors %*% diag(inv.sqrt.ev) : non-conformable arguments
	#
	# We can also try to do EFA with log2FC, this would tell us which
	# groups of genes are important to define each specific difference.
	#
	# Left for later.


	##########################################################################################
	#
    # Correspondence Analysis
	#
	##########################################################################################
	#
	# CA
	#
	# PCA maximizes the variance explained by each axis, while CA 
	# maximizes the correspondence between species and sample scores
	# (in our case the correspondence between cell types and genes).
	#
	# it should be as simple as issuing
	# data.ca <- CA(data)
	#
	# problem here is we need data to be a contingency table, i.e., the
	# sample names should be columns, and we cannot have duplicate names
	# so this might be better applied to log2FC data, although we can
	# make a contingency matrix by estimating the population counts from
	# the average of the three replicas of each sample.
	#
	### for CA using factominer
	#
	mtc <- aggregate( tc[ , 2:dim(tc)[2] ], by=list(tc$sample), FUN=mean )
	rownames(mtc) <- mtc$Group.1
	mtc$Group.1 <- NULL 
	#
	# CA will produce a plot, we will save it
	as.png(
		CA(mtc),
		'ca/mtc.ca.plot.png')
    mtc.ca <- CA(mtc),
	print(mtc.ca)
	# gives five dimensions (we have six levels)
	n.dim <- dim(mtc.ca$eig)[[1]]
	
	row <- get_ca_row(mtc.ca)
	as.png(corrplot(row$contrib, is.corr = FALSE), 
	       'impgenes/mtc_ca_contrib.png')
    
	# get eigenvalues / variances
	eig.val <- get_eigenvalue(mtc.ca)
	print(eig.val)
	as.png(
        fviz_screeplot(mtc.ca, addlabels = TRUE, ylim = c(0, 50)),
		'ca/mtc.ca.scree.png')
	as.png(
		fviz_screeplot(mtc.ca) + 
		geom_hline(yintercept = mean(eig.val), linetype = 2, color = "red"),
		'ca/mtc.ca.scree.ave.png')

	# biplot
	as.png(
    	fviz_ca_biplot(mtc.ca, repel = TRUE),
		'ca/mtc.ca.biplot.png')
	
	# graph of row variables
	# 	rows
	row <- get_ca_row(mtc.ca)
	row
    # 	row coordinates
	head(row$coord)
	#	quality
	head(row$cos2)
	#	contribution
	head(row$contrib)
	# 	plot
	as.png(
		fviz_ca_row(mtc.ca, repel = TRUE),
		'ca/mtc.ca.rows.png')
	as.png(
		fviz_ca_row(mtc.ca, col.row = "steelblue", shape.row = 15),
		'ca/mtc.ca.rows.tn.png')
	
	# quallity of representation of rows
	as.png(
		fviz_ca_row(mtc.ca, col.row = "cos2", 
					gradient.cols = c("#00AFBB", "#E7B800","#FC4E07"), 
					repel = TRUE),
		'ca/mtc.ca.rows.qc.png')
    # using transparency
	as.png(
		fviz_ca_row(mtc.ca, alpha.row = "cos2"),
		'ca/mtc.ca.rows.qc.tn.png')
	# visualize rows as big dots
	as.png(
		corrplot(row$cos2, is.corr = FALSE),
		'ca/mtc.ca.dot.rows.dot.png')
	# see as boxplot
	as.png(
		fviz_cos2(mtc.ca, choice='row', axes=1:2),
		'ca/mtc.ca.dot.rows.dot.tn.png')
	
	# contributions of rows to the dimensions
    as.png(
		corrplot(row$contrib, is.corr = FALSE),
		'ca/mtc.ca.dot.rows.contrib.png')
	# contribution of rows to dimension i
	# we shouldn't use top
	for (i in 1:n.dim) {
		as.png(
			fviz_contrib(mtc.ca, choice = "row", axes = i, top = 10),
			paste('ca/mtc.ca.dot.rows.contrib.PC', i, '.png', sep=''))
	}
	#fviz_contrib(mtc.ca, choice = "row", axes = 1, top = 10)	# DF1.I
	#fviz_contrib(mtc.ca, choice = "row", axes = 2, top = 10)	# DF1.PC DF1.PC.I
	#fviz_contrib(mtc.ca, choice = "row", axes = 3, top = 10)	# DF1
	#fviz_contrib(mtc.ca, choice = "row", axes = 4, top = 10)	# DF1.PC DF1.PC.I
	#fviz_contrib(mtc.ca, choice = "row", axes = 5, top = 10)	# DF1.P DF1.P.I
	#
	# PC2 and PC4 identify PC and PC.I cells
	# PC5 identifies P and P.I cells
	as.png(
	    fviz_contrib(mtc.ca, choice = "row", axes = c(2,4), top = 10),
	    'ca/mtc.ca.rows.contrib.dim24.png')
	# but PC2 and PC4 together also identify P cells in addition to PC and PC.I
	#
	# highlight most contributing rows
	as.png(
		fviz_ca_row(mtc.ca, col.row = "contrib", 
			gradient.cols = c("#00AFBB", "#E7B800","#FC4E07"), repel = TRUE),
		'ca/mtc.ca.rows.contrib.top.png')
	# with transparency
	as.png(
		fviz_ca_row(mtc.ca, alpha.row = "contrib", repel = TRUE),
		'ca/mtc.ca.rows.contrib.top.tn.png')
	
	# graphs of column variables
	#	NOTE:	we have too many genes to see anything clearly
	#  	As soon as we include columns (genes), we should consider
	# filtering plots by quality (adding 'select.col=list(cos2=0.8)'
	# or higher to display only most significant genes)
	#	NOTE: select.col=list(cos2=0.995|0.999) seems to allow some
	# seeing
	#
	col <- get_ca_col(mtc.ca)
	col
	head(col$coord)
	head(col$cos2)
	head(col$contrib)
	as.png(
		fviz_ca_col(mtc.ca),
		'ca/mtc.ca.cols.png')
	as.png(
		fviz_ca_col(mtc.ca, col.col = "cos2", 
			gradient.cols = c("#00AFBB", "#E7B800","#FC4E07"), repel = TRUE),
		'ca/mtc.ca.cols.qc.png')
	as.png(
		fviz_ca_col(mtc.ca, col.col = "cos2", 
		    select.col=list(cos2=0.995),
			gradient.cols = c("#00AFBB", "#E7B800","#FC4E07"), repel = TRUE),
		'ca/mtc.ca.cols.qc.dim12.0.995.png')
	as.png(
		fviz_ca_col(mtc.ca, col.col = "cos2", 
		    select.col=list(cos2=0.99), axes=c(2,4),
			gradient.cols = c("#00AFBB", "#E7B800","#FC4E07"), repel = TRUE),
		'ca/mtc.ca.cols.qc.dim24.0.99.png')
	as.png(
		fviz_cos2(mtc.ca, choice = "col", axes = 1:2),
		'ca/mtc.ca.cols.qc.tn.png')
	as.png(
		fviz_contrib(mtc.ca, choice = "col", axes = 1:2),
		'ca/mtc.ca.cols.contrib.png')

	# symmetric biplot
	as.png(
		fviz_ca_biplot(mtc.ca, repel = TRUE),
		'ca/mtc.ca.biplot.symm.png')
	# asymmetric biplot (using arrows)
	as.png(
		fviz_ca_biplot(mtc.ca, map = "rowprincipal", 
					   arrows = c(TRUE, TRUE), repel = TRUE),
		'ca/mtc.ca.biplot.asymm.arrows.png')
	# contribution biplot
	as.png(
		fviz_ca_biplot(mtc.ca, map = "colgreen", arrows = c(TRUE, FALSE), 
	            	   repel = TRUE),
		'ca/mtc.ca.biplot.contrib.png')
	
	# dimension description
	mtc.desc <- dimdesc(mtc.ca, axes = c(1,2))
	# dimension 1 on rows (cells)
	head(mtc.desc[[1]]$row, 4)
	# dimension 1 on columns (genes)
	head(mtc.desc[[1]]$col, 4)
	#
	# etc for the rest
	
	# having row%coord and col$coord we could
	#	gw <- list(names(row.coord)=c())
	#   min_dist_i_j <- +Inf
	#	for (i in names(col.coord) {	# for each gene
	#		for (j in names(row.coord) { # for each sample type
	#			dist_ij = calc_dist_ij(i, j, row.coord, col.coord, dims)
	#			if (dist_ij < min_dist_i_j) min_dist_i_j <- dist_i_j
	#		}
	#		#add gene i to gravity well of cell j
	#		gw[[j]] <- c(gw[[j]], i)
	#	}
	# 	calc_dist_ij(i, j, row.coord, col.coord, dims) {
	#		coords_i <- col.coord[ , dims]
	#		coords_j <- row.coord[ , dims]
	#		d <- sqrt(sum(sqrt(coords_i - coords_j)))
	#		return(d)
	#	}
	#
	
	
	#
	### for CA using vegan
	#
	#	with vegan we can use the same function as for CCA specifying
	# only the X matrix (if no Y matrix then CA is performed).
	mtc <- aggregate( tc[ , 2:dim(tc)[2] ], by=list(tc$sample), FUN=mean )
	rownames(mtc) <- mtc$Group.1
	mtc$Group.1 <- NULL 
	
	mtc.vca <- cca(X=mtc)
	# vegan is convenient because it facilitates access to results
	summary(mtc.vca)
	eigenvals(mtc.vca)
	summary(eigenvals(mtc.vca))
	# site (group) scores
	summary(mtc.vca)$sites 	# %>% round(3) %>% head()
	# species (genes) scores
	#	Note that the summary only returns the first 6 scores for each gene
	#   and site/sample, they can be extracted directly from mtc.vca
	#	but in our case that is not needed
	datos <- summary(mtc.vca)$species %>% round(3)
	datos_df <- as.data.frame(datos)
	df_sorted <- lapply(datos_df, function(x) head(sort(x, decreasing = TRUE), 30))
	
	
	# plot
	#	Note that as with CA the sheer amount of genes hides the
	#	BLUE dots for classes/samples
	as.png(
		plot(mtc.vca),				# same as for CA
		'ca/mtc.vca.plot.png')
	#
	# Gavin Simpson has written functions to organize the objects created 
	# by vegan in a consistent manner for graphing via ggplot2.
	# install.packages('remotes')
	# remotes::install_github("gavinsimpson/ggvegan")
	library(ggvegan)
	as.png(
		autoplot(mtc.vca),
		'ca/mtc.vca.aplot.png')
	
	# show sites names in the autoplot
	as.png(autoplot(mtc.vca, geom = "point") +
  		   geom_text(data = site_scores_df, 
		   aes(x = CA1, y = CA2, label = rownames(mtc)), 
		   vjust = -0.5, hjust = 0.7), 'impgenes/mtc.vca.aplot.rnam.png')

	# fortify() creates a dataframe with all the data needed to be
	# able to make ggplot2 plots
	mtc.vca.gg <- fortify(mtc.vca)
	as.png(
		ggplot(data = mtc.vca.gg,
    		   aes(x = CA1, y = CA2)) +
    		   geom_point(aes(shape = Score, colour = Score)) +
    		   facet_grid(facets = . ~ Score) +
    		   guides(shape = FALSE, colour = FALSE) +
    		   theme_bw(),
		'ca/mtc.vca.gplot.png')
	
	
	
	##########################################################################################
	#
	# Detrended Correspondence Analysis (DCA)
	#
	##########################################################################################
	#
	### DCA
	# Detrended Correspondence Analysis (DCA) is an extension to CA
	# that also uses one matrix of explanatory variables, and is 
	# thus also an indirect gradient analysis method. It may be
	# useful to correct some CA distortions, but it makes more difficult
	# to interpret the meaning of axes and eigenvalues as they no longer
	# represent the fraction of variance explained. For this reason
	# it is often not recommended.
	#
	mtc.dca <- decorana(mtc)
	as.png( {
		plot(mtc.dca)
		text(mtc.dca)
	}, "dca/mtc.dca.png")
	#
	# actually this brings all the groups together, which makes it
	# unattractive in our case.
	#


	##########################################################################################
	#
	# Canonical Correspondence Analysis (CCA)
	#
	##########################################################################################
	#
	# Canonical Correspondence Analysis (CCA) is another extension 
	# to CA, but uses a second matrix of explanatory variables and
	# so, it is a constrained or direct analysis method. In this case
	# we cannot use it unless we can come up with some metrics that
	# characterize each cell type/status.
	#
	# We need to think about this... but what would be the use of
	# establishing associations between groups of genes and "groups"
	# of cell types? Maybe we can make a table 
	# 		sample, replica, cell-type, infection,
	# but it would have to be numerical, and see if the method groups
	# samples according to cell type, replica or to infection. but as 
	# these are actually categorical, CCA should not work.
	#
	# On the other hand, as the number of variables increases, the
	# results of CA and CCA become similar, so its advantages may be 
	# limited in our case.

}
