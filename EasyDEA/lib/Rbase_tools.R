
#' getScriptPath()
#'
#' Get the name of the main script being run. This is FAR from ideal: it relies
#' on a series of heuristics: first it tries a method that works if the
#' script was run with 'source()', then a method that returns the first
#' --file= argument to Rscript, then returns the first -f argument to R
#' and finally tries to get the file name from the call stack.
#' 
#' @param	none
#' 
#' @return	full name of the script as given in the command line or
#'		source()
#'
#' @usage	getScriptPath()
#' 
#' @examples
#'		getScriptPath()
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
getScriptPath <- function() {
     
    # this works if we were called with 'source()'
    src.path <- getSrcFilename(function() {}, full.names=T)
    if (! is.null(src.path) && src.path != ""){
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
        #print(i); print(cmd.args[i])
        if (cmd.args[i] == '-f') return(cmd.args[i+1])
    }
    
    # if we arrive here, we didn't match anything, turn to last resort
    return (sys.frame(1)$ofile)
}


#' myname()
#'
#' myname() obtains the name of the calling function
#'
#' This function will return the name of the calling function.
#' The CALLING function !!!
#' THE CALLING FUNCTION !!!!!!
#' It is intended as a way to obtain a function's own name for error
#' reporting. 
#  We could do it directly, within a function using 
#  myname <- deparse(sys.call())
#  but it would be difficult to understand. Isolating this in a function
#  call makes it more readable.
#'
#' 
#' @return	the name of the calling function (one up in the call stack)
#'
#' @usage	me <- myname()
#'
#' @examples
#'	f <- function() { print(myname()) }
#'	f()
#'	##[1] "print(myname())"
#'	## myname() is being passed to print() as an argument and is thus
#'	## called by 'print', hence the output
#'	##
#'	f <- function() { me <- myname() ; print(me) }
#'	f()
#'	##[1] "f()"
#'	## myname() is called by f() to obtain its own name and then the name
#'	## is handled to print().
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
myname <- function() { 
    deparse(sys.calls()[[sys.nframe()-1]]) 
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
cat.err <- function(abort=FALSE, ...) { 
    #caller <- deparse( sys.calls()[[sys.nframe()-1]] )
    caller <- deparse( sys.call(-1) )
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
cat.warn <- function(...) { 
     #caller <- deparse( sys.calls()[[sys.nframe()-1]] )
     #me <- as.list(sys.call())[[1]]
     #parent <- as.list(sys.call(-1))[[1]]
     caller <- deparse( sys.call(-1) )
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

cat.info <- function(...) { 
     #caller <- deparse( sys.calls()[[sys.nframe()-1]] )
     #me <- as.list(sys.call())[[1]]
     #parent <- as.list(sys.call(-1))[[1]]
     caller <- deparse( sys.call(-1) )
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
use.package <- function(p, silent = F, VERBOSE = F) {
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
    if (! VERBOSE)
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



#' sourceDir
#'
#' sourceDir() sources all Q/R/S files in a directory
#'
#'	This may be convenient to source more than one file at once
#'	i.e. when getting dynamic constraints as a function,
#'	although it would likely make more sense if the file defining
#'	the function did include (source) itself all other required
#'	files explicitly.
#'	It may also be helpful to load all library files from a directory.
#'
#' @param	path	the path to the directory containing the files to source
#' @param	trace	Whether we want tracing output printed (default=TRUE)
#' @param	\dots	Any other parameters to pass on to source()
#'
#' @return	On end, any file ending in .(RrSsQq) present in the
#'		directory will have been sourced and the corresponding
#'		operations applied
#'
#' @usage	sourceDir('some/source/directory/')
#'
#' @examples
#'		sourceDir('./library/.')
#' 
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
sourceDir <- function(path, trace = TRUE, VERBOSE = F, ...) {
   for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
      if (VERBOSE) cat(nm,":", sep="")
      source(file.path(path, nm), ...)
      if (VERBOSE) cat(" LOADED\n")
   }
}


#' charAt
#'
#'	Extract character at a given position from a string (requires stringr)
#' 
#' @param	string	a character string
#' @param	pos	the position to extract
#' 
#' @return	the character at that position or '' if pos > nchar(string)
#'
#' @usage	charAt(string, pos)
#' 
#' @examples	charAt("Hi", 1)
#'		charAt("Hi", 2)
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
charAt <- function(string, pos) {
    str_sub(string, pos, pos)
}


#' short.title
#'
#'	A convenience function to print short titles using banner() from
#' package banneR if available or cat (requires charAt() which depends on
#' package stringr. This is intended to be used for short titles, but may
#' be used for longer titles at an aesthetic penalty.
#' 
#' @param	text	the text to print
#' 
#' @return	nothing
#'
#' @usage	shor.title(text)
#' 
#' @examples	short.title("Hello world")
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
short.title <- function(text='', level=1) {
    # we could modify this so level=1 produces capitals/==, level=2
    # produces small caps/--, level=3 small caps/ -, level=4 small
    # caps/ .
    # for now, let us ignore the level argument
    t.len <- nchar(text)
    if (exists("banner")) {
        if (t.len <= 10) {
            banner(text)
        } else {
            pieces <- ceiling(t.len / 10)
            for (i in (1:pieces)) {
                offset <- ((i - 1) * 10)
                banner(str_sub(text, offset + 1, offset + 10)) 
            }
        }
    } else {
        cat("\n============================================================\n")
        if (t.len < 30) {
            space <- (60 - (t.len * 2)) / 2
                    # output leading blanks to get centered text
            for (i in (1:space)) {
                cat(" ")
            }
        }
        for (i in (1:t.len)) {
            cat(charAt(text, i), " ", sep="")
            if ((i %% 30) == 0) cat('\n')
        }
        cat(  "\n============================================================\n")
    }
}


#'
#'
#'
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
#' @export
#
cat.line <- function(char, len=80) {
    for (i in (1:len)) {
        cat(char)
    }
    cat('\n')
}



#' title
#'
#'	A convenience function to print titles using banner() from
#' package banneR if available or cat (requires charAt() which depends on
#' package stringr. This is intended to be used for leveled titles.
#' 
#' @param	text	the text to print
#' @param	level	the title level (1:6, default = 1)
#' 
#' @return	nothing
#'
#' @usage	title(text)
#' 
#' @examples	title("Hello world")
#'		title("Hello", 2)
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
title <- function(text='', level=1) {
    line.len <- 80
    t.len <- nchar(text)

    if (level == 1) {
        text <- toupper(text)
        upchar <- '='
        dnchar <- '='
    } else if (level == 2) {
        upchar <- '-'
        dnchar <- '-'
    } else if (level == 3) {
        text <- tolower(text)
        upchar <- '.'
        dnchar <- '.'
    } else if (level == 4) {
        text <- tolower(text)
        upchar <- ' '
        dnchar <- '='
    } else if (level == 5) {
        text <- tolower(text)
        upchar <- ' '
        dnchar <- '-'
    } else {
        text <- tolower(text)
        upchar <- ' '
        dnchar <- '.'
    }

    cat('\n')
    cat.line(upchar, line.len)

    if (exists("banner")) {
        ch.per.line <- line.len / 8	# 8 is the character width in banner
        if (t.len <= ch.per.line) {
            banner(text)
        } else {
            pieces <- ceiling(t.len / ch.per.line)
            for (i in (1:pieces)) {
                offset <- ((i - 1) * ch.per.line)
                banner(str_sub(text, offset + 1, offset + ch.per.line)) 
            }
        }
    } else {
        if ( t.len < (line.len / 2) ) {
            space <- (line.len - (t.len * 2)) / 2
                    # output leading blanks to get centered text
            for (i in (1:space)) {
                cat(" ")
            }
        }
        for (i in (1:t.len)) {
            cat(charAt(text, i), " ", sep="")
            if ((i %% (line.len / 2)) == 0) cat('\n')
        }
    }
    
    cat.line(dnchar, line.len)
    cat('\n')
}

#' openLogFile()
#' 
#' openLogFile() opens a file to save all output.
#' The file will be removed if it exists if either 'remove=TRUE' or
#' 'overwrite=TRUE' is specified (any of them). Otherwise, behavior
#' depends on the value of 'append': when 'append=TRUE', the existing
#' log file is used and new log added at the end, otherwise, a new
#' version numbered log file is created so old logs are not lost. At
#' the end we return a handle to a writable file that can be used to
#' intersperse output should we be crazy enough to want to.
#'
#' @param	file	the name of the log file
#' @param	remove	whether an existing log file should be removed
#' @param	new.version	whether a new version should be used
#' @param	append	whether output should be added to an existing log
#'
#' @return	a reference to the log file that can be used to close it later
#'
#' @usage	log <- openLogFile('output.log')
#'
#' @examples
#'		log <- openLogFile('output.log')
#'		print(log)
#'		closeLogFile(log)
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#'
openLogFile <- function(file='zygR.log', overwrite = F, append = F, new.version=T) {
    if (overwrite) {
        if (file.exists(file)) file.remove(file)

        #logFile <- file(file, open="wt");
    } else if (new.version == TRUE) {
        i <- 0
        name <- sprintf("%s.%03d", file, i)
        # check the existence of increasing numbers starting at 0 of log file versions
        while (file.exists(name)) {
            i <- i + 1
            name <- sprintf("%s.%03d", file, i)
        }
        # when we find a free number, use it
        file=name
    }
    
    # open a writable connection and use it with sink()
    log.file <- file(file, open="wt");
    
    # split=T means that the output will go to the screen and to a file
    # since we use a new file, it does not make sense to append
    sink(file=log.file, type=c("output", "message"), split=T, append=append)

    # return the writable connection so it can be closed with close()
    return(log.file)
}


#' closeLogFile
#'
#' closeLogFile() - close a log file previously opened by openLogFile
#' 
#' @param	logFile	a handle to the log file returned by a previous call
#'			to openLogFile
#' 
#' @return	on end, the log file is guaranteed to be closed
#'
#' @usage	closeLogFile(log)
#' 
#' @examples
#'		log <- openLogFile('output.log')
#'		print(log)
#'		closeLogFile(log)
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
closeLogFile <- function(logFile)  {
    close(logFile)
    #while (sink.number() > 0) { sink() }
}



#' openlog
#'
#' A simple function to sink all output to a versioned log file
#' 
#' @param	logfile	the name of the log file
#' 
#' @return	nothing
#'
#' @usage	openlog(logfile)
#' 
#' @examples	log.file <- 'zygR.log'
#'		openlog(logile)
#'		### do something
#'		sink()	# close last log file opened
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
openlog <- function(logfile) {
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


# close all sinks
#' sink.titanic()
#'
#' Close all sinks.
#' 
#' @param	none
#' 
#' @return	none
#'
#' @usage	sink.titanic()
#' 
#' @examples
#' 		sink('zyg.1.log', split=T)
#' 		sink('zyg.2.log', split=T)
#'		sink.number()
#'		sink.titanic()
#'		sink.number()
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
sink.titanic <- function() {
    while (sink.number() >= 1) { sink() }
}


#' file.extension
#'
#' file.extension() returns the file extension suffix part of a file name
#' 
#' @param	path	the file name path
#' 
#' @return	the extension of the file (i.e. anything following after the
#'		last '.'
#'
#' @usage	ext <- file.extension('file/path.ext'
#' 
#' @examples
#'		ext <- file.extension('/some/path/to/a/file.ext')
#'		print(ext)
#'		## [1] "ext"
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
file.extension <- function (path) 
{
    parts <- strsplit(path, ".", fixed=T)[[1]]
    last <- parts[length(parts)]
    return(last)
}



#' loadTableFile
#'
#' loadTableFile() loads a file containing a table in any supported format
#'
#'	This is a general routine. It will first detect which formats are
#' supported in your R installation by trying to load various foreign format
#' access packages. If any of these packages is found, then all corresponding
#' file extensions are considered.
#'	Once the supported formats (and extensions) have been detected, the
#' name provided will be successively combined with each of the extensions and
#' if a file named 'file'+'.'+ ext is found, it will be read using the 
#' appropriate method. If no combination of 'file' + '.' + any of the supported
#' file extensions is found, then we will try to read the file using its
#' name as provided (considering its own extension if any). If no extension
#' match is found, then we will simply try to load it with 'read.table()'.
#'	A special provision is made for reading .R or .Rscript files: if the
#' successfully found file extension corresponds to an R script, then it will
#' be sourced, and the file will be searched for either, a table that has
#' the same name as the file (without extension), which will then be considered
#' the table to read, or a function with the same name as the file (without
#' extension), which will then be considered a setup function to create the
#' desired table, on demand. Note that this function, in turn, might as well 
#' read an external file to load the table internally, that would be hidden 
#' from us and would not make any difference.
#'	
#' XXX JR XXX NOTE: I need to run more exhaustive tests on this function
#
# This one could be simplified by ignoring the extension and smply
# trying all known formats in order
#	That'd be safer but "dirtier"
#' 
#' @param	file	the name (with ot without extension) of the file
#'			that contains the table we want to load
#'		stringsAsFactors	flag indicating whether strings should
#'			be read as factors or not (default FALSE)
#' @param	header	boolean flag indicating if the table contains a header
#'			(defaults to FALSE)
#' @param	check.names boolean flag to tell if table names should be 
#'			checked for R conformance (defaults to FALSE)
#'
#' 
#' @return	the table read
#'
#' @usage	table <- loadTableFile(file="myTable")
#' 
#' @examples
#'		## this should be enough
#'		table <- loadTableFile("myTable")	
#'		## but any of these should also work if you want to be specific
#'		table <- loadTableFile("myTable.tsv")
#'		table <- loadTableFile("myTable.csv")
#'		table <- loadTableFile("myTable.dat")
#'		table <- loadTableFile("myTable.txt")
#'		table <- loadTableFile("myTable.xlsx")
#'		table <- loadTableFile("myTable.arff")
#'		table <- loadTableFile("myTable.dbf")
#'		table <- loadTableFile("myTable.dta")
#'		table <- loadTableFile("myTable.rec")
#'		table <- loadTableFile("myTable.mtp")
#'		table <- loadTableFile("myTable.mat")
#'		table <- loadTableFile("myTable.sav")
#'		table <- loadTableFile("myTable.ssd")
#'		table <- loadTableFile("myTable.sas")
#'		table <- loadTableFile("myTable.sys")
#'		table <- loadTableFile("myTable.syd")
#'		table <- loadTableFile("myTable.xpt")
#'		table <- loadTableFile("myTable.dump")
#'		table <- loadTableFile("myTable.R")
#'		table <- loadTableFile("myTable.Rscript")
#'		table <- loadTableFile("myTable.Rdata")
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
loadTableFile <- function(file="", 
                          stringsAsFactors=FALSE, 
                          header=FALSE,
                          check.names=FALSE)
{
    if (file == "") {
        return(NULL)
    }
    
    options(stringsAsFactors=stringsAsFactors)
    
    # make a vector with all acceptable suffixes 
    suff <- c('tsv', 'csv', 'tab', 'dat', 'txt', 'R', 'Rscript', 'Rdata')
    if ((is.element('xlsx', installed.packages()[,1]))) {
        suff <- c(suff, 'xlsx')
    }
    if ((is.element('foreign', installed.packages()[,1]))) {
        suff <- c(suff, 'arff', 'dbf', 'dta', 'rec', 'mtp', 'mat', 
                  'sav', 'ssd', 'sas', 'sys', 'syd', 'xpt', 'dump')
    }
    fnames <- c(file)
    fnames <- c(fnames, paste(file, suff, sep='.'))
    avail <- file.exists(fnames)

    # if none is available, then ask for a file name
    if (all(avail == FALSE)) {
        file <- file.choose()
        fnames <- c(file)
        avail <- c(TRUE)
	suff <- c(file_ext(file))
    }
        
    table <- data.frame()
    for (i in 1:length(fnames)) {
        if (! avail[i]) {
            next
        }
        # get suffix
        ext <- file_ext(fnames[i])
        if (ext == 'tsv') {
            table <- read.table(fnames[i], sep="\t", 
                                header=header, check.names=check.names,
            			fill=T, strip.white=T, blank.lines.skip=T)
            if (! empty(table)) { return(table) }
        } else if (ext == 'tab') {
            table <- read.table(fnames[i], sep="\t", 
                                header=header, check.names=check.names,
            			fill=T, strip.white=T, blank.lines.skip=T)
            if (! empty(table)) { return(table) }
        } else if (ext == 'csv') {
            table <- read.table(fnames[i], sep=",", 
                                header=header, check.names=check.names,
            			fill=T, strip.white=T, blank.lines.skip=T)
            if (! empty(table)) { return(table) }
        } else if (ext == 'txt') {
            table <- read.table(fnames[i], 
                                header=header, check.names=check.names,
            			fill=T, strip.white=T, blank.lines.skip=T)
            if (! empty(table)) { return(table) }
        } else if (ext == 'dat') {
            table <- read.table(fnames[i], 
            			header=header, check.names=check.names,
            			fill=T, strip.white=T, blank.lines.skip=T)
            if (! empty(table)) { return(table) }
        } else if ((is.element('xlsx', installed.packages()[,1]))) {
	    if (ext == 'xlsx') {
                library(xlsx)
                # use xlsx2 to convert dates to POSIXct
                table <- read.xlsx2(fnames[i], 
                                    header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            }            
        } else if ((is.element('foreign', installed.packages()[,1]))) {
            library(foreign)
            if (ext == 'arff') {
                table <- read.arff(fnames[i], 
                                   header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if (ext == 'dbf') {
                table <- read.dbf(fnames[i], 
                                  as.is=stringsAsFactors, 
                                  header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if (ext == 'dta') {
                table <- read.dta(fnames[i], 
                         convert.factors=stringsAsFactors,
                         convert.dates=TRUE,
                         convert.underscore=TRUE, 
                         header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if (ext == 'rec') {
                table <- read.epiinfo(fnames[i], 
                                      get.broken.dates=TRUE,
                                      header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if (ext == 'mtp') {
                table <- read.mtp(fnames[i], 
                                  header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if (ext == 'mat') {
                table <- read.mtp(fnames[i], 
                                  header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if (ext == 'sav') {
                table <- read.spss(fnames[i], 
                                   use.value.abels=stringsAsFactors,
                                   to.data.frame=TRUE,
                                   header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if ((ext == 'ssd') || (ext == 'sas')) {
                table <- read.ssd(fnames[i], 
                                  header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if ((ext == 'sys') || (ext == 'syd')) {
                table <- read.systat(fnames[i], 
                                  to.data.frame=TRUE,
                                  header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if ((ext == 'xpt')) {
                table <- read.xport(fnames[i], 
                                    header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if (ext == 'dump') {
                #table <- data.restore(fnames[i]) 
                table <- read.S(fnames[i], 
                                header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            }
        } else if ((ext == 'R') || (ext == 'Rscript') ||
                   (ext == 'r') || (ext == 'rscript')) {
            # it is an R file that should contain a function named like the file
            source(fnames[i])
            # we could have used return(dget(data'.R') but then the file should
            # contain only an anonymous function, this is more versatile.
            # Check if an object named as data exists and if it is a function
            data <- file_path_sans_ext(basename(fname[i]))
            
            # if the content was a function
            if (exists(data) && is.function(get(data))) {
                # in principle, this is the same as return(get(data))
                # this way, we have together both ways of doing it
                return(eval(parse(text=data)))
            } 
            # if the content was an array or a data.frame
            else if (exists(data) && 
                       (is.array(get(data)) || is.data.frame(get(data)))) {
                # in case that instead of a function we got a data definition
                # (matrix is a kind of array)
                table <- data
            }else {
                cat('loadTableFile: ERROR - ', data, '.R{script} should contain R code\n',
                    '    as a function ', data, '(model, concentrations, fluxes, step)\n',
                    '    or a "', data, '" array or data.frame', sep='')
                return(NULL)
            }
            ### NOTE ### NEEDS TO BE TESTED
        } else if (((ext == "Rdata") || (ext == "rdata")) && exists(data)) {
            return(load(data))
        } else {
            # At this point we have tried all known suffixes. As a last resort, 
            # we will try to read it as a text file
            table <- read.table(fnames[i], header=header, check.names=check.names)
	    # and leave it as default if nothing else works
            if (! empty(table)) { return(table) }
        }
    }
    # if we reach here it is because no approach has worked
    cat("\n\nloadTableFile: File", data, "was not valid!\n\n\n")

    return(table)    
}




#'
#'
#'
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
continue.on.key <- function(interactive=F) {
    #usePackage('keypress')		# install 'keypress' if not available
    
    if (! interactive)
        return('')

    cat("\nPress any key to continue: ")
    key <- keypress()
    return(key)
}


#'
#'
#'
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
continue.on.enter <- function(prompt='Press ENTER to continue: ', interactive=F) {
    if (! interactive)
        return('')
    
    return(readline(prompt))
}


#'
#'
#'
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
more.columns <- function (data, columns=c(1:dim(data)[2]), lines=20, header='', prompt="More? ") {
    use.package('keypress')
#library(keypress)

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

#' get.name.of
#'
#' get a variable or any known R object name as a character string
#' 
#' @param v	the R object
#' 
#' @return	the R object's name as a character string
#'
#' @usage	get.name.of(v)
#' 
#' @examples
#'			my.var <- 1
#'			get.name.of(my.var)
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
get.name.of <- function(v) {
  deparse(substitute(v))
}


#' show.data.frame
#'
#' show the contents of a data.frame using gWidgets
#' 
#' @param	df	the data.frame to show
#' @param	window.name	the name for the display window (defaults to
#'			the name of the data frame
#' @param	visible	whether the window should be shown or only created
#'			to be displayed at a later time
#' 
#' @return	the widget for the newly created window
#'
#' @usage	show.data.frame(df, window.name, visible)
#' 
#' @examples	my.df <- data.frame(X=c(1,2,3,4), Y=c(5,6,7,8))
#'		show.data.frame(my.df)
#'		show.data.frame(my.df, "Have a look at the results")
#'		w <- show.data.frame(my.df, visible=F)
#'		visible(w) <- TRUE
#'		Sys.sleep(5)
#'		visible(w) <- FALSE
#'		Sys.sleep(5)
#'		visible(w) <- TRUE
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
show.data.frame <- function(df, 
                            window.name=get.name.of(df), 
                            visible=TRUE,
                            interactive=F) {
                            
    #usePackage('gWidgets2')
    
    window.name <- window.name
    window <- gwindow(title=window.name, visible=visible)
    tab <- gtable(df,
           container=window)
    if (interactive)
        visible(window) <- TRUE	# setting it to FALSE removes window
    return(window)
}




#'
#'
#'
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
                   overwrite=TRUE, VERBOSE = T) {
    if (is.null(file)) {
        tryCatch(print(PLOT))
    } else if ( overwrite || ! file.exists(file) ) {
        if (VERBOSE){
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
