get.options <- function()
{
	## 1) Default values for EasyDEA
    # we will also use them as argument 'default' in make_option()
    # (possibly modified after reading a defaul configuration file)
	
    PAIRED                  <- FALSE
    ALIGN 					<- FALSE
    BOTH 					<- FALSE
    USE.ONLINE.ANNOTATION 	<- TRUE
    USE.EDGER 				<- TRUE
    USE.DESEQ2 				<- TRUE
    reference 				<- "./ref/Gg_GRGc7b.fna"
    annotation 				<- "./ref/Gg_GRGc7b.gtf"
    release 				<- "GRCg7b"
    target.organism 		<- 'Gallus gallus'
    ens.version 			<- '112'
    mart.name 				<- 'ggallus_gene_ensembl'
    org.package 			<- "org.Gg.eg.db"
    ncbi.taxid              <- 9031         # human: 9606
	kegg.organism			<- "gga"
    n.genes 				<- 1000
    fastq.dir 				<- 'fastq_dir'
    alignment.dir 			<- 'STAR-aln'
	feature.count.dir		<- "feature_counts"
    rnaseq.out 				<- 'rnaseq.out'
    my.name 				<- 'J. R. Valverde'
    my.email 				<- '<jrvalverde@cnb.csic.es>'
    my.user 				<- 'sci'
    my.password 			<- "password"
    metadata 				<- './SampleInfo.tab'
    cpm.threshold 			<- 0.5
    significance.threshold 	<- 1
    design.column 			<- "sample"
    config.file 			<- NULL
    INTERACTIVE 			<- FALSE
    VERBOSE 				<- TRUE

	## 2) User-defined default values file defined by the user
	#
    # This is a default options file that the user may have with
    # base options that will be used to override the arbitrary
    # defaults we have defined here
	
    default.config.file <- '.EasyDEA.rc'
    
	# It should be formatted as option=value
    # with comment starting with a #
    # therefore it would be a valid R input file

    # Check if there is a user-defined global default options file that
    # will override global defaults within the user's HOME directory 
    if (file.exists(paste('~/', default.config.file, sep='')))
        source(paste('~', default.config.file, sep='/'), local=TRUE)

    # Check if there is a user-defined default options file within the
    # EasyDEA execution directory (the directory from which we call the
    # program)
    if (file.exists(paste('./', default.config.file, sep = "")))
        source(default.config.file, local=TRUE)
    
    # There remain still two chances for modification:
    # -- as a specific commmand line specified config file
    # -- as separate specific command line arguments
    #
    #
    # EXPLANATION:
    # ============
    #   This arrangement allows us to have 
    #       a user-defined default set of options (useful, e.g. in 
    # shared group accounts or as generic preferences) -- ~/.EasyDEA.rc
    #
    #       a project-defined set of options (e.g. tailored to a specific
    # organism or read type) -- ./.EasyDEA.rc
    #
    #       a set of files with analysis-specific options (e.g. pivoting
    # through different metadata columns/variables) -- and specify each
    # of them in the command line
    #
    #       a 'spur of the moment' modification to test options as command
    # line options
    #
    #
    # As an example, at SCS-CNB
    # -------------------------
    #   - we run analyses in a specific account, and often select similar
    # options in all our projects for consistence (e.g. output arrangement
    # and hierarchy, etc.), although we recognize we prefer self-contained
    # project-specific or analysis-specific configuration files
    #   - we work on projects related to different organisms, so we need
    # to indicate the organism (H. sapiens, C. japonica, G. gallus...), or
    # with paired or unpaired reads (and we need to state that too...), etc.
    # on a per-project basis; we keep each project in a separate directory,
    # so having a project-config helps
    #   - when we analyze an organism we often want to compare using various
    # metadata variables (e.g. "sample" will compare all samples against
    # each other, but as it considers only two each time, many reported genes 
    # will not be sample-specific, but comparison-specific, so we also compare 
    # by e.g. "wt" to compare 'wt' against all other groups pooled together 
    # to obtain only the genes that are specific to wt and not expressed in 
    # any of the others, etc.); ee store each of these as a separate 
    # configuration file in a 'cfg' subdirectory, saving the results in 
    # separate output directories so that we can keep a record of the 
    # parameters used for each  analysis (this is specilly noticeable in
    # running-times and results: the more groups to compare, the longer
    # the running time and the larger the number of DE genes detected in
    # each comparison, when pooling, only reference-specific genes are
    # reported, which are less and analyses run faster)
    #   - specially at the beginning of a project, we may need to test 
    # alternate values for an option and use the command line with a 
    # "generic name" output directory, together with interactive mode, to 
    # check values and stop the program, so we can find the appropriate 
    # value for a given option (e.g. cpm or significance threshold); for 
    # this, where we try options to see if they work, the command line  
    # is usually very convenient, specially to generate the config file
    # (which can simply be copied, once the run is seen to work as we want, 
    # from the saved config in the output directory).
    #
    #
    # _In any case_, the options used for each run are saved as a config
    # file in the output directory, just in case, but if you re-use the 
    # same output directory for different analyses (e.g. first with edgeR  
    # and then with DESeq2), only the parameters of the last analysis run
    # will be kept, hence our emphasis on config files for the record and
    # future reference, documentation and repeatability
    #
    #

###############################################################################

    # Now we define all the command line options:
    # we do not want to overwrite values not specified in the command line,
    # for that we need to be able to tell which values were not specified,
    # and this is easier if we leave all options with their default value
    # as NULL (by leaving it unspecified), although it is less convenient
    # for users
    # NOTE: we could use as default value the latest value of the option
    # (default plus config file modifications) and then simply reassign
    # them all
    cmd.line.defaults.null <- TRUE
    # NOTE: if we changed the name of the options to match those produced
    # by optparse, and used the config-modified defaults for all of them,
    # then we could simply return 'opt' and use 'attach(opt)'
    # WHY NOT? two reasons, we favor explicit assignments over attach, and
    # we added option-handling late in the development of the script, when
    # we realized we had too many configuration options and the current 
    # option variable names were already dispersed all over the program
    
	option_list <- list(
	
        #-h / --help is added by default

       make_option(c("-o", "--output-dir"), 
            type="character", 
            help=paste("the name of the folder to save the output [default", rnaseq.out, "]")),

       make_option(c("-a", "--alignment-dir"), 
            type="character", 
            help=paste("directory containing alignments SAM files [default", alignment.dir, "]")),
			
        make_option(c("-A", "--align"), 
            action="store_true", 
            help=paste("whether to align fastq reads to the reference [default", ALIGN, "]")),

        make_option(c("-B", "--both-ends"), 
            action="store",
            type="logical", 
            help=paste("(T or F) whether both ends of paired reads must match [default", BOTH, "]")),

       make_option(c("-C", "--feature-count-dir"), 
            type="character", 
            help=paste("directory containing the feature count data [default", feature.count.dir, "]")),

       make_option(c("-c", "--config-file"), 
            type="character", 
            help=paste("configuration file with option=value definitions [default", config.file, "]")),

       make_option(c("-D", "--DESeq2"), 
            action="store_true", 
            help=paste("whether to apply the DESeq2 analysis [default", USE.DESEQ2, "]")),

       make_option(c("-E", "--edgeR"), 
            action="store_true", 
            help=paste("whether to apply the edgeR analysis [default", USE.EDGER, "]")),

       make_option(c("-f", "--fastq-dir"), 
            type="character", 
            help=paste("directory containing the FASTq files [default", fastq.dir, "]")),

       make_option(c("-g", "--gtf"), 
            type="character", 
            help=paste("GTF/GFF file containing genome annoation [default", annotation, "]")),

       make_option(c("-i", "--interactive"), 
            action="store_true", 
            help=paste("whether to activate interactive mode [default", INTERACTIVE, "]")),

       make_option(c("-k", "--comparison-column"), 
            type="character", 
            help=paste("column in the metadata file to compare by [default", design.column, "]")),

       make_option(c("-K", "--KEGG-organism"), 
            type="character", 
            help=paste("the KEGG basename of the target organism (e.g. 'hsa', 'gga', 'cja') [default", kegg.organism, "]")),

        make_option(c("-L", "--onLine-annotation"), 
            action="store_true", 
            help=paste("whether to use online annotation [default", USE.ONLINE.ANNOTATION, "]")),

       make_option(c("-m", "--mart-name"), 
            type="character", 
            help=paste("the name of the MART to use for annotation[default", mart.name, "]")),

       make_option(c("-M", "--my-email"), 
            type="character", 
            help=paste("the email of the user [default", my.email, "]")),

       make_option(c("-n", "--n-genes"), 
            type="integer", 
            help=paste("the maximum number of top genes to save [default", n.genes, "]")),
       
       make_option(c("-N", "--my-name"), 
            type="character", 
            help=paste("the name of the user [default", my.name, "]")),

       make_option(c("-O", "--Org-package"), 
            type="character", 
            help=paste("the name of the Org package to use for annotation[default", org.package, "]")),

       make_option(c("-p", "--paired"),
            action="store", 
            type='logical',
            help=paste("(T or F) whether we are working with paired reads (unsets both-ends) [default", PAIRED, "]")),

       make_option(c("-P", "--my-sql-password"), 
            type="character", 
            help=paste("the local MySQL password [default", my.password, "]")),

       make_option(c("-r", "--reference"), 
            type="character", 
            help=paste("reference name (without extensions) [default", reference, "]")),

       make_option(c("-R", "--release"), 
            type="character", 
            help=paste("the release name of the reference [default", release, "]")),

       make_option(c("-s", "--significance-threshold"), 
            type="double", 
            help=paste("log2 threshold to consider a difference significant [default", significance.threshold, "]")),

       make_option(c("-t", "--target-organism"), 
            type="character", 
            help=paste("the target organism name [default", target.organism, "]")),

       make_option(c("-T", "--cpm-threshold"), 
            type="double", 
            help=paste("minimum number of CPM to keep a sample [default", cpm.threshold,"]")),

       make_option(c("-U", "--my-sql-username"), 
            type="character", 
            help=paste("the local MySQL Username [default", my.user, "]")),

       make_option(c("-v", "--verbose"), 
            action="store_true", 
            help=paste("whether to output more information [default", VERBOSE, "]")),

       make_option(c("-V", "--ENSEMBL-version"), 
            type="character", 
            help=paste("the Version of ENSEMBL to use for annotation[default", ens.version, "]")),

       make_option(c("-x", "--metadata"), 
            type="character", 
            help=paste("metadata file with eXtra eXperimental information [default", metadata, "]")),

       make_option(c("-X", "--ncbi-taxid"), 
            type="character", 
            help=paste("NCBI TaxID [default", ncbi.taxid, "]"))

    )
    # get command line options, if help option encountered print help and exit,
    # otherwise if options not found on command line then set defaults, 
    parser <- OptionParser(option_list=option_list, add_help_option=TRUE)
    opt <- parse_args(parser, convert_hyphens_to_underscores=TRUE)

    ## User-defined configuration file with arguments and options
	# Check if a config file was specified in the command line, and
    # source it before assigning the other command line variables
    if (! is.null(opt$config_file)) source(opt$config_file, local=TRUE)

    # Command-line options will take precedence over any file-defined options
    # The order of priority is, from lowest to highest:
    #
	#	1)Global options defined here (global, in this function) are opaque
    #		to the user
    #
	#	2)User-defined options specified at the global ($HOME) level 
    #		set by the user in an "invisible" (hidden) file
    #
	#	3)User-defined options specified at a directory level are more specific
    # 		and apply only when we run in that directory, we expect them
    #		to have been written more recently and be better remembered
    #		by the user, but they are still in a hidden file
    #
	#	4)User-defined options in a file specified in the command line will 
    #		typically be used to reduce command line length, the file will
    #		likely not be a hidden file, they will enjoy permanence, but 
    #		the contents will be visible only when the file is displayed, 
    #		otherwise they must be remembered
    #
	#	5) Command-line options are always evident, and are impermanent, so
    #		we assume that the user resorts to them for volatile changes and
    # 		tests when s/he wants to override other existing defaults
    #
	#	Hierarchy of options:
    #	1.Command-line-options > 2.Command-line-file >
	#	3.Project-defaults(current dir) > 4.User-defaults(HOME dir) >
	#	5.Global-defaults
    #		
    
    if (cmd.line.defaults.null == TRUE) {
        # Now, set any options given in the command line to override
        # default or config file values
        # we  know they were given in the command line if they are not NULL
        if (! is.null(opt$align)) ALIGN <- opt$align
        if (! is.null(opt$paired)) PAIRED <- opt$paired
        if (! is.null(opt$both_ends)) BOTH <- opt$both_ends
        if (! is.null(opt$online_annotation)) USE.ONLINE.ANNOTATION <- opt$online_annotation
        if (! is.null(opt$edgeR)) USE.EDGER <- opt$edgeR
        if (! is.null(opt$DESeq2)) USE.DESEQ2 <- opt$DESeq2
        if (! is.null(opt$reference)) reference <- opt$reference
        if (! is.null(opt$annotation)) annotation <- opt$annotation
        if (! is.null(opt$release)) release <- opt$release
        if (! is.null(opt$target_organism)) target.organism <- opt$target_organism
        if (! is.null(opt$ENSEMBL_version)) ens.version <- opt$ENSEMBL_version
        if (! is.null(opt$mart_name)) mart.name <- opt$mart_name
        if (! is.null(opt$Org_package)) org.package <- opt$Org_package
        if (! is.null(opt$ncbi_taxid)) ncbi.taxid <- opt$ncbi_taxid
	    if (! is.null(opt$KEGG_organism)) kegg.organism <- opt$KEGG_organism
        if (! is.null(opt$n_genes)) n.genes <- opt$n_genes
        if (! is.null(opt$fastq_dir)) fastq.dir <- opt$fastq_dir
        if (! is.null(opt$alignment_dir)) alignment.data <- opt$alignment_dir
        if (! is.null(opt$feature_count_dir)) feature.count.dir <- opt$feature_count_dir
        if (! is.null(opt$output_dir)) rnaseq.out <- opt$output_dir
        if (! is.null(opt$my_name)) my.name <- opt$my_name
        if (! is.null(opt$my_email)) my.email <- opt$my_email
        if (! is.null(opt$my_sql_username)) my.user <- opt$my_sql_username
        if (! is.null(opt$my_sql_password)) my.password <- opt$my_password
        if (! is.null(opt$metadata)) metadata <- opt$metadata
        if (! is.null(opt$cpm_threshold)) cpm.threshold <- opt$cpm_threshold
        if (! is.null(opt$comparison_column)) design.column <- opt$comparison_column
        if (! is.null(opt$config_file)) config.file <- opt$config_file
        if (! is.null(opt$verbose)) VERBOSE <- opt$verbose
        if (! is.null(opt$interactive)) INTERACTIVE <- opt$interactive
    } else {
        # command line options defaults were set to default+config-file values
        # this assume ALL options have a valid default
        ALIGN <- opt$align
        PAIRED <- opt$paired
        BOTH <- opt$both_ends
        USE.ONLINE.ANNOTATION <- opt$online_annotation
        USE.EDGER <- opt$edgeR
        USE.DESEQ2 <- opt$DESeq2
        reference <- opt$reference
        annotation <- opt$annotation
        release <- opt$release
        target.organism <- opt$target_organism
        ens.version <- opt$ENSEMBL_version
        mart.name <- opt$mart_name
        org.package <- opt$Org_package
        ncbi.taxid <- opt$ncbi_taxid
	    kegg.organism <- opt$KEGG_organism
        n.genes <- opt$n_genes
        fastq.dir <- opt$fastq_dir
        alignment.data <- opt$alignment_dir
        feature.count.dir <- opt$feature_count_dir
        rnaseq.out <- opt$output_dir
        my.name <- opt$my_name
        my.email <- opt$my_email
        my.user <- opt$my_sql_username
        my.password <- opt$my_password
        metadata <- opt$metadata
        cpm.threshold <- opt$cpm_threshold
        design.column <- opt$comparison_column
        config.file <- opt$config_file
        VERBOSE <- opt$verbose
        INTERACTIVE <- opt$interactive
    }
    # if we are working with unpaired reads, asking for both
    # ends to match is absurd
    if (PAIRED == FALSE) BOTH <- FALSE

    # return options list
    options <-list(
                   ALIGN=ALIGN,
                   PAIRED=PAIRED,
                   BOTH=BOTH,
                   USE.ONLINE.ANNOTATION=USE.ONLINE.ANNOTATION,
                   USE.EDGER=USE.EDGER,
                   USE.DESEQ2=USE.DESEQ2,
                   reference=reference,
				   annotation=annotation,
                   release=release,
                   target.organism=target.organism,
                   ens.version=ens.version,
                   mart.name=mart.name,
                   org.package=org.package,
                   ncbi.taxid=ncbi.taxid,
				   kegg.organism = kegg.organism,
                   n.genes=n.genes,
                   fastq.dir=fastq.dir,
                   alignment.dir=alignment.dir,
				   feature.count.dir=feature.count.dir,
                   rnaseq.out=rnaseq.out,
                   my.name=my.name,
                   my.user=my.user,
                   my.password=my.password,
                   metadata=metadata,
                   cpm.threshold=cpm.threshold,
                   design.column=design.column,
                   config.file=config.file,		# may be NULL if untouched
                   INTERACTIVE=INTERACTIVE,
                   VERBOSE=VERBOSE
                  )

    return(options)
}    

save.options <- function(options, file='')
{
    if (is.null(file))
        return(FALSE)
    
    if (! is.list(options))
        return(FALSE)
        
    if (file != '') {
        # we need to open a connection, otherwise only the last line 
        # will be saved
        out <- file(file, open="w+")
        if (! exists('out')) return(FALSE)
    } else {
        # we need to use '' so cat output to the console/sink
        out <- file
    }
    
    for (i in 1:length(options)) 
        cat(names(options)[i], '=', options[[i]], '\n', file=out)
    
    if (out != '') 
        close(out)
    
    return(TRUE)
}
