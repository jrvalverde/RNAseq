get.options <- function()
{
	## 1) Default values for EasyDEA
    # we could set them as argument 'default' in make_option()
	
    ALIGN 					<- FALSE
    BOTH 					<- TRUE
    USE.ONLINE.ANNOTATION 	<- TRUE
    USE.EDGER 				<- TRUE
    USE.DESEQ2 				<- TRUE
    reference 				<- "./ref/Gg_GRGc7b.fna"
    annotation 				<- "./ref/Gg_GRGc7b.gtf"
    release 				<- "GRCg6a"
    target.organism 		<- 'Gallus gallus'
    ens.version 			<- '106'
    mart.name 				<- 'ggallus_gene_ensembl'
    org.package 			<- "org.Gg.eg.db"
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
    # will override global defaults within the user HOME directory 
    if (file.exists(paste('~/', default.config.file, sep='')))
        source(paste('~/', default.config.file, sep=''))
		
    # Check if there is a user-defined default options file within the
    # EasyDEA installation directory
    if (file.exists(paste('../', default.config.file, sep = "")))
        source(paste('../', default.config.file, sep = ""))

###############################################################################

    # Now we define all the command line options
    # we do not want to overwrite values not specified in the command line,
    # for that we need to be able to tell which values were not specified,
    # and this is easier if we leave all options with their default value
    # as NULL (by leaving it unspecified)
    
	option_list <- list(
	
        #-h / --help is added by default

        make_option(c("-A", "--align"), 
            action="store_true", 
            help=paste("whether to align fastq reads to the reference [default", ALIGN, "]")),

        make_option(c("-B", "--both-ends"), 
            action="store_true", 
            help=paste("whether both ends of paired reads must match [default", BOTH, "]")),

        make_option(c("-L", "--onLine-annotation"), 
            action="store_true", 
            help=paste("whether to use online annotation [default", USE.ONLINE.ANNOTATION, "]")),

        make_option(c("-E", "--edgeR"), 
            action="store_true", 
            help=paste("whether to apply the edgeR analysis [default", USE.EDGER, "]")),

        make_option(c("-D", "--DESeq2"), 
            action="store_true", 
            help=paste("whether to apply the DESeq2 analysis [default", USE.DESEQ2, "]")),

       make_option(c("-r", "--reference"), 
            type="character", 
            help=paste("reference name (without extensions) [default", reference, "]")),

       make_option(c("-a", "--annotation"), 
            type="character", 
            help=paste("GTF/GFF file containing genome annoation [default", annotation, "]")),

       make_option(c("-R", "--release"), 
            type="character", 
            help=paste("the release name of the reference [default", release, "]")),

       make_option(c("-t", "--target-organism"), 
            type="character", 
            help=paste("the target organism name [default", target.organism, "]")),

       make_option(c("-V", "--ENSEMBL-version"), 
            type="character", 
            help=paste("the version of ENSEMBL to use for annotation[default", ens.version, "]")),

       make_option(c("-m", "--mart-name"), 
            type="character", 
            help=paste("the name of the MART to use for annotation[default", mart.name, "]")),

       make_option(c("-O", "--Org-package"), 
            type="character", 
            help=paste("the name of the Org package to use for annotation[default", org.package, "]")),

       make_option(c("-K", "--KEGG-organism"), 
            type="character", 
            help=paste("the KEGG basename of the target organism (e.g. 'hsa', 'gga', 'cja') [default", kegg.organism, "]")),

       make_option(c("-n", "--n-genes"), 
            type="integer", 
            help=paste("the maximum number of top genes to save [default", n.genes, "]")),

       make_option(c("-f", "--fastq-dir"), 
            type="character", 
            help=paste("directory containing the FASTq files [default", fastq.dir, "]")),

       make_option(c("-a", "--alignment-dir"), 
            type="character", 
            help=paste("directory containing alignments SAM files [default", alignment.dir, "]")),
			
       make_option(c("-c", "--feature-count-dir"), 
            type="character", 
            help=paste("directory containing the feature count data [default", feature.count.dir, "]")),

       make_option(c("-o", "--output-dir"), 
            type="character", 
            help=paste("the name of the folder to save the output [default", rnaseq.out, "]")),

       make_option(c("-N", "--my-name"), 
            type="character", 
            help=paste("the name of the user [default", my.name, "]")),

       make_option(c("-M", "--my-email"), 
            type="character", 
            help=paste("the email of the user [default", my.email, "]")),

       make_option(c("-U", "--my-sql-username"), 
            type="character", 
            help=paste("the local MySQL username [default", my.user, "]")),

       make_option(c("-P", "--my-sql-password"), 
            type="character", 
            help=paste("the local MySQL password [default", my.password, "]")),

       make_option(c("-x", "--metadata"), 
            type="character", 
            help=paste("metadata file with eXtra eXperiment information [default", metadata, "]")),

       make_option(c("-t", "--cpm-threshold"), 
            type="double", 
            help=paste("minimum number of CPM to keep a sample [default", cpm.threshold,"]")),

       make_option(c("-s", "--significance-threshold"), 
            type="double", 
            help=paste("log2 threshold to consider a difference significant [default", significance.threshold, "]")),

       make_option(c("-k", "--comparison-column"), 
            type="character", 
            help=paste("column in the metadata file to compare by [default", design.column, "]")),

        make_option(c("-v", "--verbose"), 
            action="store_true", 
            help=paste("whether to output more information [default", VERBOSE, "]")),

        make_option(c("-i", "--interactive"), 
            action="store_true", 
            help=paste("whether to activate interactive mode [default", INTERACTIVE, "]")),

       make_option(c("--config-file"), 
            type="character", 
            help=paste("configuration file with option=value definitions [default", config.file, "]"))
    )
    
    # get command line options, if help option encountered print help and exit,
    # otherwise if options not found on command line then set defaults, 
    parser <- OptionParser(option_list=option_list)
    opt <- parse_args(parser, convert_hyphens_to_underscores=TRUE)


    ## User-defined configutaion file with arguments and options
	# Check if a config file was specified in the command line, and
    # source it before assigning the other command line variables
	
    if (! is.null(opt$config_file)) source(opt$config_file)

    # Command-line options will take precedence over any file-defined options
    # The order of priority is, from lowest to highest:
    #
	#	1)Global options defined here (global, in this function) are opaque
    #		to the user
    #
	#	2)User-defined options specified at the global ($HOME) level have
    #		been set by the user but are "invisible" (hidden file)
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
	#	3.User-defaults(HOME dir) > 4.User-defaults(EasyDEA dir) >
	#	5.Global-defaults
    #		
    
    # Now, set any options given in the command line to override
    # default or config file values
    # we  know they were given in the command line if they are not NULL
    if (! is.null(opt$align)) ALIGN <- opt$align
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

    # return options list
    options <-list(
                   ALIGN=ALIGN,
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
