##==========================##
##   ANNOTATION FUNCTIONS   ##
##==========================##


#' build.offline.annotation
#'
#' Try different strategies to buid the refernece genome annotation 
#' database online.
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
build.offline.annotation <- function(
                                     folder,
                                     db.dir,
                                     target.organism,
                                     reference.gtf,
                                     release,
                                     ens.version,
                                     user, 
                                     pass,
                                     author,
                                     maintainer,
                                     license
                                    ) 
{

    # We'll try to build an EnsDb Package
    # so we can keep it locally for future use
    #EnsDbPackageDir <- paste(folder, 'EnsDb.Ggallus.v106', sep='/')
    #EnsDbPackageDir <- paste(folder, 'EnsDb.Cjaponica.v105', sep='/')
    EnsDbPackageDir <- paste(folder, db.dir, sep='/')

    # source script included with package 'ensembldb'
    scr <- system.file("scripts/generate-EnsDbs.R", package = "ensembldb")
    source(scr)

    # use our local GTF file if that didn't work
    if (dir.exists(EnsDbPackageDir) ) {
        sqlite <- paste(target.organism, ".", release, ".", ens.version, ".sqlite", sep='/')
        DBFile <- paste(folder, '/', sqlite)
    } else {
        # generate SQLite database in place from GTF. Produces e.g.
        #    ./gallus_gallus.GRCg6a.106.sqlite
        gtf.edb <- ensDbFromGtf(gtf=reference.gtf, 
	        organism=target.organism,
                genomeVersion=release,
                version=ens.version,
                destDir=folder)			# lacks entrezid
        #dir()
        
        # file name of the SQLite database created
        #DBFile <- paste(folder, 'gallus_gallus.GRCg6a.106.sqlite', sep='/')
        sqlite <- paste(target.organism, ".", release, ".", ens.version, ".sqlite", sep='/')
        DBFile <- paste(folder, '/', sqlite)

        # we'll select the SQLite database file generated from the GTF to
        # create an EnsDb R package 
        system(paste("mv", sqlite, DBFile))

        makeEnsembldbPackage(DBFile, version=ens.version, 
    		             maintainer=maintainer, 
    		             author=author,
                             destDir=folder, license=license)
    }

    # retry with GFF if that didn't work
    if ( ! dir.exists(EnsDbPackageDir) ) {
        # As an alternative, we may use the GFF3 file (which should contain
        # the same information
        gff.edb <- ensDbFromGff(gff=reference.gff, 
                organism=target.organism,
                genomeVersion=release,
                version=ens.version,
                destDir=folder)
        # but it seems to lack entrezid and transcript_id

        #dir()
        # file name of the SQLite database created
        sqlite <- paste(target.organism, ".", release, ".", ens.version, ".sqlite", sep='/')
        DBFile <- paste(folder, '/', sqlite)

        # we'll select the SQLite database file generated from the GFF to
        # create an EnsDb R package 
        system(paste("mv", sqlite, DBFile))

        makeEnsembldbPackage(DBFile, version=ens.version, 
    		             maintainer=maintainer, 
    		             author=author,
                             destDir=folder, license=license)
    } else {
        sqlite <- paste(target.organism, ".", release, ".", ens.version, ".sqlite", sep='/')
        DBFile <- paste(folder, '/', sqlite)
    }

    # we still have a third option, build it from a MySQL file
    if ( ! file.exists(EnsDbPackageDir) ) {
        # This shouln't be needed because we have already done it above
        # using the GTF/GFF3 file.
        # This is an alternate way to do generate the EnsDb package from
        # MySQL data downloaded from ENSEMBL
        local.mysql.db <- paste("mysql/", 
        	tolower(target.organism), 
                "_core_", ens.version, "_", release, sep='')
        #local.mysql.db <- "mysql/coturnix_japonica_core_104_2"
        #local.mysql.db <- "mysql/gallus_gallus_core_105_6"
        createEnsDbForSpecies(ens_version = ens.version,
                user = user, 
                pass=pass,
                host = "localhost", 
                local_tmp=local.mysql.db, 
                species=target.organism,
                dropDb=FALSE)
        # This should create a .sqlite file like ensDbFromGtf above,
        # named, e.g.
        #	Gallus_gallus.GRCg6a.106.sqlite
        #       ^
        # with similar contents to the one produced from ensDBFromGtf/Gff;
        # if none of these does work, then tweak the process by hand in
        # source("scripts/generate-EnsDBs.R")
        sqlite <- paste(target.organism, ".", release, ".", ens.version, ".sqlite", sep='/')
        DBFile <- paste(folder, '/', sqlite)
        system(paste("mv", sqlite, DBFile))
        makeEnsembldbPackage(DBFile, version=ens.version, 
    		             maintainer=maintainer, 
    		             author=author,
                             destDir=folder, license=license)
                             
    } 
    return(DBFile)
}



#' get.biomart.ensembl.annotation
#'
#' obtain ENSEMBL-related annotation using BiomaRt either from a local
#' cache file previously saved or directly from the network
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
get.biomart.ensembl.annotation <- function(mart.db, folder) {

    if ( ! file.exists( paste(folder, 'biomart.ensembl.tab', sep='/')) ) {
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
    return(bm.ensembl.annot)
}


#' get.biomart.entrez.annotation
#'
#' obtain ENTREZ/NCBI-related annotation from BiomaRt, using either local
#' cached data or directly from the network
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
get.biomart.entrez.annotation <- function(mart.db, folder) {
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
    return(bm.entrez.annot)
}



#' get.biomart.go.annotation
#'
#' obtain GO annotation from BiomaRt
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
get.biomart.go.annotation <- function(mart.db, folder) {
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
    return(bm.go.annot)
}



#' get.biomart.goslim.annotation
#'
#' obtain GOslim annotatin from BiomaRt
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
get.biomart.goslim.annotation <- function(mart.db, folder) {
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
}


#' get.biomart.family.annotation
#'
#' obtain protein family (PFAM, TIGRFam, etc...) annotation from BiomaRt
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
get.biomart.family.annotation <- function(mart.db, folder) {
    if ( ! file.exists(paste(folder, '/biomart.fam.tab', sep='')) ) {
        bm.fam.annot <- getBM(attributes=c(
                               "ensembl_gene_id",
                               "pfam",
                               "pirsf",
                               "prints",
                               "tigrfam"
                               ), 
                           mart=mart.db)
        write.table(bm.fam.annot, paste(folder, '/biomart.fam.tab', sep=''), 
	        row.names=T, col.names=T, sep='\t')
    } else {
        bm.fam.annot <- read.table(paste(folder, '/biomart.fam.tab', sep=''), 
	        header=T, sep='\t')
    }
    return(bm.fam.annot)
}


#' get.biomart.prosite.annotation
#'
#' obtain PROSITE annotation from BiomaRt
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
get.biomart.prosite.annotation <- function(mart.db, folder) {
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
    return(bm.prosite.annot)
}



#' get.biomart.superfamily.annotation
#'
#' obtain SuperFam annotation from BiomaRt
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
get.biomart.superfamily.annotation <- function(mart.db, folder) {
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
    return(bm.sfam.annot)
}


#' get.biomart.extra.annotation
#'
#' obtain additional miscellaneous annotation from BiomaRt
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
get.biomart.extra.annotation <- function(mart.db, folder) {
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
    return(bm.extra.annot)
}


#' biomart.merge.annotations
#'
#' merge together various annotation datasets into a single, consolidated one
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
biomart.merge.annotations <- function(ann, by="ensembl_gene_id") {
    if (! file.exists(paste(folder, '/biomaRt.annotation.txt', sep=''))) {
        # ann is a list of annotations to merge
        if (length(ann) == 1) return(ann)		# nothing to merge

        # we have at least two
        bm.annot <- merge(ann[[1]], ann[[2]], by=by)

        if (length(ann == 2)) return(bm.annot)

        # we have more annotations to merge
        for (i in 3:length(ann))
            bm.annot <- merge(bm.annot, ann[[i]], by=by)
        write.table(bm.annot, file=paste(folder, '/biomaRt.annotation.txt', sep=''), 
                    sep='\t', row.names=T, col.names=T)
    } else {
        bm.annot <- read.table(paste(folder, '/biomaRt.annotation.txt', sep=''), 
	        header=T, sep='\t')
    }
    return(bm.annot)
}

get.biomart.annotation <- function(annotation.dir) {
    # get all the annotation
	
	## Set up connection to ensembl database
	marts <- listMarts()
	if (VERBOSE) head(marts)
	bm.ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
	bm.ens.datasets <- listDatasets(bm.ensembl)
	if (VERBOSE) head(bm.ens.datasets)
		
	## Get BiomaRt dataset name using RegExp (mart.name is global)
	target.ds <- subset(bm.ens.datasets, grepl(mart.name, dataset))
	biomart.ds.name <- target.ds$dataset		
	
	## Load dataset
	cat("\nLoading dataset", biomart.ds.name, "\n")
	if ( ! nrow(target.ds) > 0) { stop("BiomaRt Dataset Not Found!") }
	mart.db <- useDataset(biomart.ds.name, mart = bm.ensembl)

    if (VERBOSE) {
        #listDatasets(ensembl) %>%  filter(str_detect(description, release))
        cat("\n   > BiomaRt Attributes")
		listAttributes(mart.db)[,1:2] # %>% head(20)
        cat("\n   > BiomaRt Filters")
        listFilters(mart.db) # %>% head(20)
    }
    
    # we cannot get all the annotation at once because it times out
    #full.annot <- getBM(attributes=
    #                       c("ensembl_gene_id", "ensembl_transcript_id", 
    #			   				"start_position", "end_position", 
    #                          "chromosome_name", "gene_biotype", 
    #                          "description", 
    #                          "entrezgene_id", "entrezgene_accession", "entrezgene_description", 
    #                          "go_id", "name_1006", "definition_1006", "go_linkage_type", "namespace_1003", 
    #                          "goslim_goa_accession", "goslim_goa_description", 
    #                          "pdb", 
    #                          "reactome", "uniprotswissprot"), 
    #                       mart=mart.db)
    #
    # so we will retrieve the data in pieces, including ensembl_gene_id in
    # each piece so we can use it as key for merging the annotation and
    # saving it in a file to avoid future downloads
    # each of these calls will cache locally the annotation ro speed up
    # subsequent accesses
    bm.ensembl.annot <- get.biomart.ensembl.annotation(mart.db, annotation.dir)

    bm.entrez.annot <- get.biomart.entrez.annotation(mart.db, annotation.dir)

    bm.go.annot <- get.biomart.go.annotation(mart.db, annotation.dir)

    bm.goslim.annot <- get.biomart.goslim.annotation(mart.db, annotation.dir)

    bm.fam.annot <- get.biomart.family.annotation(mart.db, annotation.dir)

    bm.prosite.annot <- get.biomart.prosite.annotation(mart.db, annotation.dir)

    bm.sfam.annot <- get.biomart.superfamily.annotation(mart.db, annotation.dir)

    bm.extra.annot <- get.biomart.extra.annotation(mart.db, annotation.dir)

    # Now that we have all the pieces, merge them all together
    # into a single annotation variable
    #	THIS TAKES TOO LONG AND TOO MUCH MEMORY, COMMENTED FOR NOW
    #bm.annot <- biomart.merge.annotations(list(
    #                  bm.ensembl.annot,
    #                  bm.entrez.annot,
    #                  bm.go.annot,
    #                  bm.goslim.annot,
    #                  bm.fam.annot,
    #                  bm.prosite.annot,
    #                  bm.sfam.annot,
    #                  bm.extra.annot
    #                ),
    #                by="ensembl_gene_id", folder)
    #
    #write.table(bm.annot, file=paste(folder, '/biomaRt.annotation.txt', sep=''), 
    #        sep='\t', row.names=T, col.names=T)


    # Now that we have the annotation we can select unique entries
    # for our dataset

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

    # or we could deduplicate everything first and match afterwards
    # this has the advantage that we keep all the annotation at hand and
    # can reuse it for any gene dataset instead of annotating each
    # specific dataset separately
    bm.ensembl.annot.1 <- bm.ensembl.annot[ ! duplicated(bm.ensembl.annot$ensembl_gene_id), ]
    bm.entrez.annot.1 <- bm.entrez.annot[ ! duplicated(bm.entrez.annot$ensembl_gene_id), ]
    bm.go.annot.1 <- bm.go.annot[ ! duplicated(bm.go.annot$ensembl_gene_id), ]
    bm.goslim.annot.1 <- bm.goslim.annot[ ! duplicated(bm.goslim.annot$ensembl_gene_id), ]
    bm.fam.annot.1 <- bm.fam.annot[ ! duplicated(bm.fam.annot$ensembl_gene_id), ]
    bm.sfam.annot.1 <- bm.sfam.annot[ ! duplicated(bm.sfam.annot$ensembl_gene_id), ]
    bm.prosite.annot.1 <- bm.prosite.annot[ ! duplicated(bm.prosite.annot$ensembl_gene_id), ]
    bm.extra.annot.1 <- bm.extra.annot[ ! duplicated(bm.extra.annot$ensembl_gene_id), ]

    # this should now be manageable (many entries should have been removed)
    bm.annot.1 <- merge(bm.ensembl.annot.1, bm.entrez.annot.1, by="ensembl_gene_id")
    bm.annot.1 <- merge(bm.annot.1, bm.go.annot.1, by = "ensembl_gene_id")
    bm.annot.1 <- merge(bm.annot.1, bm.goslim.annot.1, by = "ensembl_gene_id")
    bm.annot.1 <- merge(bm.annot.1, bm.fam.annot.1, by = "ensembl_gene_id")
    bm.annot.1 <- merge(bm.annot.1, bm.sfam.annot.1, by = "ensembl_gene_id")
    bm.annot.1 <- merge(bm.annot.1, bm.prosite.annot.1, by = "ensembl_gene_id")
    bm.annot.1 <- merge(bm.annot.1, bm.extra.annot.1,  by = "ensembl_gene_id")

    write.table(bm.annot.1, file = paste(annotation.dir, 'biomaRt.annotation.1st.txt', sep='/'), 
				sep='\t', row.names=T, col.names=T)
    # we save it to avoid repeating this in the future

    # or even do it all at once? if we had enough power for building bm.annot:
    ##bm.annot.1 <- bm.annot[ ! duplicated(bm.annot$ensembl_gene_id), ]
    ##write.table(bm.annot.1, 
    ##            file=paste(folder, '/biomart.annotation.1st.txt', sep=''), 
    ##            sep='\t', row.names=T, col.names=T)
    return(bm.annot.1)
}


prepare.Ensembl.Db <- function(	target.organism,
								annotation.dir,
								use.online.annotation = TRUE,
								annotation,
								ensembl.version){

	    ############################################
	    ## USE ONLINE ENSEMBL ANNOTATION DATABASE ##
	    ############################################
	    
	# SELECT ENSEMBL ANNOTATION DATABASE AVAILABLE ON-LINE
    if (use.online.annotation == TRUE) {
		
		cat("\n							\n")
		cat("\tFETCHING ENSEMBL DATABASE\n")
	    cat("\t=========================\n")

		# Create output directory
		net.ens.out <- paste(annotation.dir, 'net.ens.Db', sep = '/')
	    dir.create(net.ens.out, showWarnings=FALSE)

	    #if ( ! file.exists(paste(net.ens.out, 'net.ens.db.rds', sep='/')) ) {

		# Get ENSEMBL annotations for our target
            ah <- AnnotationHub()
            qr <- query(ah, c("EnsDb", target.organism))

		# Choose the most recent (last) one: AH98040 Ensembl Db version 105
		last.ref <- tail(names(qr), n = 1)
		net.ens.db <- qr[[last.ref]]

		# Save online ensemble annotation
		saveRDS(net.ens.db, file=paste(net.ens.out, 'net.ens.db.rds', sep='/'))
		save(net.ens.db, file=paste(net.ens.out, 'net.ens.db.RData', sep='/'))


        ## ERROR LOADING:
        ## Extarnal pointer is not valid. ¿Can online DBs be saved as objects? 
        #} else {
	    #
        #    net.ens.db <- readRDS(file=paste(net.ens.out, 'net.ens.db.rds', sep='/'))
        #}

	    # Check it
	    print(net.ens.db)
	    columns(net.ens.db)
	    #head(keys(net.ens.db, 'GENEID'))

	    return(net.ens.db)
    }


	    #######################################
	    ## BUILD ENSEMBL ANNOTATION DATABASE ##
	    #######################################

    ### BUILD DATABASE FROM GTF/GFF,

    if (use.online.annotation == FALSE) {
		
		cat("\n                            \n")
		cat("\tBUILDING ENSEMBL DATABASE\n")
	    cat("\t=========================\n")

	    # Load required packages
        source(system.file("scripts/generate-EnsDBs.R", package = "ensembldb"))

        # Use our local GTF file to generate ENSEMBL Database
	    gtf <- paste(sub('\\.g(f|t)f$', '', annotation), 'gtf', sep='.')
	    gff <- paste(sub('\\.g(f|t)f$', '', annotation), 'gff', sep='.')
	    ref.base <- sub('\\.g(f|t)f$', '', basename(annotation))

	    # Create output directory	
        ens.pkg.out <- paste(annotation.dir, 'ens.pkg.Db', sep='/')
	    dir.create(ens.pkg.out, showWarnings = FALSE)

	    ##############
	    ##	ADRIAN  ##
	    ##############

	    # WARNING MESSAGES:
	    # I'm missing column(s): 'gene_name','entrezid'. The corresponding database column(s) will be empty!
	    # No column 'exon_id' present, created artificial exon IDs by concatenating the transcript ID and the exon number.
	    # Could not determine length for all seqnames.  Unable to retrieve sequence lengths from Ensembl.

	    sqlite <- paste(basename(gtf), 'sqlite', sep = '.')
	    #sqlite <- paste(target.organism, basename(gtf), ensembl.version, 'sqlite', sep='.')
	    ens.pkg.db <- paste(ens.pkg.out, sqlite, sep = '/') 
		
        if ( ! file.exists(ens.pkg.db) ) {

		    # Generate SQLite database in place from GTF
		    #gtf.edb <- ensDbFromGtf(gtf=gtf, 
	        #	    organism=target.organism,
            #       genomeVersion=basename(gtf),
            #       version=ensembl.version,
            #       destDir=ens.pkg.out)			# lacks entrezid

            ensDbFromGtf(	gtf = gtf, 
        				    organism= target.organism, 
        				    genomeVersion = basename(gtf),
        				    version = ensembl.version,
        				    outfile = ens.pkg.db)

		    # Move the database to the output directory
            # system(paste("mv", sqlite, ens.pkg.db))	

		    # Generate Ensembl DB Package
            makeEnsembldbPackage(	ensdb = ens.pkg.db,
        						    version = ensembl.version, 
    		             		    maintainer="J. R. Valverde <jrvalverde@cnb.csic.es>", 
    		             		    author="J. R. Valverde",
                         		    destDir=ens.pkg.out,
                         		    license="Artistic-2.0")
	    }

	    # If building DB from GTF file fails, we may use the GFF3 file (which should contain
	    # the same information

	    ##############
	    ##  ADRIAN  ##
	    ##############

	    # DB FROM GFF DOES NOT WORK. ERROR MESSAGE:
	    # Required columns/fields gene_id;exon_id;biotype not present in the GFF file!

	    sqlite <- paste(basename(gff), 'sqlite', sep = '.')
	    #sqlite <- paste(target.organism, basename(gff), ensembl.version, 'sqlite', sep='.')
	    ens.pkg.db <- paste(ens.pkg.out, sqlite, sep = '/') 
		

        if ( ! file.exists(ens.pkg.db) ) {

		    #gff.edb <- ensDbFromGff(gff=gff, 
		    #	organism=target.organism,
		    #	genomeVersion=basename(gff),
		    #	version=ensembl.version,
		    #	destDir=ens.pkg.out)			# lacks entrezid

		    ensDbFromGff(	gff = gff, 
        				    organism= target.organism, 
        				    genomeVersion = basename(gff),
        				    version = ensembl.version,
        				    outfile = ens.pkg.db)

	    # Move the database to the output directory
            # system(paste("mv", sqlite, ens.pkg.db))	

	    # Generate Ensembl DB Package
            makeEnsembldbPackage(	ensdb = ens.pkg.db,
        		            	version = ensembl.version, 
    		             	    maintainer="J. R. Valverde <jrvalverde@cnb.csic.es>", 
    		             	    author="J. R. Valverde",
                         	    destDir=ens.pkg.out,
                         	    license="Artistic-2.0")


            # Move the database to the output directory
            system(paste("mv", sqlite, ens.pkg.db))				

            makeEnsembldbPackage(ens.pkg.db, version=ensembl.version, 
							    maintainer="J. R. Valverde <jrvalverde@cnb.csic.es>", 
							    author="J. R. Valverde",
							    destDir=ens.pkg.out, license="Artistic-2.0")
        }
		return(ens.pkg.db)
    }
}



#' get_saved_annotation
#'
#
# e.g.
#   if (exists(quote(bm.annot.1)))
#       bio.ann <- bm.annot.1
#   else
#       bio.ann <- get_saved_annotation(ensembl.db=NULL, 
#                                       biomart.db='./biomaRt.annotation.1st.txt')$biomart
#
#   if (is.null(bio.ann)) cat.err("Could not read BIOMART database\n")
#
get_saved_annotation <- function(ensembl.db=NULL, biomart.db=NULL, verbose=F) {
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

	###################################################
	## EXTRACT GENE ANNOTATION FROM ENSEMBL DATABASE ##
	###################################################


#' enrich_genes
#'
# add the annotation as additional columns to data, using as the
# column data.genes from data as index into ann.genes to find the
# corresponding annotation.
#
# This function was designed to be used with get_saved_annotation()
# abovem which returns a list with (possibly empty) ENSEMBL and BIOMART 
# annotation tables as members.
#
# use as
#
#   annotation <- get_saved_annotation()
#   ensembl.ann <- annotation$ensembl
#   biomart.ann <- annotation$biomart
#	annotated.data <- enrich_genes(data, "ensembl.gene.id", ensembl.ann, "GENEID")
#	annotated.data <- enrich_genes(data, "ensembl.gene.id", biomart.ann, "ensembl_gene_id")
#
# for ENTREZGENE data
#     annotated.data <- enrich_genes(data, "gene.id", biomart.ann, "ensembl_gene_id")
#
enrich_genes <- function(data, ann, data.genes='ensembl_gene_id', ann.genes=data.genes) {
	
    enriched_genes <- cbind(data, ann[match(data[ , data.genes], ann[ ,ann.genes]), ]) 

    return(enriched_genes)
}


# Annotate genes with ENSEMBL.
#
annotateGenesFromENSEMBL <- function(ensembl.db, gene.ids, save = TRUE, out.dir) {
	
	if ( missing(gene.ids) ) {
		error("Provide either GeneIds OR GeneNames as input")}

	if (all(str_sub(gene.ids, 1, 3) == "ENS")){
        by='GENEID'
    }else{
    	by='GENENAME'}

	# Retrieve Information by GeneName
	ens.ann <- ensembldb::select(	ensembl.db, 
				      	keytype= by, keys = as.character(gene.ids), 
				      	columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION',
						'GENENAME', 'GENEID', 'ENTREZID', 
						'TXNAME', 'TXBIOTYPE', 
						'PROTEINID', 'UNIPROTID'))
		
	# Check if the amount of genes we have is the same as the number of 
	# the annotations that we have extracted
	if ( length(ens.ann[,by]) != length(gene.ids) ) {
		cat(">>> ENSEMBL Annotation does not uniquely match genes\n")
		cat(">>>     Using only one entry per gene\n")
		ens.ann.1 <- ens.ann[ ! duplicated(ens.ann$GENEID), ]

	} else {
		ens.ann.1 <- ens.ann
	}

	if (save == TRUE){

	# SAVE ANNOTATION
	# ---------------
		write.table(ens.ann, file = paste(out.dir, 'ensembl.annotation.txt', sep = '/'), 
					sep = '\t', row.names = T, col.names = T)
		
		write.table(ens.ann.1, file = paste(out.dir, 'ensembl.annotation.uniq.txt', sep = '/'), 
			sep = '\t', row.names = T, col.names = T)
	}
	
	return(ens.ann.1)
}

	###############################################
	## EXTRACT GENE ANNOTATION FROM ORG DATABASE ##
	###############################################

annotateGenesFromORG <- function(org.db, gene.ids, save = TRUE, out.dir) {
	
    if (VERBOSE) cat.nl("Annotating genes from ORG database")
	if ( missing(gene.ids) ) {
		error("Provide either GeneIds OR GeneNames as input")}

	if (all(str_sub(gene.ids, 1, 3) == "ENS")){
        by = 'ENSEMBL'
    } else if (! any(is.na(as.numeric(gene.ids)))){
    	by = 'ENTREZID'
	} else {
		by = 'SYMBOL'}

	# Retrieve Org Package Annotation.
	# We retrieve annotation in 3 different times and then we merge it.
	# This way we reduce computing times.
	org.ann <- AnnotationDbi::select(org.db, 
				    keys = as.character(gene.ids), 
				    columns = c("ENTREZID", "GENENAME", "GENETYPE", "SYMBOL",
								#"ALIAS", "ENZYME", "ACCNUM", "PMID",
								#"UNIPROT", "PROSITE", "PFAM",
								"GO", "GOALL","ONTOLOGYALL", "PATH"
								), 
				    keytype = by, 
				    multiVals = "CharacterList")
					
	org.ann.2 <- AnnotationDbi::select(org.db, 
				    keys = as.character(gene.ids), 
				    columns = c("ENTREZID", #"GENENAME", "GENETYPE", "SYMBOL",
								"ALIAS", "ENZYME", "ACCNUM", "PMID"
								#"UNIPROT", "PROSITE", "PFAM",
								#"GO", "GOALL","ONTOLOGYALL", "PATH"
								), 
				    keytype = by, 
				    multiVals = "CharacterList")
						
	org.ann.3 <- AnnotationDbi::select(org.db, 
				    keys = as.character(gene.ids), 
				    columns = c("ENTREZID", #"GENENAME", "GENETYPE", "SYMBOL",
								#"ALIAS", "ENZYME", "ACCNUM", "PMID",
								"UNIPROT", "PROSITE", "PFAM"
								#"GO", "GOALL","ONTOLOGYALL", "PATH"
								), 
				    keytype = by, 
				    multiVals = "CharacterList")					
					
	for (i in colnames(org.ann.2)){
		org.ann[i] <- org.ann.2[match(org.ann$ENTREZID, org.ann.2$ENTREZID), i]
	}
	for (i in colnames(org.ann.3)){
		org.ann[i] <- org.ann.3[match(org.ann$ENTREZID, org.ann.3$ENTREZID), i]
	}

	# Retrieve Gene Ontology descriptions from GO.db Package
	GO.ann <- AnnotationDbi::select(GO.db,
					keys = as.character(org.ann$GO), 
					columns = c(	"GOID", "TERM", "DEFINITION", "ONTOLOGY"),
					keytype = "GOID")
	
	# Merge both datasets by GO ID.
	GO.ann <- rename(GO.ann, "GO" = "GOID" )
	
	for (i in colnames(GO.ann)){
		org.ann[i] <- GO.ann[match(org.ann$GO, GO.ann$GO), i]
	}
	
	# Check if the amount of genes we have is the same as the number of 
	# the annotations that we have extracted
	org.ann.1 <- org.ann[ ! duplicated(org.ann$ENTREZID), ]
	
	if (save == TRUE){

	# SAVE ANNOTATION
	# ---------------
		write.table(org.ann, file = paste(out.dir, 'orgDb.GO.annotation.txt', sep = '/'), 
					sep = '\t', row.names = T, col.names = T)
				
		write.table(org.ann.1, file = paste(out.dir, 'orgDb.GO.annotation.1st.txt', sep = '/'), 
			sep = '\t', row.names = T, col.names = T)
	}
	
	return(org.ann.1)
}




ds2.analyze.go.representation <- function(ds.data, bm.go.annot, folder) 
{ 
    # Analyze GO representation
    # -------------------------
    #
    # we use the biomaRt database from above
    #	mart.db might have been already assigned above, but not
    #	necessarily
    #ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
    #mart.db <- useMart("ensembl", mart.name)
    #listAttributes(mart.db)

    for (n in names(ds.data)) {
        cat("Processing GOs for", n, '\n')
        sig <- as.data.frame(ds.data[[n]]$signif)
        sig.a <- as.data.frame(ds.data[[n]]$signif.ann)
        filterValues <- rownames(sig)
		
		if (dim(sig)[1] < 1) next	# Skip if no significant changes detected

		if (all(str_sub(filterValues , 1, 3) == "ENS")){
        	by = 'ensembl_gene_id'
    	} else if (! any(is.na(as.numeric(filterValues)))){
    		by = 'entrezgene_id'
		} else {
			by = 'entrezgene_accession'}

        sig[,by] <- filterValues
		
        # this will give us ALL GO annotations
        # we should use the saved searches from above to save bandwidth.
    #    gos <- getBM(attributes=c("ensembl_gene_id", 
    #                              "go_id", "name_1006", 
    #                              "definition_1006" ), 
    #                 mart=mart.db,
    #                 filters=ourFilterType,
    #                 values=filterValues)
        gos <- bm.go.annot[ bm.go.annot[,by] %in% filterValues,
                            c(	"ensembl_gene_id",
								"entrezgene_id",
								"entrezgene_accession",
                              	"go_id", "name_1006", 
                              	"definition_1006" ) ]
		
        sig.a <- merge(sig, gos, by = by)
        table(sig.a$go_id)	# this gives counts, but we'd like to multiply those
                                # counts by log2FC
        
		# Aggregate log2FC by go_id and sum the values
        cat("Top 10 aggregated GO IDs\n")
        print(head(aggregate(sig.a$log2FoldChange, 
                             by = list(Category=sig.a$go_id),
                             FUN=sum), 10))
       
	    # Or with formula interface
        gosums <- aggregate(log2FoldChange ~ go_id, sig.a, sum)
        gosums <- gosums[ order(gosums$log2FoldChange), ]
        
		# Annotate gosums
        gosums.a <- cbind(gosums, bm.annot.1[ match(gosums$go_id, bm.annot.1$go_id, nomatch = NA), ])

        out.dir<- paste(folder, '/DESeq2/go', sep='')
        dir.create(out.dir, showWarning=FALSE)
        out.file <- paste(out.dir, '/GO_', n, "_sum.tab", sep='')
        write.table(gosums.a, out.file, sep='\t')

        goavgs <- aggregate(log2FoldChange ~ go_id, sig.a, mean)
        goavgs <- gosums[ order(gosums$log2FoldChange), ]
        
		# annotate gosums
        goavgs.a <- cbind(goavgs, bm.annot.1[ match(goavgs$go_id, bm.annot.1$go_id, nomatch = NA), ])
        out.file <- paste(out.dir, '/GO_', n, "_average.tab", sep='')
        write.table(goavgs.a, out.file, sep='\t')
    }

}


annotate.go.ancestry <- function (ann.results, l2fc.threshold, comparison, out.dir) {
    res <- ann.results
	
	## Biological Process
    goBPanc <- as.list(GOBPANCESTOR)
    # remove GO terms that do not have any ancestor
    goBPanc <- goBPanc[ ! is.na(goBPanc) ]
	
	## Cellular Component
    goCCanc <- as.list(GOCCANCESTOR)
    # remove GO terms that do not have any ancestor
    goCCanc <- goCCanc[ ! is.na(goCCanc) ]
	
	## Mollecular Fucntion
    goMFanc <- as.list(GOMFANCESTOR)
    # remove GO terms that do not have any ancestor
    goMFanc <- goMFanc[ ! is.na(goMFanc) ]

    res <- res[ abs(res$log2FoldChange) > l2fc.threshold, ] %>% as.data.frame
    go.l2fc <- res[ , c("ensembl_gene_id", "log2FoldChange", "GENENAME", "go_id") ] %>% as.data.frame
	
	if (nrow(res) == 0) {
		cat.warn("\n\tNo significant ontologies found\n\n")
		return()
	}
	
    cat('\n')
    for (i in 1:nrow(res)) {
    	    
        # Skip if ancestry is null
        if ((i %% 100) == 0) cat(".")
        if (is.null(res[i, "go_id"])) next
        if (length(res[i, "go_id"]) == 0) next
		if (is.na(res[i, "go_id"])) next
        if (res[i, "go_id"] == '') next
        
		cat(res[i, "go_id"]," ")
		
        anc <- goBPanc[[ res[i, "go_id"] ]]
		
        #cat(anc)
		if (is.null(anc)) next
		if (length(anc) == 1 && anc == "") next 
        if (length(anc) == 0) next
		#if (is.na(anc)) next
  
        # add  these GoIds to the table of gene/l2FC/GoIds 
        for (anc.goid in anc) {
			if (anc.goid == "" || anc.goid == "all") next 
            go.l2fc <- rbind(go.l2fc, data.frame(
                  ensembl_gene_id = res[i, "ensembl_gene_id"],
                  log2FoldChange = res[i, "log2FoldChange"],
                  GENENAME = res[i, "GENENAME"],
                  go_id = anc.goid))
        }
        
        anc <- goCCanc[[ res[i, "go_id"] ]]
        
		if ( is.null(anc) ) next
        if ( is.na(anc) ) next
		if ( length(anc) == 1 && anc == "" ) next 
        
		# add  these GoIds to the table of gene/l2FC/GoIds 
        for (anc.goid in anc) {
        	if (anc.goid == "" || anc.goid == "all") next
			go.l2fc <- rbind(go.l2fc, data.frame(
                	   ensembl_gene_id=res[i, "ensembl_gene_id"],
                	   log2FoldChange=res[i, "log2FoldChange"],
                	   GENENAME=res[i, "GENENAME"],
                	   go_id=anc.goid))
        }
        
        
        anc <- goMFanc[[ res[i, "go_id"] ]]
        
		if ( is.null(anc) ) next
        if ( is.na(anc) ) next
		if ( length(anc) == 1 && anc == "" ) next
	 
        # add  these GoIds to the table of gene/l2FC/GoIds 
        for (anc.goid in anc) {
        	if (anc.goid == "" || anc.goid == "all") next
            go.l2fc <- rbind(go.l2fc, data.frame(
                  ensembl_gene_id=res[i, "ensembl_gene_id"],
                  log2FoldChange=res[i, "log2FoldChange"],
                  GENENAME=res[i, "GENENAME"],
                  go_id=anc.goid))
        }
    }
	cat('\n')
    
    # Aggregate data by GOID (will result in c(Category, x) columns
    go.l2fc.sum <- aggregate(go.l2fc$log2FoldChange, 
               by=list(Category=go.l2fc$go_id), FUN=sum)
    go.l2fc.avg <- aggregate(go.l2fc$log2FoldChange, 
               by=list(Category=go.l2fc$go_id), FUN=mean)

    # Rename Category, x to go_id, l2fc
    colnames(go.l2fc.sum) <- c('go_id', 'sum.l2fc')
    colnames(go.l2fc.avg) <- c('go_id', 'avg.l2fc')


    # Annotate
    godesc <- AnnotationDbi::select(GO.db, keys = go.l2fc$go_id, 
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


    # Sort by l2FC
    go.l2fc <- go.l2fc[ order(go.l2fc$log2FoldChange, decreasing=T), ]
    go.l2fc.sum <- go.l2fc.sum[ order(go.l2fc.sum$sum.l2fc, decreasing=T), ]
    go.l2fc.avg <- go.l2fc.avg[ order(go.l2fc.avg$avg.l2fc, decreasing=T), ]
	
    # Sort by abs(log2FC)
    go.l2fc.abs <- go.l2fc[ order(abs(go.l2fc$log2FoldChange), decreasing=T), ]
	go.l2fc.sum.abs <- go.l2fc.sum[ order(abs(go.l2fc.sum$sum.l2fc), decreasing=T), ]
    go.l2fc.avg.abs <- go.l2fc.avg[ order(abs(go.l2fc.avg$avg.l2fc), decreasing=T), ]

    # Save
    out.file <- paste(out.dir, '/GOanc_', comparison, ".tab", sep='')
    write.table(go.l2fc, out.file, sep='\t')
    out.file <- paste(out.dir, '/GOanc_', comparison, "_sum.tab", sep='')
    write.table(go.l2fc.sum, out.file, sep='\t')
    out.file <- paste(out.dir, '/GOanc_', comparison, "_average.tab", sep='')
    write.table(go.l2fc.avg, out.file, sep='\t')
	out.file <- paste(out.dir, '/GOanc_', comparison, "_abs.tab", sep='')
    write.table(go.l2fc.abs, out.file, sep='\t')
    out.file <- paste(out.dir, '/GOanc_', comparison, "_sum_abs.tab", sep='')
    write.table(go.l2fc.sum.abs, out.file, sep='\t')
    out.file <- paste(out.dir, '/GOanc_', comparison, "_average_abs.tab", sep='')
    write.table(go.l2fc.avg.abs, out.file, sep='\t')
}



ds2.analyze.go.ancestry <- function(ds.data, l2fc.threshold = 0, folder)
{
    # Reannotate GO with ancestry
    # library(GO.db)

    out.dir <- paste(folder, '/DESeq2/GO_anc', sep = "")
    dir.create(out.dir, showWarning=FALSE)

    for (n in names(ds.data)) {
        cat("\nTracing GO ancestry for", n, '\n')
        res <- ds.data[[n]]$signif.ann

        annotate.go.ancestry(res, l2fc.threshold, comparison = n, out.dir)

        #ans <- continue.on.enter("Press RETURN to continue: ")
        #if (ans == "q") break
    }
}



ds2.analyze.pfam.representation <- function(ds.data, bm.fam.annot, folder)
{
   # Analyze PFAM representation
    # ---------------------------
    #
    # for PFAM, we can use
    #
    # Get PFAM database and use AC (accession) -> DE (definition) mapping
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
               bm.fam.annot[ match(sig$ensembl_gene_id, bm.fam.annot$ensembl_gene_id, nomatch = NA), ])

        pfam.ids <- sig.a$pfam[ ! is.na(sig.a$pfam) ]

        aggregate(sig.a$log2FoldChange, by=list(Category=sig.a$pfam), FUN=sum)

        sig.a <- merge(sig, bm.fam.annot, by="ensembl_gene_id")
        cat("10 most significant genes with PFAM annotation:\n")
        print(head(sig.a, 10))
        
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
        out.file <- paste(folder, '/DESeq2/PFAM_', n, "_sum.tab", sep='')
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
        out.file <- paste(folder, '/DESeq2/PFAM_', n, "_average.tab", sep='')
        write.table(sig.avg, out.file, sep='\t')

        # save also the files sorted by abs(l2fc)
        sig.sum.abs <- sig.sum[ order(abs(sig.sum$sum.l2FC), decreasing=T), ]
        sig.avg.abs <- sig.avg[ order(abs(sig.avg$avg.l2FC), decreasing=T), ]
        out.file <- paste(folder, '/DESeq2/PFAM_', n, "_sum_abs.tab", sep='')
        write.table(sig.sum.abs, out.file, sep='\t')
        out.file <- paste(folder, '/DESeq2/PFAM_', n, "_average_abs.tab", sep='')
        write.table(sig.avg.abs, out.file, sep='\t')

        if (INTERACTIVE) {
            #continue.on.key()
            ans <- continue.on.enter(prompt="Press return to continue ")
            if (ans == "q") break
        }
    }
    cat('\n')
    options(width=80)
}


### XXX JR XXX This should be two separate functions (GO and KEGG)
GO_KEGG_clusterProfiler <- function(ann.data, 
                      ranking.column='lfc',
		     		  max.size=250,
                      out.dir=paste(folder, 'cProf', sep='/'), 
                      out.name='cProf',
                      use.description=TRUE,
                      OrgDb = NULL,
                      kegg_organism = NULL,			# (cjo = coturnix japonica)
                                             # (gga = gallus gallus)
                      top.n = 10,			# Top number of GO/KEGG categories to represent
                      top.biblio = 5,
					  overwrite = FALSE,
                      verbose=FALSE) {
    
    gseaDat <- subset(ann.data, !is.na(ENTREZID))
    ranking <- gseaDat[ , ranking.column]
    names(ranking) <- gseaDat$ENTREZID
    ranking <- na.omit(ranking)

    s.ranking <- sort(ranking, decreasing=T)
    
    if (max.size < 250) min.size=3 else min.size=10
    ### GO annotation
	if (! is.null(OrgDb)){
	
		gse.out.file <- paste(out.dir, '/', out.name, '.topGO.RData', sep='')
        gse.raw.out.file <- paste(out.dir, '/', out.name, '.raw_p.topGO.RData', sep='')
    	
		if (VERBOSE) cat("Doing GO GSEA with Cluster profiler\n")
    	
		if (file.exists(gse.out.file) && ! overwrite) {
            
			# Default is adjusted p
            gse <- readRDS(gse.out.file)
			out.base <- out.name
			
		} else if (file.exists(gse.raw.out.file) && ! overwrite) {

            # try with raw p data 
            gse <- readRDS(gse.raw.out.file)
            out.base <- paste(out.name, '.raw_p', sep='')
			gse.out.file <- gse.raw.out.file
        
		} else {
            # use an empty table so next check can be done
            gse <- table(c())
		
        	gse <- gseGO(geneList=s.ranking, 
                    	 ont ="ALL", 
                    	 keyType = "ENTREZID", 
                    	 minGSSize = min.size, 
                    	 maxGSSize = max.size, 
                    	 pvalueCutoff = 0.05, 
                    	 verbose = TRUE, 
                    	 OrgDb = OrgDb, 
                    	 pAdjustMethod = "fdr")
                    	 #pAdjustMethod = "none")
        	
				out.base <- out.name
			
			if (dim(gse)[1] == 0) {
            	# p.adjusting may have failed to produce any result:
            	# Try without correction issuing a warning
            	gse <- gseGO(geneList=s.ranking, 
                    	 ont ="ALL", 
                    	 keyType = "ENTREZID", 
                    	 minGSSize = min.size, 
                    	 maxGSSize = max.size, 
                    	 pvalueCutoff = 0.05, 
                    	 verbose = TRUE, 
                    	 OrgDb = OrgDb, 
                    	 #pAdjustMethod = "fdr")
                    	 pAdjustMethod = "none")
						 
            	out.base <- paste(out.name, '.raw_p', sep='')
				gse.out.file <- gse.raw.out.file
        	}
			
        	if (dim(gse)[1] > 1) {
            	saveRDS(gse, file = gse.out.file)
            	gse.tab.file <- paste(out.dir, '/', out.base, '.topGO.tab', sep='')
            	write.table(gse, file = gse.tab.file, row.names = T, col.names=T, sep='\t')
        	}
    	}
		
    	# if dim(gse)[1] == 0 then we won't save anything
    	#	hopefully, if re-run again it might work next time by chance
    	# else
    	if (dim(gse)[1] > 0) {
		if (VERBOSE) cat("Plotting GO GSEA with Cluster profiler\n")
        	# DO GENERIC PLOTS
        	# **Dotplot**: For each group shows if it is up or down 
        	# regulated and to which extent. The circle size is proportional to the
        	# size (the number of genes contained) of the group, and the color to the
        	# p-value.
        	out.png <- paste(out.dir, '/', out.base, '.GOdotplot.png', sep='')
        	as.png(
            	dotplot(gse, showCategory = top.n, split=".sign") + facet_grid(.~.sign)
            	, out.png, overwrite = overwrite)
 			
			#dotplot(gse, showCategory = top.n, split=".sign",
			#		label_format = 50, font.size = 10,
			#		title = "GO Enrichment Analysis") + facet_grid(.~.sign)
        
			# **Enrichment map**: organizes terms in a network with edges 
        	# connecting overlapping gene sets (i.e. shows which gene sets contain common
        	# genes). Dot diameter represents the size of the gene set (the number of genes
        	# it contains) and color represents the adjusted probability.
        	out.png <- paste(out.dir, '/', out.base, '.GOemmaplot.png', sep='')
        	as.png(
            	emapplot(pairwise_termsim(gse), showCategory = top.n)
            	, out.png, overwrite = overwrite)

        	# **Ridgeplot**: density plots grouped by gene set depicting
        	# the frequency of fold change values per gene within each set. Helps 
        	# interpret up/down-regulated pathways.
        	out.png <- paste(out.dir, '/', out.base, '.GOridgeplot.png', sep='')
        	as.png(
            	ridgeplot(gse) + labs(x = "Enrichment distribution")
            	, out.png, overwrite = overwrite)

        	# **PubMed trend** of enriched terms
        	# Plots the number/proportion of publications trend based on 
        	# the query result from PubMed Central.
        	out.png <- paste(out.dir, '/', out.base, '.GOpmcplot.png', sep='')
        	cur.year <- as.integer(format(Sys.Date(), "%Y"))
        	terms <- gse$Description[1:top.biblio]
        	as.png(
            	pmcplot(terms, (cur.year-10):(cur.year-1), proportion=FALSE)
            	, out.png, overwrite = overwrite)

        	# DO DETAILED PLOTS
        	for (i in 1:dim(gse)[1] ) {
            	# **GSEA plot**: __Plot of the Running Enrichment Score__ (green
            	# line) for a gene set as the analysis walks down the ranked gene list,
            	# including the location of the maximum enrichment score (the red line).
            	# The black lines in the Running Enrichment Score show where the members of the
            	# gene set appear in the ranked list of genes, indicating the leading edge
            	# subset.
            	# 
            	# __Ranked list metric__ shows the value of the ranking
            	# metric (log2 fold change) as you move down the list of ranked genes. The
            	# ranking metric measures a gene’s correlation with a phenotype.

            	# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
            	out.png <- paste(out.dir, '/', 
                    	out.base, '.GOgseaplot.', i, '.', gse$ID[i], '.png', sep='')
            	as.png(
                	gseaplot(gse, by = "all", title = gse$Description[i], geneSetID = i)
                	, out.png, overwrite = overwrite)

            	# **Category Netplot**: shows the
            	# linkage between genes and biological concepts as a network (helpful to see
            	# which genes are involved in enriched pathways and genes that may belong to
            	# multiple annotation categories).
            	out.png <- paste(out.dir, '/', 
                    	out.base, '.GOcnetplot.', i, '.', gse$ID[i], '.png', sep='')
            	# categorySize can be either 'pvalue' or 'geneNum'
            	as.png(
                	cnetplot(gse, categorySize="pvalue", foldChange=s.ranking, 
                        	 showCategory=i, cex_label_gene = 0.5)
                	, out.png, overwrite = overwrite)
        	}
      	}
    }
    
	### K E G G annotation
    # let's try with KEGG (from the ENTREZID which is the same a ncbi-genid)
    
	if (! is.null(kegg_organism)){
    	
		if (VERBOSE) cat("Doing KEGG GSEA with ClusterProfiler\n")
    	kse.out.file <- paste(out.dir, '/', out.name, '.topKEGG.RData', sep='')
        kse.raw.out.file <- paste(out.dir, '/', out.name, '.raw_p.topKEGG.RData', sep='')
    	
        if (file.exists(kse.out.file) && ! overwrite) {
		
            # default is adjusted p
            kse <- readRDS(kse.out.file)
			out.base <- out.name
			
        } else if (file.exists(kse.raw.out.file) && ! overwrite) {
		
            # try raw p
            kse <- readRDS(kse.raw.out.file)
           	out.base <- paste(out.name, '.raw_p', sep='')
			kse.out.file <- kse.raw.out.file
 
        } else {

            # use an empty table so next check can be done
            kse <- table(c())

        	kse <- gseKEGG(geneList     = s.ranking,
                	   organism     = kegg_organism,
                	   #nPerm        = 10000,
                	   minGSSize    = min.size,
                	   maxGSSize    = max.size,
                	   pvalueCutoff = 0.05,
                	   pAdjustMethod = "fdr",
                	   keyType       = "ncbi-geneid",
                	   nPermSimple = 100000)
			
			out.base <- out.name

        	if (dim(kse)[1] == 0) {
            	# Try without correction issuing a warning
            	kse <- gseKEGG(geneList     = s.ranking,
                    	   organism     = kegg_organism,
                    	   #nPerm        = 10000,
                    	   minGSSize    = min.size,
                    	   maxGSSize    = max.size,
                    	   pvalueCutoff = 0.05,
                    	   pAdjustMethod = "none",
                    	   keyType       = "ncbi-geneid",
                    	   nPermSimple = 100000)
						   
            	out.base <- paste(out.name, '.raw_p', sep='')
				kse.out.file <- kse.raw.out.file
			}
			
			if (dim(kse)[1] > 0) {
            	saveRDS(kse, file = kse.out.file)
            	kse.tab.file <- paste(out.dir, '/', out.base, '.topKEGG.tab', sep='')
            	write.table(kse, file=kse.tab.file, row.names = T, col.names=T, sep='\t')
        	}
    	}
    
	
    	# if dim(kse)[1] == 0 then we won't save anything
    	#	hopefully, if re-run again it might work next time by chance
    	# else
    	if (dim(kse)[1] > 0) {
        	if (VERBOSE) cat("Plotting KEGG GSEA with ClusterProfiler\n")
        	# A *Dotplot*: For each group shows if it is up or down 
        	# regulated and to which extent. The circle size is proportional to the
        	# size (the number of genes contained) of the group, and the color to the
        	# p-value.
        	out.png <- paste(out.dir, '/', out.base, '.KEGGdotplot.png', sep='')
        	as.png( 
            	dotplot(kse, showCategory = top.n, title = "Enriched Pathways" , 
                    	color = "pvalue", split=".sign") + facet_grid(.~.sign)
            	, out.png, overwrite = overwrite)
			
			#dotplot(kse, showCategory = top.n, split=".sign", color = "pvalue", 
			#		label_format = 50, font.size = 10,
			#		title = "KEGG Enrichment Analysis") + facet_grid(.~.sign)
        

        	# B *Enrichment map*: organizes terms in a network with edges 
        	# connecting overlapping gene sets (i.e. shows which gene sets contain common
        	# genes). Dot diameter represents the size of the gene set (the number of genes
        	# it contains) and color represents the adjusted probability.
        	out.png <- paste(out.dir, '/', out.base, '.KEGGemapplot.png', sep='')
        	as.png(
            	emapplot(pairwise_termsim(kse))
            	, out.png, overwrite = overwrite)


        	# C *Ridgeplot*: density plots grouped by gene set depicting
        	# the frequency of fold change values per gene within each set. Helps 
        	# interpret up/down-regulated pathways.
        	out.png <- paste(out.dir, '/', out.base, '.KEGGridgeplot.png', sep='')
        	as.png(
            	ridgeplot(kse, showCategory = top.n, fill = "pvalue") + labs(x = "enrichment distribution")
            	, out.png, overwrite = overwrite)
 			#ridgeplot(kse, showCategory = top.n, fill = "pvalue") + labs(x = "enrichment distribution", font.size = 10, title = "KEGG Enrichment Analysis")
        	
        	# D *Category Netplot*: shows the
        	# linkage between genes and biological concepts as a network (helpful to see
        	# which genes are involved in enriched pathways and genes that may belong to
        	# multiple annotation categories).
        	# categorySize can be either 'pvalue' or 'geneNum'
        	out.png <- paste(out.dir, '/', out.base, '.KEGGcnetplot.png', sep='')
        	as.png(
            	cnetplot(	kse, showCategory = top.n, categorySize = "pvalue",
							cex_label_gene = 0.5, foldChange = s.ranking)
            	, out.png, overwrite = overwrite)

			# PubMed trend of enriched terms
        	# 
        	# Plots the number/proportion of publications trend based on the query result
        	# from PubMed Central.

        	out.png <- paste(out.dir, '/', out.base, '.KEGGpmcplot.png', sep='')
        	cur.year <- as.integer(format(Sys.Date(), "%Y"))
	    	terms <- kse@result$Description[1:top.biblio]
			as.png(
				pmcplot(terms, (cur.year-10):(cur.year-1),
						proportion = FALSE)
            	, out.png, overwrite = overwrite)

        	cur.dir <- getwd()
        	for (i in 1:dim(kse)[1]) {
            	# for each of the pathways in kse
            	# Use the `Gene Set` param for the index in the title, and as the value for geneSetId

            	# GSEA plot 
            	# Plot of the Running Enrichment Score (green
            	# line) for a gene set as the analysis walks down the ranked 
            	# gene list, including the location of the maximum enrichment 
            	# score (the red line). The black lines in the Running Enrichment 
            	# Score show where the members of the gene set appear in the 
            	# ranked list of genes, indicating the leading edge subset.
            	#
            	# The Ranked list metric shows the value of
            	#  the ranking metric (log2 fold change) as you move down the 
            	# list of ranked genes. The ranking metric measures a gene’s 
            	# correlation with a phenotype.
            	out.png <- paste(out.dir, '/', 
                	out.base, '.KEGGgseaplot.', i, '.', kse$ID[i], '.png', sep='')
            	as.png(
            	  gseaplot(kse, by = "all", title = kse$Description[i], geneSetID = i)
            	  , out.png, overwrite = overwrite)

            	# **Pathview**
            	# This will create a PNG and a __different__ PDF of the enriched 
            	# KEGG pathway in the current working directory.
            	wd <- setwd(out.dir)	# change to the appropriate directory
				
            	# Produce the native KEGG plot (PNG)
				if (! file.exists(paste(kse$ID[i], "pathview", "png", sep = "."))){
            	  tryCatch(
                    { dme <- pathview(gene.data=s.ranking, 
                                    pathway.id=kse$ID[i], 
                                    species = kegg_organism) },
                    error=function(err) { cat.err("Error in pathview(PNG):", i,"id =",kse$ID[i],'\n'); return}
                  )
	            	# Produce a PDF plot with the graphviz layout engine
            	  tryCatch(
    	        	{ dme <- pathview(gene.data=s.ranking, 
                                    pathway.id=kse$ID[i], 
                                    species = kegg_organism, 
                                    kegg.native = F) },
                    error=function(err) { cat.err("Error in pathview(PDF):", i,"id =",kse$ID[i],'\n'); return }
                  )
                }
				wd <- setwd(cur.dir) # and back
        	}
    	}
	}
}



# -------------------------------
# Do Gene Set Enrichment Analysis
# -------------------------------

GO_fgsea <- function (ann.data, 
		      		  ranking.column = 'lfc',
                      max.size = 250,
                      out.dir = paste(folder, 'go_fgsea', sep='/'), 
                      out.name = 'fgsea',
                      use.description = TRUE,
                      top.n = 10,
                      top.biblio = 5,
                      verbose = FALSE,
					  interactive = FALSE) {
    
	## Do GSEA on GO terms using fgsea

    # Rank all genes on their fold change.
    #	Here we exclude genes for which we have no EntrezID and
    #	use shrunk LFC values
    gseaDat <- filter(ann.data, !is.na(ENTREZID))
    #gseaDat <- filter(ann.shrunk.lfc, !is.na(GENEID))

    ranking <- gseaDat[ , ranking.column]
    #names(ranking) <- gseaDat$ENTREZID
    names(ranking) <- gseaDat$GENEID
    cat("Top significant genes by", ranking.column, ":\n")
    print(head(ranking))
    #uranks <- ranking[!duplicated(sort(names(ranking)))]

    # Plot all the ranked fold changes
    out.png <- paste(out.dir, '/', out.name, '.barplot.png', sep='')
    # BARPLOT return a numeric vector (or matrix) gicing the 
    # coordinates of ALL the bar midpoints drawn for adding to the
    # graph, this will show when as.png calls print(barcode...), so
    # we need to eat its return values
    # We might try invisible(barplot(...)) instead or maybe
    # capture.output(barplot(...), file='/dev/null'))
    as.png( {x<-barplot(sort(ranking, decreasing=T), 
                    ylab = "Log 2 Fold Change",
					main = "Distribution of LFC",
					xlab = "Genes", 
                    density = TRUE,
					cex.names = 0.5); rm(x)}, file=out.png)

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
    if (verbose) cat.info("Loading GO table\n")
    if ( ! exists(substitute(bm.go.annot)) ) {
        bm.go.annot <- read.table(paste(folder, '/annotation/biomart.go.tab', sep=''), 
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
        if (verbose) cat('Loading precomputed FGSEA data\n')
        load(file=out.file)
    } else {
        if (verbose) cat('Computing and saving FGSEA data\n')
        fgseaRes <- fgsea::fgsea(pathways=pathways.go, 
                        		 stats=ranking, 
				        		 minSize=10, 
				        		 maxSize=max.size, 
                        		 nPermSimple=100000,
                        		 scoreType='std'
                        		 )
        save(fgseaRes, file=out.file)
    }
    if (verbose)
        cat.info("First 10 GO terms")
        print(head(fgseaRes[order(padj, -abs(NES)), ], n=10))

    if (interactive) {
        # plot enrichment score
        sorted.fgsea.res <- fgseaRes[order(padj, -abs(NES)), ]
        sfr.names <- sorted.fgsea.res$pathway
        for (i in 1:top.n) {
            if (use.description == FALSE)
                descr <- bm.go.annot[bm.go.annot$go_id == sfr.names[i], "name_1006"]
            else
                descr <- sfr.names[i]
            print(
                fgsea::plotEnrichment(pathways.go[[ sfr.names[i] ]], ranking) +
                    labs(title=descr)
                )
            ans <- readline("Press RETURN to continue: ")
            if (ans == "q") break
        }
    }
    if (verbose) cat('Selecting top up- and down-regulated sets\n')
    # gsea table plot of top.n gene families
    #   top_n() is now deprecated
    topUp <- fgseaRes %>%
        filter(ES > 0) %>%
        top_n(n = top.n, col="padj")

    topDown <- fgseaRes %>%
        filter(ES < 0) %>%
        top_n(n = top.n, col="padj")

    topPathways <- bind_rows(topUp, topDown) %>%
        arrange(-ES)
    # last resort when "pos" is used    
    #topPathways <- sorted.fgsea.res[1:top.n, ]

    # do the plots and save descriptions
    if (verbose) cat('Plotting top', top.n, 'up- and down-regulated sets\n')
    out.file <- paste(out.dir, "/topUp.", top.n, '.txt', sep='')
    
	for (i in 1:dim(topUp)[1]) {
        
		name <- topUp$pathway[i]
        
		if (use.description == FALSE){
            descr <- unique(bm.go.annot[bm.go.annot$go_id == name, "name_1006"])
        } else {
            descr <- name
        }
		cat(i, descr, file=out.file, '\n', sep='\t', append = TRUE)
 
        #print(name)
        #out.png <- paste(out.dir, '/', out.name, '.topUp', i, '.enrichment.png', sep='')
        out.png <- sprintf("%s/%s.topUp.%03d.enrichment.png", 
        	out.dir, out.name, i)
        out.plt <-
            fgsea::plotEnrichment(pathways.go[[ name ]] , ranking) +
            labs(title=descr)
        ggsave(out.png, out.plt)
    }
    
	out.file <- paste(out.dir, "/topDn.", top.n, '.txt', sep='')
    
	for (i in 1:dim(topDown)[1]) {
        
		name <- topDown$pathway[i]
		
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
            fgsea::plotEnrichment(pathways.go[[ name ]] , ranking) +
            labs(title=descr)
        ggsave(out.png, out.plt)
    }

    out.png <- paste(out.dir, '/', out.name, '.GSEAtable.png', sep='')
    as.png(
        fgsea::plotGseaTable(pathways.go[topPathways$pathway], 
                      ranking, 
                      fgseaRes, 
                      gseaParam = 0.5)
    , out.png, width=1024, height=100*top.n, overwrite=TRUE)
    
    
    # and now do a plot of the interest in citations during
    # the last ten years for the top 5 sets
    if (verbose) cat('Getting and plotting PubMed citations\n')
    out.png <- paste(out.dir, '/', out.name, '.EUPMC.png', sep='')
    cur.year <- as.integer(format(Sys.Date(), "%Y"))
    terms <- topPathways$pathway[1:top.biblio]
    as.png(
        pmcplot(terms, (cur.year-10):(cur.year-1), proportion=FALSE),
        out.png)

    
    # finally, return fgseaRes
    return(fgseaRes)
}



#' get.stringsdb.ppi
#'
#' get a full PPI network for a given organism
#' the network will have rownames and colnames set to ENSEMBL_gene_ids
#' and so it is subsettable selecting detected genes
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
get.stringsdb.ppi <- function(ncbi.taxid=9606, mart.name='hsapiens_gene_ensembl', save.dir='.')
{
    ppi.file=file.path(save.dir, '/STRINGdb.PPI.rds')
    if (file.exists(ppi.file)) {
        ppi=readRDS(ppi.file)
    } else {
        # 1. getSTRINGdb for the organism
        string_db <- STRINGdb$new(species=as.numeric(ncbi.taxid))
        organism_graph <- string_db$get_graph()

        # 2. get edges with high confidence score
        edge.scores <- E(organism_graph)$combined_score
        ninetieth.percentile <- quantile(edge.scores, 0.9)
        thresh <- data.frame(name='90th percentile',
                             val=ninetieth.percentile)
        organism_graph <- subgraph_from_edges(organism_graph,
                  E(organism_graph)[combined_score > ninetieth.percentile])

        # 3. create adjacency matrix
        adj_matrix <- as_adjacency_matrix(organism_graph)

        ### extract protein ids from the organism network
        protein_ids <- sapply(strsplit(rownames(adj_matrix), '\\.'),
                              function(x) x[2])


        # 4. map gene ids to protein ids

        ### get gene/protein ids via Biomart
        mart=useMart(biomart='ENSEMBL_MART_ENSEMBL',
                     dataset=mart.name)     # note that there is also a global mart.name

        ### get protein to gene id mappings
        mart_results <- getBM(attributes = c("ensembl_gene_id",
                                             "ensembl_peptide_id"),
                              filters = "ensembl_peptide_id", values = protein_ids,
                              mart = mart)    

        ### replace protein ids with gene ids
        ix <- match(protein_ids, mart_results$ensembl_peptide_id)
        ix <- ix[!is.na(ix)]
        # we could also use the ann.data columne ensembl_peptide_id
        # and ensembl_gne_id instead
        # ix <- match(protein_ids, ann.data$ensembl_peptide_id)
        # ix <- ix[! is.na(ix) ]

        # convert from protein_ids to ensembl_gene_ids
        # and use to label rows and cols of adj_matrix
        newnames <- protein_ids
        newnames[match(mart_results[ix,'ensembl_peptide_id'], newnames)] <-
            mart_results[ix, 'ensembl_gene_id']
        rownames(adj_matrix) <- newnames
        colnames(adj_matrix) <- newnames

        # remove duplicate names and names with no interactions
        ppi <- adj_matrix[!duplicated(newnames), !duplicated(newnames)]
        nullrows <- Matrix::rowSums(ppi)==0
        ppi <- ppi[!nullrows,!nullrows] ## ppi is the network with ensembl gene ids}

        saveRDS(ppi, file=file.path(save.dir, '/STRINGdb.PPI.rds'))
        # reload with ppi <- readRDS(file=...)
        save(ppi, file=file.path(save.dir, 'STRINGdb.PPI.RData'))
        # reload with read(file=...)
        write.csv(as.matrix(ppi), file=file.path(save.dir, 'STRINGdb.ppi.csv'),
	      row.names=T, col.names=T, quote=F)
        # reload with read.csv()
        write.table(as.matrix(ppi), file=file.path(save.dir, 'STRINGdb.PPI.tab'),
	      sep='\t', row.names=T, col.names=T, quote=F)
        #readr::write_delim(as.matrix(ppi), 
        #    file=file.path(save.dir, 'STRINGdb.PPI.txt'), 
	    #    delim='\t')
        # reload with read.table(row.names=1, header=T,file=...) or read_delim()
    }
    return(ppi)
}
