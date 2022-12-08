#- Load processed ENSEMBL annotation
#-----------------------------------

._EA_exn_env <- NULL #- exons we use
._EA_cds_env <- NULL #- cds sequences we use
._EA_spl_env <- NULL #- splice regions we use
._EA_spl_grl <- NULL #- splice regions we use as grl
._EA_exn_grl <- NULL #- exons we use as grl
._EA_txs_grl <- NULL #- transcripts we use as grl
._EA_snv_tri <- NULL #- SNVs that caust stop-gains in transcripts we use
._EA_set_env <- NULL #- list of single exon transcripts ; TODO: make me a tri

#- loading this data package will actually populate the objects above, 
#  which are needed for the aenmd package to function.
._EA_dataPackage_name <- 'aenmd.data.ensdb.v105'

.onLoad <- function(libname, pkgname) {
 
    tmp <- paste("-> using data from package:", ._EA_dataPackage_name)
    lne <- paste(rep('-',max(nchar(tmp), 33)),collapse="")

    #- load envs containing data we need
    packageStartupMessage("")
    packageStartupMessage(lne)
    packageStartupMessage("Package: aenmd /e\u026A-i\u02D0-\u025Bn-\u025Bm-di\u02D0/")
    packageStartupMessage("(annotating escape from NMD)")
    packageStartupMessage(paste0("Version: ", utils::packageVersion('aenmd')))
    packageStartupMessage(lne)

    packageStartupMessage(tmp) 
   
    ._EA_exn_env <<- get0( '._EA_exn_env' , envir = asNamespace(._EA_dataPackage_name))
    ._EA_cds_env <<- get0( '._EA_cds_env' , envir = asNamespace(._EA_dataPackage_name))
    ._EA_spl_env <<- get0( '._EA_spl_env' , envir = asNamespace(._EA_dataPackage_name))
    ._EA_spl_grl <<- get0( '._EA_spl_grl' , envir = asNamespace(._EA_dataPackage_name))
    ._EA_exn_grl <<- get0( '._EA_exn_grl' , envir = asNamespace(._EA_dataPackage_name))
    ._EA_txs_grl <<- get0( '._EA_txs_grl' , envir = asNamespace(._EA_dataPackage_name))
    ._EA_snv_tri <<- get0( '._EA_snv_tri' , envir = asNamespace(._EA_dataPackage_name))
    ._EA_set_env <<- get0( '._EA_set_env' , envir = asNamespace(._EA_dataPackage_name))
    packageStartupMessage(lne)
    rm(tmp)
}

