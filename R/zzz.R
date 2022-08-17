#- Load processed ENSEMBL annotation
#-----------------------------------

._EA_exn_env <- NULL
._EA_cds_env <- NULL
._EA_spl_env <- NULL
._EA_spl_grl <- NULL
._EA_txs_grl <- NULL

.onLoad <- function(libname, pkgname) {
    #- load envs containing data we need
    packageStartupMessage("\nPackage aenmd /aɪ-iː-ɛn-ɛm-diː/")
    packageStartupMessage("Annotating escape from NMD")
    packageStartupMessage("==========================\n")
    packageStartupMessage("Loading transcript annotation...")
    prefix       <- system.file("extdata", package = "aenmd")
    ._EA_exn_env <<- readRDS(paste0(prefix,'/','env_ensdb_v105_exns_byTx_fil.rds'))
    ._EA_cds_env <<- readRDS(paste0(prefix,'/','env_ensdb_v105_seqs_byTx_fil.rds'))
    ._EA_spl_env <<- readRDS(paste0(prefix,'/','env_ensdb_v105_splc_byTx_fil.rds'))
    ._EA_spl_grl <<- readRDS(paste0(prefix,'/','grl_ensdb_v105_splc_byTx_fil.rds'))
    ._EA_txs_grl <<- readRDS(paste0(prefix,'/','grl_ensdb_v105_trnscrpts_fil.rds'))

    packageStartupMessage("done.\n")
}

