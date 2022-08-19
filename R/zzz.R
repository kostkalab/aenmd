#- Load processed ENSEMBL annotation
#-----------------------------------

._EA_exn_env <- NULL #- exons we use
._EA_cds_env <- NULL #- cds sequences we use
._EA_spl_env <- NULL #- splice regions we use
._EA_spl_grl <- NULL #- splice regions we use as grl
._EA_exn_grl <- NULL #- exons we use as grl
._EA_txs_grl <- NULL #- transcripts we use as grl
._EA_snv_env <- NULL #- SNVs that caust stop-gains in transcripts we use

.onLoad <- function(libname, pkgname) {
    #- load envs containing data we need
    packageStartupMessage("\n--------------------------------")
    packageStartupMessage("Package: aenmd /e\u026A-i\u02D0-\u025Bn-\u025Bm-di\u02D0/")
    packageStartupMessage("(annotating escape from NMD)")
    packageStartupMessage("--------------------------------\n")
    prefix       <- system.file("extdata", package = "aenmd")
    ._EA_exn_env <<- future::future({readRDS(paste0(prefix,'/','env_ensdb_v105_exns_byTx_fil.rds'))}, lazy = TRUE)
    ._EA_cds_env <<- future::future({readRDS(paste0(prefix,'/','env_ensdb_v105_seqs_byTx_fil.rds'))}, lazy = TRUE)
    ._EA_spl_env <<- future::future({readRDS(paste0(prefix,'/','env_ensdb_v105_splc_byTx_fil.rds'))}, lazy = TRUE)
    ._EA_spl_grl <<- future::future({readRDS(paste0(prefix,'/','grl_ensdb_v105_splc_byTx_fil.rds'))}, lazy = TRUE)
    ._EA_exn_grl <<- future::future({readRDS(paste0(prefix,'/','grl_ensdb_v105_exns_byTx_fil.rds'))}, lazy = TRUE)
    ._EA_txs_grl <<- future::future({readRDS(paste0(prefix,'/','grl_ensdb_v105_trnscrpts_fil.rds'))}, lazy = TRUE)
    ._EA_snv_env <<- future::future({readRDS(paste0(prefix,'/','env_ensdb_v105_fil_all-stop-making-snvs.rds'))}, lazy = TRUE)
}

