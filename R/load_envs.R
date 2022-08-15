
.onLoad <- function(libname, pkgname) {
    #- load envs containing data we need
    prefix     <- system.file("extdata", package = "aenmd")
    exon_env   <- readRDS(paste0(prefix,'/','env_ensdb_v105_exns_byTx_fil.rds'))
    cds_env    <- readRDS(paste0(prefix,'/','env_ensdb_v105_seqs_byTx_fil.rds'))
    splice_env <- readRDS(paste0(prefix,'/','env_ensdb_v105_splc_byTx_fil.rds'))
}
