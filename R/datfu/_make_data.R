
#- THIS FILE MAKES DATA aenmd needs
#- It is only needed to create the package.
#- Currently using:
#
#  AnnotationHub / Ensembl v. 105 ("AH98047")

library(Biostrings)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(ensembldb)
library(AnnotationHub)

AH_VERSION <- 'AH98047'

#- Get the data
ah   <- AnnotationHub()
edb  <- ah[[AH_VERSION]]

#- COMPILE TRANSCRIPT SET
#------------------------

#- all protein-coding txs on standard chromosomes
FilterList <- list(     SeqNameFilter(as.character(c(1:22,"X","Y","MT"))),
                        TxBiotypeFilter("protein_coding"))
txs <- ensembldb::transcripts(edb, filter = FilterList,
                               columns = c( "tx_id_version","gene_id", "tx_support_level", "tx_is_canonical"))

#- corresponding exon ranges
all_exn_rng   <- ensembldb::cdsBy(edb, "tx", filter = TxIdFilter(names(txs)));
txs           <- txs[names(all_exn_rng),]
txs$num_exons <- elementNROWS(all_exn_rng)

#- tsl 1 transcripts plus
#- recover 1-exon transcripts withou transcript support level info
ind <- txs$tx_support_level == 1
ind <- ind | ((txs$num_exons == 1) & (is.na(txs$tx_support_level)))

txs_used     <- txs[ sort(which(ind))  ]
exn_rng_used <- all_exn_rng[names(txs_used)]
cds_used     <- getSeq(Hsapiens, exn_rng_used)

#- FILTER CDS must be divisible by three
#- some of these sequences are sort of incomplete
#  in this case the cds is not divisible by three.
#  we will exclude these but keep track.
out <- which(width(exn_rng_used) |> lapply(sum) |> unlist() %% 3 != 0) |> sort()
# length(out)
# 927
write.table(names(out) |> sort(), file="../../inst/extdata/enst_length_excluded.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
txs_used     <- txs_used[ -out  ]
exn_rng_used <- exn_rng_used[ -out  ]
cds_used     <- cds_used[ -out  ]

#- TODO:
#- FILTER CDS must have at least 3 codons (start, content, stop)

#- annotate splice regions based on:
#  3bp into the exon and 8bp into the intron;

#- exons (with utrs) ; but for ANY protein-coding transcript on standard chromosomes
ex   <- ensembldb::exonsBy(edb, "tx", filter = FilterList, columns = c("exon_id", "tx_name"))
ex_s <- resize(ex, width=3, fix = "start")
ex_e <- resize(ex, width=3, fix = "end")
sr_s <- resize(ex_s, width = 11, fix = "end")
sr_e <- resize(ex_e, width = 11, fix = "start")
exn_spl_rng <-  pc(sr_s,sr_e)  %>% sort()

#- aggregate for all txs_used ALL overlapping splice regions
ov <- findOverlaps(txs_used, exn_spl_rng)

afu <- function(ind){
	#tx_nme                <- names(txs_used)[ind]
    sh                    <- subjectHits(ov)[which(queryHits(ov)==ind)]
	res                   <- exn_spl_rng[sh] %>% unlist() %>% sort()
	mcols(res)$tx_biotype <- NULL
	return(res)
}

#- FIXME: slow, there needs to be a better way.
spl_rng_used        <- sapply(seq_len(length(txs_used)), afu) %>% GRangesList()
names(spl_rng_used) <- names(txs_used)

#- Could keep these for documentation puproses; don't need them for the package
if(FALSE){
    saveRDS(txs_used,         "../../dat/renvs/ensdb_v105_txs-fil.rds") #- 31,506 ENSEMBL 105
    saveRDS(exn_rng_used,     "../../dat/renvs/ensdb_v105_exns-rng-fil.rds")
    saveRDS(spl_rng_used,     "../../dat/renvs/ensdb_v105_spl-rng-fil.rds")
    saveRDS(cds_used,         "../../dat/renvs/ensdb_v105_cds-fil.rds")
}

#- NOTE
#======
#  It appears several orders of magnitude faster to
#  look up GRanges in an environment compared with a GRangesList.
#  Therefore we make these.
#  Also, we pre-make all the reference sequences.
#  Takes long to make, though.

#- Empty envs
exon_env   <- new.env(hash=TRUE, parent = emptyenv())
cds_env    <- new.env(hash=TRUE, parent = emptyenv())
splice_env <- new.env(hash=TRUE, parent = emptyenv())

#- Fill envs with content
for( i in seq_len(length(txs_used))){
        if(i %% 10 == 0) message( i )
        txn               <- names(txs_used)[i]
        exns              <- exn_rng_used[[txn]]
        exon_env[[txn]]   <- exns
        seq_ref_exns      <- cds_used[[txn]]
        seq_ref_c         <- unlist(seq_ref_exns)
        cds_env[[txn]]    <- seq_ref_c
    	spl               <- spl_rng_used[[txn]]
	    splice_env[[txn]] <- spl
}

#- we put the environments into inst/extdata
if(TRUE){
    saveRDS(exon_env,     file = "../../inst/extdata/env_ensdb_v105_exns_byTx_fil.rds")
    saveRDS(cds_env,      file = "../../inst/extdata/env_ensdb_v105_seqs_byTx_fil.rds")
    saveRDS(splice_env,   file = "../../inst/extdata/env_ensdb_v105_splc_byTx_fil.rds")
    saveRDS(spl_rng_used, file = "../../inst/extdata/grl_ensdb_v105_splc_byTx_fil.rds")
    saveRDS(exn_rng_used, file = "../../inst/extdata/grl_ensdb_v105_exns_byTx_fil.rds")
    saveRDS(txs_used,     file = "../../inst/extdata/grl_ensdb_v105_trnscrpts_fil.rds")
}





