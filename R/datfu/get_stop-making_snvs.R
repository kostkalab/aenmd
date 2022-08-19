
library(AnnotationHub)
library(BSgenome.Hsapiens.NCBI.GRCh38)
AH_VERSION <- 'AH98047'


#- Get the data
ah   <- AnnotationHub()
edb  <- ah[[AH_VERSION]]

txs <- transcripts(edb,TxBiotypeFilter("protein_coding"),
                   columns = c("protein_id", "tx_cds_seq_start", "tx_cds_seq_end"))

get_stop_making_snvs <- function(txname, txs, verb = FALSE){
#-----------------------------------------------------------

    if(verb) message(txname)
    #- only use transcripts we care about
    if( is.null( seq <- get0(txname, future::value(._EA_cds_env)) )){
        return(NULL)
    }

    #- character vector with codons in cds
    cc_vec <- Biostrings::codons(seq) |> as.character()

    #- patterns for stop-codon-creations snvs
    # cdn_stop   <- c("TAA","TAG", "TGA")
    pat_new_stop_a_2 <- c("^T[C,G,T]A$|^T[C,G,T]G$")
    pat_new_stop_a_3 <- c("^TA[C,G,T]$|^TG[C,G,T]$")
    pat_new_stop_g_2 <- c("^T[A,C,T]A$")
    pat_new_stop_g_3 <- c("^TA[A,C,T]$")
    pat_new_stop_t_1 <- c("^[A,C,G]AA$|^[A,C,G]AG$|^[A,C,G]GA$")

    #- coordinates of new stop codons in cds-codon space, with protein name
    ind_a_2 <- which(! (stringi::stri_match_all_regex(cc_vec, pat_new_stop_a_2) |> unlist() |> is.na() )) |> sort() |> unique()
    ind_a_3 <- which(! (stringi::stri_match_all_regex(cc_vec, pat_new_stop_a_3) |> unlist() |> is.na() )) |> sort() |> unique()
    ind_g_2 <- which(! (stringi::stri_match_all_regex(cc_vec, pat_new_stop_g_2) |> unlist() |> is.na() )) |> sort() |> unique()
    ind_g_3 <- which(! (stringi::stri_match_all_regex(cc_vec, pat_new_stop_g_3) |> unlist() |> is.na() )) |> sort() |> unique()
    ind_t_1 <- which(! (stringi::stri_match_all_regex(cc_vec, pat_new_stop_t_1) |> unlist() |> is.na() )) |> sort() |> unique()
    rng_a_2 <- IRanges(ind_a_2, ind_a_2)
    rng_a_3 <- IRanges(ind_a_3, ind_a_3)
    rng_g_2 <- IRanges(ind_g_2, ind_g_2)
    rng_g_3 <- IRanges(ind_g_3, ind_g_3)
    rng_t_1 <- IRanges(ind_t_1, ind_t_1)

    rng <- c(rng_a_2, rng_a_3, rng_g_2, rng_g_3, rng_t_1)
    mcols(rng)$protein_id <- txs[txname]$protein_id

    #- We are done here if we don't find stop-creating snvs
    if(length(rng) == 0) return(NULL)

    #- remove the "orginal" stop coding, if we have it
    fix   <- c(rep('center',length(rng_a_2)),
               rep('end',   length(rng_a_3)),
               rep('center',length(rng_g_2)),
               rep('end'   ,length(rng_g_3)),
               rep('start', length(rng_t_1)))

    alt_nuc <- c(rep('A', length(rng_a_2)),
                 rep('A', length(rng_a_3)),
                 rep('G', length(rng_g_2)),
                 rep('G', length(rng_g_3)),
                 rep('T', length(rng_t_1)))
    cmp <- c("T","G","C","A")
    names(cmp) <- c("A","C","G","T")
    if(strand(txs[txname]) |> as.character() == '-'){
        alt_nuc <- cmp[alt_nuc]
        names(alt_nuc) <- NULL
    }
    mcols(rng)$alt_nuc <- alt_nuc

    if(max(start(rng))==length(cc_vec)){
        ind <- which(start(rng) == length(cc_vec))
        rng <- rng[-ind]
        fix <- fix[-ind]
    }

    #- Again done, if the only stop-creating subs overlapped the original stop
    if(length(rng) ==0) return(NULL)

    #- get tx coordinates
    rng_tx <- ensembldb:::.proteinCoordsToTx(rng)
    #- don't need to add 5'UTR, because our exons are gotten as CDS
    #if(strand(txs[txname]) |> as.character() == '-'){
    #    offset <- end(txs[txname]) - txs[txname]$tx_cds_seq_end
    #} else {
    #    offset <- txs[txname]$tx_cds_seq_start - start(txs[txname])
    #}

    offset <- 0
    rng_tx <- shift(rng_tx, shift = offset)
    rng_tx <- IRanges::resize(rng_tx,width=1,fix=fix)
    mcols(rng_tx)$tx_id <- txname

    tmp <- paste(cc_vec, collapse="") |> strsplit(split="")
    ref_nuc <- tmp[[1]][start(rng_tx)]
    if(strand(txs[txname]) |> as.character() == '-'){
        ref_nuc <- cmp[ref_nuc]
        names(alt_nuc) <- NULL
    }
    mcols(rng_tx)$ref_nuc <- ref_nuc

    #- map transcript coordinates to genome
    #  adapted from ensembldb package, but this is way easier as width(rng_tx) == 1 always
    exn     <- future::value(._EA_exn_env)[[txname]]
    exn_spl <- ensembldb:::.splice(exn)
    ov      <- IRanges::findOverlaps(rng_tx, exn_spl)
    qH      <- S4Vectors::queryHits(ov)
    sH      <- S4Vectors::subjectHits(ov)

    #- how deep into each exon are we?
    offset <- start(rng_tx[qH]) - start(exn_spl[sH])
    if(strand(txs[txname]) |> as.character() == '-'){
       g_pos <- end(exn[sH]) - offset
    } else {
       g_pos <- start(exn[sH]) + offset
    }
    res_grng <- GRanges( seqnames(txs[txname]) |> as.character(),IRanges(g_pos,g_pos))
    #- check the reference nucleotides match
    ref_nuc  <- getSeq(Hsapiens,res_grng) |> as.character()
    if(! all(ref_nuc == mcols(rng_tx)$ref_nuc)) stop("internal error.")

    keys <- paste0( seqnames(res_grng) |> as.character(),":",
                    start(res_grng),"|",
                    mcols(rng_tx)$ref_nuc,"|",
                    mcols(rng_tx)$alt_nuc)

    return(keys)

    ###- we can't test like so, because of the UTRs:
    ##mpd_g <- transcriptToGenome(rng_tx, edb, id='tx_id') |> GRangesList() |> unlist()
    ##ref_nuc2  <- getSeq(Hsapiens,mpd_g) |> as.character()

}


get_stop_making_snvs(names(txs)[2],txs)

tictoc::tic(); tmp <- sapply(names(txs)[1:100], function(x) get_stop_making_snvs(x,txs)) ; tictoc::toc()
#- ~ 8h on one core
res <- tmp[!lapply(tmp,is.null) |> unlist()]

library(BiocParallel)
tictoc::tic(); tmp <- bplapply(names(txs)[1:100],function(x) get_stop_making_snvs(x,txs), BPPARAM = MulticoreParam(workers = 12,progressbar = TRUE)); tictoc::toc()
#- ~2h on eight cores
tictoc::tic(); tmp <- bplapply(names(txs),function(x) get_stop_making_snvs(x,txs), BPPARAM = MulticoreParam(workers = 12,progressbar = TRUE)); tictoc::toc()
saveRDS(tmp, file="./foo.rds")

table(lapply(tmp, is.null) |> unlist())                              #- 30,572 transcripts with stp-gains
res <- tmp[!lapply(tmp,is.null) |> unlist()] |> unlist() |> unique() #- 3,505,838 stop-gain-making snvs

#- just put keys in an environment for lookup now
snv_ptc_env  <- new.env(hash=TRUE, parent = emptyenv(), size = floor(length(res)*1.1))
pbapply::pbsapply(res, function(key) { snv_ptc_env[[key]] <- TRUE; return(invisible()) })

#- takes long to save...
saveRDS(snv_ptc_env, file="./inst/extdata/env_ensdb_v105_fil_all-stop-making-snvs.rds")

#- TODO: add rule annotation here; then we can have a for each variant with ENST and names and rules as values.




