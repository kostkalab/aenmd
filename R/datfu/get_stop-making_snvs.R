library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
seqlevelsStyle(Hsapiens) <- 'NCBI'

#- collect futures b/c we parallelize later on.
#. might want to get rid of fturures; they are also slow.

future::value(._EA_cds_env)
future::value(._EA_exn_env)

get_stop_making_snvs <- function(txname, verb = FALSE){
    #- get sequence and exons
    #------------------------
    if( is.null( seq <- get0(txname, future::value(._EA_cds_env)) )){
        return(NULL)
    }

    if( is.null( exn <- get0(txname, future::value(._EA_exn_env)) )){
        return(NULL)
    }

    #- vector of codons / nucleotides
    #--------------------------------
    cdn_vec <- Biostrings::codons(seq) |> as.character() #- this is a bottleneck currently (?!)
    nuc_vec <- as.character(seq) |> strsplit(split="") |> unlist()

    #- match potential new stop codons
    #---------------------------------
    #- patterns for stop-codon-creations snvs
    # cdn_stop   <- c("TAA","TAG", "TGA")
    pat_new_stop_a_2 <- c("^T[C,G,T]A$|^T[C,G,T]G$")
    pat_new_stop_a_3 <- c("^TA[C,G,T]$|^TG[C,G,T]$")
    pat_new_stop_g_2 <- c("^T[A,C,T]A$")
    pat_new_stop_g_3 <- c("^TA[A,C,T]$")
    pat_new_stop_t_1 <- c("^[A,C,G]AA$|^[A,C,G]AG$|^[A,C,G]GA$")

    #- coordinates of new stop codons in cds-codon space, with protein name
    ind_a_2 <- which(! (stringi::stri_match_all_regex(cdn_vec, pat_new_stop_a_2) |> unlist() |> is.na() )) |> sort() |> unique()
    ind_a_3 <- which(! (stringi::stri_match_all_regex(cdn_vec, pat_new_stop_a_3) |> unlist() |> is.na() )) |> sort() |> unique()
    ind_g_2 <- which(! (stringi::stri_match_all_regex(cdn_vec, pat_new_stop_g_2) |> unlist() |> is.na() )) |> sort() |> unique()
    ind_g_3 <- which(! (stringi::stri_match_all_regex(cdn_vec, pat_new_stop_g_3) |> unlist() |> is.na() )) |> sort() |> unique()
    ind_t_1 <- which(! (stringi::stri_match_all_regex(cdn_vec, pat_new_stop_t_1) |> unlist() |> is.na() )) |> sort() |> unique()

    #- in one vector
    iinds <- c(ind_a_2,
               ind_a_3,
               ind_g_2,
               ind_g_3,
               ind_t_1)

    #- reference and alternative nucleotides (FORWARD strand)
    #--------------------------------------------------------
    ioffs <- c( rep( 1, length(ind_a_2) ),
                rep( 2, length(ind_a_3) ),
                rep( 1, length(ind_g_2) ),
                rep( 2, length(ind_g_3) ),
                rep( 0, length(ind_t_1) ))

    ref_nuc <- nuc_vec[(iinds - 1) * 3 +1 + ioffs]
    alt_nuc <- c(rep('A', length(ind_a_2)),
                 rep('A', length(ind_a_3)),
                 rep('G', length(ind_g_2)),
                 rep('G', length(ind_g_3)),
                 rep('T', length(ind_t_1)))
    cmp <- c("T","G","C","A")
    names(cmp) <- c("A","C","G","T")
    if(all(GenomicRanges::strand(exn) == '-')){
        alt_nuc <- cmp[alt_nuc]
        ref_nuc <- cmp[ref_nuc]
        names(alt_nuc) <- NULL
        names(ref_nuc) <- NULL
    }

    #- starts in genome coordinates
    #------------------------------
    #- we map each nucleotide in bp_map
    se <- GenomicRanges::start(exn)
    ee <- GenomicRanges::end(exn)
    if(all(GenomicRanges::strand(exn)=='+')) bp_map <- lapply(seq_len(length(se)),function(ind) seq(se[ind],ee[ind]))
    if(all(GenomicRanges::strand(exn)=='-')) bp_map <- lapply(seq_len(length(se)),function(ind) seq(se[ind],ee[ind]) |> rev())

    #- collect them in a matrix (each codon is one column)
    bp_map               <- bp_map |> unlist()
    bp_map_cdn           <- matrix(bp_map, nrow=3)
    colnames(bp_map_cdn) <- cdn_vec

    #- pull out the coordinates of the stop-making nucleotide
    starts <- bp_map_cdn[cbind(ioffs+1, iinds)]

    #- remove start in the last (stop) codon, if applicable
    last_cdn_inds <- bp_map_cdn[,ncol(bp_map_cdn)]
    inds_out      <- which(starts %in% last_cdn_inds) |> sort()

    #- make a key; we zero-pad so that we *could* search better in a trie (probably won't)
    starts <- stringr::str_pad(starts,9L,pad="0")
    kys <- paste0(as.character(seqnames(exn))[1], ":", starts, "|", ref_nuc, "|", alt_nuc)
    if(length(inds_out)>0) kys <- kys[-inds_out]

    return(kys)
}

#- genome-wide stop-making SNVs
#------------------------------
our_tx_ids <- future::value(._EA_txs_grl)$tx_id
tictoc::tic()
tmp <- pbapply::pblapply(our_tx_ids, function(x) get_stop_making_snvs(x))
tictoc::toc()

#- What do we find? ~3.5 million ptc-making snvs.
table(lapply(tmp, is.null) |> unlist())                              #- 31,302 transcripts with stp-gains
kys <- tmp[! (lapply(tmp,is.null) |> unlist()) ] |> unlist() |> unique() #- 3,641,384 stop-gain-making

#- check: does the reference in the key match the reference in the genome?
#-------------------------------------------------------------------------
set.seed(82522)
kys_chk <- kys[sample(seq_len(length(kys)),1000)]
chr_chk <- strsplit(kys_chk,split=":")
chr_chk <- lapply(chr_chk, "[[", 1) |> unlist()
sta_chk <- strsplit(kys_chk,split="\\|")
ref_chk <- lapply(sta_chk, "[[", 2) |> unlist()
sta_chk <- lapply(sta_chk, "[[", 1) |> unlist()
sta_chk <- strsplit(sta_chk,split=":")
sta_chk <- lapply(sta_chk, "[[", 2) |> unlist()

seq <- getSeq(Hsapiens, GRanges(paste0(chr_chk,":", sta_chk)))
if(!all(as.character(seq) == ref_chk)) error("reference allele does not match genome.")

#- save them in a trie to look up; (works faster than hash, did a benchmark)

m_keys <- kys
m_vals <- rep(TRUE, length(m_keys)) #- dummy value - could use rules in second step.
m_trie <- triebeard::trie(keys = m_keys, values = m_vals)

#- since triebeard does not have serialize, we'll do it by hand:
#  we save:
saveRDS(m_keys, file="./inst/extdata/tri-keys_ensdb_v105_fil_all-stop-making-snvs.rds")
saveRDS(m_vals, file="./inst/extdata/tri-vals_ensdb_v105_fil_all-stop-making-snvs.rds")
#  and read in by hand:
m_keys <- readRDS(file="./inst/extdata/tri-keys_ensdb_v105_fil_all-stop-making-snvs.rds")
m_vals <- readRDS(file="./inst/extdata/tri-vals_ensdb_v105_fil_all-stop-making-snvs.rds")
m_trie <- triebeard::trie(keys = m_keys, values = m_vals)

