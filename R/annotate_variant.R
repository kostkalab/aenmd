#' Workhorse function annotating a variant with escape from nonsense mediated decay
#' @param rng GRanges object. Start is first base of variant, end is last base of variant.
#'            Needs an ``ref`` and ``alt`` metadata column.
#' @param tx String. The transcript ID for which the variant is annotated.
#' @return List.
#' @importFrom utils head tail
#' @examples
#' @export
annotate_variant <- function(rng, tx){
#=====================================

    #- refrence exons and reference sequence from their environments
    #- envs are in ./inst/extdata and loaded in zzz.R.
    exn     <- ._EA_exn_env[[tx]] #- exons of the transcript;
    seq_ref <- ._EA_cds_env[[tx]] #- reference allele sequence;

    #- add the exon starts to seq_ref
    S4Vectors::metadata(seq_ref) <- list(exon_starts = c(1, cumsum(GenomicRanges::width(exn))+1) |>
                                                            head(n=-1))

    #- get first and last base of variant
    sae <- get_start_end(gr1 = rng, gr2 = exn)

    #- get the alternative version of the DNA sequence
    alt_var <- rng$alt
    if(all(GenomicRanges::strand(exn) == '-')) {
        alt_var <- Biostrings::reverseComplement(alt_var)
    }
    seq_alt <- Biostrings::xscat(Biostrings::subseq(seq_ref, 1, sae[1]-1), alt_var,
                                 Biostrings::subseq(seq_ref, sae[2]+1, length(seq_ref)))[[1]]

    #- map exon starts to alternative sequence
    es_new      <- S4Vectors::metadata(seq_ref)$exon_starts
    delta_len   <- Biostrings::width(rng$alt) - Biostrings::width(rng$ref)
    ind         <- es_new > sae[2]
    es_new[ind] <- es_new[ind] + delta_len
    #- we could have removed whole exons:
    drop                         <- es_new <= max(es_new[!ind])
    drop                         <- drop & ind
    es_new                       <- es_new[!drop]
    S4Vectors::metadata(seq_alt) <- list(exon_starts = es_new)

    #- translate into protein space
    seq_ref_p <- Biostrings::translate(seq_ref)
    seq_alt_p <- Biostrings::translate(seq_alt)
    #- add the exon boundaries (1st aa in each exon) to protein space
    S4Vectors::metadata(seq_ref_p)$exon_starts <- (S4Vectors::metadata(seq_ref)$exon_starts -1) %/% 3 + 1
    S4Vectors::metadata(seq_alt_p)$exon_starts <- (S4Vectors::metadata(seq_alt)$exon_starts -1) %/% 3 + 1

    #- get all the stop codons in the alternative sequence
    stop_pos <- Biostrings::matchPattern("*",seq_alt_p)
    if(length(stop_pos) == 0) {
        fst_stop <- NA
    } else {
        fst_stop <- stop_pos[1] |> Biostrings::start()
    }

    #- use the first stop codon to apply NMD-escape rules
    #----------------------------------------------------

    #- do we have a PTC?
    #- IF the variant does NOT overlap the orginal stop
    #  AND the first stop codon is at the end, we are not having a PTC.
    is_ptc <- FALSE
    if(!is.na(fst_stop)){

        #- we have a PTC if it is not at the end
        cnd1 <- fst_stop < length(seq_alt_p)
        #- we have a PTC if it is at the end, but the original one was (partially) deleted/expanded.
        cnd2 <- (delta_len>0) & (sae[1] <= (length(seq_ref) -2)) & (sae[2] >= length(seq_ref)-2)

        if(cnd1 | cnd2){
            is_ptc <- TRUE
        }
    }

    #- last exon rule (TRUE if first PTC overlaps last exon)
    is_last <- FALSE
    if(!is.na(fst_stop) & (length(S4Vectors::metadata(seq_alt_p)$exon_starts) > 1)){
        cnd1 <- fst_stop >=  ( S4Vectors::metadata(seq_alt_p)$exon_starts |> tail(n=1) )
        if(cnd1){
            is_last <- TRUE
        }
    }

    #- penultimate exon rule (true if PTS overlaps last 50bp of penultimate exon)
    #- note that we use 17 AA, which is 51 bp and not 50
    is_penultimate <- FALSE
    if(!is.na(fst_stop) & (length(S4Vectors::metadata(seq_alt_p)$exon_starts) > 1)){

        cnd1 <- fst_stop >=  (S4Vectors::metadata(seq_alt_p)$exon_starts |> tail(n=2))[1]
        cnd2 <- fst_stop <   (S4Vectors::metadata(seq_alt_p)$exon_starts |> tail(n=2))[2]
        cnd3 <- (S4Vectors::metadata(seq_alt_p)$exon_starts |> tail(n=2))[2] - fst_stop <= 17

        if(cnd1 & cnd2 & cnd3){
            is_penultimate <- TRUE
        }
    }

    #- first exon rule
    is_first <- FALSE
    if(!is.na(fst_stop) & (length(S4Vectors::metadata(seq_alt_p)$exon_starts) > 1)){

        cnd1 <- fst_stop <  (S4Vectors::metadata(seq_alt_p)$exon_starts |> head(n=2))[2]

        if(cnd1){
            is_first <- TRUE
        }
    }

    #- single exon rule
    is_single <- FALSE
    if(length(S4Vectors::metadata(seq_alt_p)$exon_starts) == 1){
        is_single <- TRUE
    }

    rule_res        <- c(is_ptc, is_last, is_penultimate, is_first, is_single)
    names(rule_res) <- c("is_ptc","is_last", "is_penultimate", "is_first", "is_single")

    res <- list( rules   = rule_res,
                 seq_ref = seq_ref_p,
                 seq_alt = seq_alt_p)

    return(res)
}

