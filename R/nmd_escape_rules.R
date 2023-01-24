#' Determine NMD escape rules
#'
#' @param ptc_loc Integer. Location of PTC. In CDS codon coordinates (i.e., amino acids), 1-based.
#' @param exn_ind Integer. Index of the exon the variant overlaps, 1-based.
#' @param exn_sta Integer. Start of the exon. In TX coordinates (i.e., nucleotides), 1-based, inclusive (i.e., first nucleotide of the exon).
#' @param exn_end Integer. End of the exon. In TX coordinates, 1-based, inclusive (i.e., last nucleotide of the exon).
#' @param num_exn Integer. Number of exons in the CDS of the transcript (i.e., coding exons).
#' @param txname  Character. Ensembl transcript id.
#' @param css_prox_dist Integer. Distance cutoff (strict, <) for NMD-escape region near the coding start site. Default: 150.
#' @param penultimate_prox_dist Integer. Distance cutoff (strict, <) for NMD-escape region neat the 3'-boundary of the penultimate exon. Default: 50.
#' @return Logical. Vogical vector with six named elements.
#' @details 
#' Return value (Logical vector) has six named elements.
#' - \code{is_ptc = TRUE} if there is a PTC in the alternative transcript. 
#' - \code{is_last = TRUE} if the PTC is in the last (3'-most) coding exon.
#' - \code{is_penultimate3PP = TRUE} if the PTC is within the last (i.e., downstream, 3'-most) \code{penultimate_prox_dist}bp of the penultimate exon.
#' - \code{is_css5PP  = TRUE} if the PTC is within the first (upstream, 5'-most) \code{css_prox_dist} of the translation/coding start site. 
#' - \code{is_single = TRUE} if the PTC is in a single exon transcript.
#' - \code{is_407plus = TRUE} if the PTC is in an exon that is longer than 407bp.

apply_nmd_escape_rules <- function(ptc_loc, exn_ind, exn_sta, exn_end, num_exn, txname, css_prox_dist = 150L, penultimate_prox_dist = 50L){
    #----------------------------------------------------------------------

    res <- rep(FALSE, 6)
    names(res) <- c("is_ptc","is_last", "is_penultimate", "is_cssProximal", "is_single", "is_407plus")

    #- didnt' get a PTC
    #------------------
    if(is.na(ptc_loc) || is.null(ptc_loc)) return(res)
    res["is_ptc"] <- TRUE

    #- are we in the first 150 bp proximal to the css?
    if( get_dist_upstream(ptc_loc, 1L) < css_prox_dist)  res["is_cssProximal"] <- TRUE

    #- are we in the last (coding!) exon?
    if(exn_ind == num_exn) res["is_last"] <- TRUE

    #- are we in a single exon transcrtipt
    #- need to use ._EA_set_env b/c num_exn only counts the coding exons, here we need all
    if( exists(txname, where = future::value(._EA_set_env)) ) res["is_single"] <- TRUE

    #- the exon longer than 407bp?
    if( (exn_end - exn_sta) >= 407L ) res["is_407plus"] <- TRUE

    #- are we in (the last 51bp of) the penultimate exon?
    cnd1 <- exn_ind == (num_exn-1)
    cnd2 <- get_dist_downstream(ptc_loc, exn_end) < penultimate_prox_dist
    if( cnd1 && cnd2) res["is_penultimate"] <- TRUE

    return(res)
}
