#' Get distance between a PTC (or any AA) and a downstream boundary (e.g., exon end)
#' 
#' @param aa_pos Integer. Location of the AA/PTC in the CDS. CDS coordinates (codonn), 1-based.
#' @param ds_bnd Integer. Location of the last nucleotide before the boundary (e.g., last nucleotide in exon), 1-based.
#' @return Integer. Distance of the AA/PTC (acutally, the downstream-most nucleotide covered by the AA/PTC to the boundary).
#' @details
#' \preformatted{
#'NUC:  123456789...  ( '|' = exon-exon boundary (nucleotide to the left is last of upstream exon), X = PTC-covered, O = non-PTC-covered)
#'AA:   111222333
#'PTC:  ATGOOOOOOX|XXOOOOOOOOOOOO  aa_pos: 4, ds_bnd: 10, DIST: 0
#'PTC:  ATGOOOOOOXX|XOOOOOOOOOOOO  aa_pos: 4, ds_bnd: 11, DIST: 0
#'PTC:  ATGOOOOOOXXX|OOOOOOOOOOOO  aa_pos: 4, ds_bnd: 12, DIST: 0
#'PTC:  ATGOOOOOOXXXO|OOOOOOOOOOO  aa_pos: 4, ds_bnd: 13, DIST: 1
#'PTC:  ATGOOOOOOXXXOO|OOOOOOOOOO  aa_pos: 4, ds_bnd: 14, DIST: 2
#'PTC:  ATGOOOOOOXXXOOO|OOOOOOOOO  aa_pos: 4, ds_bnd: 15, DIST: 3
#' }
get_dist_downstream <- function(aa_pos, ds_bnd){
#-------------------------------------------------

    #-nc_max <- (_pos * 3)       #- maximum nucleotide covered by PTC
    # nc_min <- (aa_pos * 3) - 2   #- minimum nucleotide covered by PTC

    #dist <- max(0, (ds_bnd - nc_min - 2))
    #- or, equivalently
    #dist <- max(0, ds_bnd - aa_pos * 3 ) 

    if( ds_bnd < aa_pos * 3L - 2) stop('Downstream boundary before minimum nucleotide covered by AA/PTC')
    return(max(0L, ds_bnd - aa_pos * 3L ))
}

#' Get distance between a PTC (or any AA) and a upstream boundary (e.g., a CSS site)
#' 
#' @param aa_pos Integer. Location of the AA/PTC in the CDS. CDS coordinates (codonn), 1-based.
#' @param us_bnd Integer. Location of the last nucleotide before the boundary (e.g., first nucleotide in exon), 1-based.
#' @return Integer. Distance of the AA/PTC (acutally, the upstream-most nucleotide covered by the AA/PTC to the boundary).
#' @details
#' \preformatted{
#'NUC:  12 3456789...  ( '|' = exon-exon boundary (nucleotide to the left is last of upstream exon), X = PTC-covered, O = non-PTC-covered)
#'PTC  |ATGXXXOOOOOOOOO  aa_pos: 2, us_bnd: 1, DIST: 3
#'PTC  |ATGOOOXXXOOOOOO  aa_pos: 3, us_bnd: 1, DIST: 6
#'PTC  |ATGOOOOOOXXXOOO  aa_pos: 4, us_bnd: 1, DIST: 9
#'PTC  |ATGOOOOOOOOOXXX  aa_pos: 5, us_bnd: 1, DIST: 12
#'
#'NUC:  123456789 
#'PTC:  ATGOOOOOOXX|XOOOOOOOOOOO  aa_pos: 4, us_bnd: 12, DIST: 0
#'PTC:  ATGOOOOOOX|XXOOOOOOOOOOO  aa_pos: 4, us_bnd: 11, DIST: 0
#'PTC:  ATGOOOOOO|XXXOOOOOOXXXOO  aa_pos: 4, us_bnd: 10, DIST: 0
#'PTC:  ATGOOOOO|OXXXOOOOOOXXXOO  aa_pos: 4, us_bnd: 9,  DIST: 1
#'PTC:  ATGOOOO|OOXXXOOOOOOXXXOO  aa_pos: 4, us_bnd: 8,  DIST: 2
#'PTC:  ATGOOO|OOOXXXOOOOOOXXXOO  aa_pos: 4, us_bnd: 7,  DIST: 3
#' }
get_dist_upstream  <- function(aa_pos, us_bnd){
#-------------------------------------------------

    #-nc_max <- (aa_pos * 3)       #- maximum nucleotide covered by PTC
    # nc_min <- (aa_pos * 3) - 2   #- minimum nucleotide covered by PTC

    #dist <- max(0, (nc_min - us_bnd ))
    #- or, equivalently
    #dist <- max(0, aa_pos * 3 - 2 - us_bnd ) 

    if( us_bnd > aa_pos * 3L) stop('Upstream boundary after maximum nucleotide covered by AA/PTC')
    return(max(0L, aa_pos * 3L - us_bnd - 2L ))
}
 


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
