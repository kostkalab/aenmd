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
#' @param aa_pos Integer. Location of the AA/PTC in the CDS. CDS coordinates (codon), 1-based.
#' @param us_bnd Integer. Location of the last nucleotide before the boundary (e.g., first nucleotide in exon), 1-based.
#' @return Integer. Distance of the AA/PTC (acutally, the upstream-most nucleotide covered by the AA/PTC to the boundary).
#' @details
#' \preformatted{
#'NUC:  123456789...  ( '|' = exon-exon boundary (nucleotide to the left is last of upstream exon), X = PTC-covered, O = non-PTC-covered)
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


#' Variant coordinates in CDS space
#'
#' Get first and last base of a variant in a transcript's CDS coordinates
#'
#'      |123|---|4567|-----|89|  <- positive strand
#'      |xxx|---|xxxx|-----|xx|
#'        ===================    gr1 => (2,8)
#'        **     ****       *    intersect(gr1,gr2)
#'
#'       |987|---|6543|-----|21|  <- negative strand
#'       |xxx|---|xxxx|-----|xx|
#'          ===============       gr1 ==> (3,7)
#'          *     ****
#' @param gr1 GRanges object. Single range, no strand. The variant. Start is first nucleotide, end is last nucleotide.
#' @param gr2 GRanges object. Exon ranges. All on the same strand. Sorted 5' to 3' (i.e., "reversed" on negative strand).
#'            Both ``gr1`` and ``gr2`` need to be on the same chromosome.
#' @return Numeric. Vector with start, end, and the overlapping exons.
#' @examples
#' gr1  <- GenomicRanges::GRanges('chr1:2-18:*')
#' gr21 <- GenomicRanges::GRanges('chr1:1-3:+')
#' gr22 <- GenomicRanges::GRanges('chr1:7-10:+')
#' gr23 <- GenomicRanges::GRanges('chr1:16-20:+')
#' gr2  <- c(gr21,gr22,gr23)
#' aenmd:::get_start_end(gr1,gr2) #- 2, 10, 1, 3
get_start_end <- function(gr1, gr2){
#====================================
        if(!( (GenomicRanges::seqnames(gr1) == GenomicRanges::seqnames(gr2)) |> all())){
            stop("chromosome mismatch; not supported.")
        }
        GenomicRanges::strand(gr1) <- '*'
        if(all(as.character(GenomicRanges::strand(gr2)) == "+")){
                strnd = '+'
        } else if(all(as.character(GenomicRanges::strand(gr2)) == "-")){
                strnd = '-'
        } else {
                stop("Consistent meaningful strand required\n")
        }
        n2 <- length(gr2)


        res        <- c(NA,NA)
        names(res) <- c("start","end")

        #gr1     <- GenomicRanges::intersect(gr1, gr2, ignore.strand = TRUE)
        #- for some reason this is a lot faster (order of magnitude)
        #  FIXME: this breaks currently if the intersection is empty.
        #         not an issue b/c we filter for cds overlap before, but not robust.
        gr1     <- GenomicRanges::GRanges(GenomicRanges::seqnames(gr1),
                                          IRanges::intersect( IRanges::IRanges(GenomicRanges::start(gr1),GenomicRanges::end(gr1)),
                                                              IRanges::IRanges(GenomicRanges::start(gr2),GenomicRanges::end(gr2))))
        gr2_hit <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(gr1, gr2))

        #- no overalp, we are done
        if(length(gr2_hit) == 0) return(res)

        ht_up <- min(gr2_hit) #- most 5'
        ht_dn <- max(gr2_hit) #- most 3' #- ht_dn >= ht_up

        up_plus <- 0
        dn_plus <- 0
        if(ht_up > 1)  {
                up_plus <- sum(GenomicRanges::width(gr2)[1:(ht_up-1)])
        }
        if(ht_dn > 1){
            dn_plus <- sum(GenomicRanges::width(gr2)[1:(ht_dn-1)])
        }

        if( strnd == '+'){
                res[1] <- min(GenomicRanges::start(gr1)) - GenomicRanges::start(gr2[ht_up]) + up_plus + 1
                res[2] <- max(GenomicRanges::end(gr1))   - GenomicRanges::start(gr2[ht_dn]) + dn_plus + 1
        } else {
                res[1] <- GenomicRanges::end(gr2[ht_up]) - max(GenomicRanges::end(gr1))   + up_plus + 1
                res[2] <- GenomicRanges::end(gr2[ht_dn]) - min(GenomicRanges::start(gr1)) + dn_plus + 1
        }

        return(c(res[1], res[2], ex_up_idx = ht_up, ex_dn_idx = ht_dn))
}


#' Update exon boundaries for alternative allele's sequence
#' @param exn_sta_t Integer. Vector, the exon starts (in transcript coordinates), reference allele
#' @param exn_end_t Integer. Same as \code{exn_sta_t}, just the exon ends
#' @param exn_inds Integer. The exons (typically one, can be more than one) the variant overlaps
#' @param d_w Integer. The diffrence in length between reference and alterntive allele.
#' @param sci Integer. The start in the alternative sequence where its start codon starts.
#' 
#' @return List. Two components, \code{exn_sta_t_alt} and \code{exn_end_t_alt}.
#' 
mev_alt_exn_bnd <- function(exn_sta_t, exn_end_t, exn_inds, d_w, sci){

        exn_sta_t_alt <- exn_sta_t
        exn_end_t_alt <- exn_end_t

        #- fuse exon boundaries for multi-exon variants
        exn_ind_5p <- min(exn_inds)
        exn_ind_3p <- max(exn_inds)
        if( exn_ind_5p != exn_ind_3p ){
                exn_in  <- exn_ind_5p
                exn_out <- (exn_ind_5p:exn_ind_3p)[-1]
                exn_end_t_alt[exn_in] <- exn_end_t[exn_out |> max()]
                exn_sta_t_alt <- exn_sta_t_alt[-exn_out]
                exn_end_t_alt <- exn_end_t_alt[-exn_out]
                exn_inds <- exn_ind_5p
        }

        #- new exon boundaries (in tx/nuc space)
        #- repair first exon boundary for 5' overhaning variants
        inds          <- exn_ind_5p : length(exn_sta_t_alt)
        if( (min(inds) == 1) && (exn_sta_t[1] != 1)){
                overhang <- exn_sta_t[1] -1 
                exn_sta_t_alt[1] <- 1
                d_w <- sign(d_w) * (abs(d_w) - overhang) #- keep outside CDS bases the same for ref and alt
        }
        delt                    <- rep(d_w, length(inds))
        exn_sta_t_alt[inds[-1]] <- exn_sta_t_alt[inds[-1]] + delt[-1] #- first exon can be wrong, adjust below, so ignore possible overhang here
        exn_end_t_alt[inds]     <- exn_end_t_alt[inds]     + delt

        #- adjust for new CDS, if applicable (first "exon" starts at CSS)
        sta_coords_alt_t  <- (sci - 1) + 1:3
        new_first_exn_ind <- ( exn_sta_t_alt <=  ( sta_coords_alt_t  |> min() ) ) |> which() |> max()
        exn_in            <- new_first_exn_ind:length(exn_sta_t_alt)
        exn_sta_t_alt     <- exn_sta_t_alt[exn_in]
        exn_end_t_alt     <- exn_end_t_alt[exn_in]
        #- delete nucleotides before new start
        del_nucs <- sci-exn_sta_t_alt[1]
        #if( (exn_ind_5p == 1) && (exn_sta_t_alt[1] > 1)){
        #    #- adjust for 5p overhang
        #    del_nucs = del_nucs + exn_sta_t_alt[1] - 1
        #}

        exn_sta_t_alt[1] <- exn_sta_t_alt[1] + del_nucs
        #- make it 1-based, because we are starting with the start codon
        tmp           <- exn_sta_t_alt[1] - 1
        exn_sta_t_alt <- exn_sta_t_alt -  tmp
        exn_end_t_alt <- exn_end_t_alt -  tmp
        
        return(list(exn_sta_t_alt = exn_sta_t_alt, exn_end_t_alt = exn_end_t_alt))
}

is_vcf_rng <- function(vcf_rng, check_key = FALSE){
        pass <- TRUE
        cn   <- S4Vectors::mcols(vcf_rng) |> colnames() #- colnames of the DataFrame

        #- need 'ref' and 'alt' columns
        if( sum(c('ref', 'alt') %in% cn) < 2 ) pass = FALSE

        #- they need to be DNAStringSets
        if(!is(S4Vectors::mcols(vcf_rng)$ref, "DNAStringSet")) pass = FALSE
        if(!is(S4Vectors::mcols(vcf_rng)$alt, "DNAStringSet")) pass = FALSE

        #- check that ranges have the correct length
        if( ! all(GenomicRanges::width(vcf_rng) == Biostrings::width(vcf_rng$ref)) ) pass = FALSE

        #- check that the keys are right (optional)
        if(check_key){
                if(! all(vcf_rng$key == make_keys(vcf_rng))) pass = FALSE
        }

        return(pass)
}

make_keys <- function(vcf_rng){
        starts <- GenomicRanges::start(vcf_rng) |> stringr::str_pad(9L, pad="0")
        keys <- paste0(GenomicRanges::seqnames(vcf_rng), ":", starts,"|" ,vcf_rng$ref, "|", vcf_rng$alt)
        return(keys)
}

filter_splice_overlap <- function(rng){
#--------------------------------------

    ovc <- GenomicRanges::countOverlaps(rng, ad_get_spl_mask(), type="any", ignore.strand = TRUE, minoverlap = 1L)
    ind <- which(ovc>0) |> sort() 
  
    if(length(ind)>0){
        rng <- rng[-ind]
    }
    return(rng)        
}

filter_in_coding_exon <- function(rng){
#--------------------------------------

    #- these are contained in the tx mask, but they might not necessarily be within a single exon
 
    ovc  <- GenomicRanges::countOverlaps(rng, ad_get_txs_mask(), type="within", ignore.strand = TRUE)
    ind  <- which(ovc > 0)
    if(length(ind)>0) rng  <- rng[ind]

    return(rng)

}

filter_ptc_generating_snvs <- function(rng){
#-------------------------------------------

    tps <- rng$type
    kys <- rng$key

    mtc <- ad_is_ptc_snv(kys[tps == 'snv']) 
    ind <- (tps != 'snv')
    ind[tps == 'snv'] <- mtc

    return(rng[ind])
}

stratify_by_tx <- function(rng){
#-------------------------------

    ov  <- GenomicRanges::findOverlaps(rng, ad_get_txs_mask())
    txs <- ad_get_txs_mask()$transcript_id[S4Vectors::subjectHits(ov)] 

    #- invert the list (why is there not a single command?)
    inds        <- rep(seq_len(length(txs)),lengths(txs))
    inds_by_tx  <- tapply(inds, unlist(txs),c)

    return(inds_by_tx)
}




