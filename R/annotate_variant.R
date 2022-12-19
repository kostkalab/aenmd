
collect_vars_by_exons <- function(txname, vars){
#===============================================

    #- get exons for the transcript
    if( is.null( exn <- get0(txname, future::value(._EA_exn_env)) )){
        stop(paste0("Cannot find exons for ", txname))
    }

    #- find exon-variant combinations
    #- up till now: only one exon per variant
    #. i.e., multiple exon spanning variants are discarded.
    ov <- GenomicRanges::findOverlaps(vars,exn) #- query = vars, subjects = exn
    qH <- S4Vectors::queryHits(ov)
    sH <- S4Vectors::subjectHits(ov)
    if( (qH |> max() |> table()) > 1) stop("multiple-exon variants are not supported.")
    #- we can now look up for each variant its exon;
    #sHm        <- sH
    #names(sHm) <- qH
    sHm <- tapply(sH, qH, \(x) x, simplify = FALSE)
    return(list(exn = exn, exn_x_vrs = sHm))
}

#- exn_x_vrs has exons indexes as values and variant indexes as names. 
annotate_vars_by_tx_snv <- function(txname, vars, exn, exn_x_vrs, css_prox_dist = 150L,
                                    penultimate_prox_dist = 50L, detailed = FALSE){
#===================================================================================

    #-exon starts and ends (genome, "transcript", protein) for reference
    #- note: we could pre-caclulate that if it turns out slow, but unlikely.
    exn_sta_g <- GenomicRanges::start(exn)
    exn_end_g <- GenomicRanges::end(exn)
    exn_sta_t <- c(1, cumsum(GenomicRanges::width(exn))[-length(exn)]+1)
    exn_end_t <- cumsum(GenomicRanges::width(exn))
    exn_sta_p <- (exn_sta_t -1) %/% 3 + 1
    exn_end_p <- (exn_end_t -1) %/% 3 + 1
    #- which nucleotide in codon: 1, 2 or 3
    exn_sta_p_nc <- (exn_sta_t -1) %% 3 + 1
    exn_end_p_nc <- (exn_end_t -1) %% 3 + 1

    #------
    #- All SNVs are PTC-generating (others have been filtered out)
    #- Exon boundaries don't change.
    evr_ind_snv <- (vars$type == "snv") |> which()
    if( length(evr_ind_snv) == 0 ){ #- already done here
        return(NULL)
    }
    #- exons for SNVs
    exn_ind_snv <- exn_x_vrs[evr_ind_snv |> as.character()]
    #- SNVS should not span multiple exons
    if( (lapply(exn_ind_snv, length) |> unlist() > 1) |> any()){
        stop("SNV spans multiple exons. Aborting.")
    }
    exn_ind_snv <- exn_ind_snv |> unlist() #- use it as vector here.

    #- get the codons
    if(all(GenomicRanges::strand(exn)=="-")){
        #- index of codon in "current" exon
        cdn_ind_curr <- ( exn_end_g[exn_ind_snv] - GenomicRanges::start(vars[evr_ind_snv]) +
                              exn_sta_p_nc[exn_ind_snv] -1 ) %/% 3 + 1
        cdn_ind      <- cdn_ind_curr + exn_sta_p[exn_ind_snv] - 1
    } else{
        cdn_ind_curr <- ( GenomicRanges::start(vars[evr_ind_snv]) - exn_sta_g[exn_ind_snv] +
                              exn_sta_p_nc[exn_ind_snv] -1 ) %/% 3 + 1
        cdn_ind      <- cdn_ind_curr + exn_sta_p[exn_ind_snv] - 1
    }

    tbl_snv <- purrr::pmap_dfr( list( ptc_loc = cdn_ind,
                                      exn_ind = exn_ind_snv,
                                      exn_sta = exn_sta_t[exn_ind_snv],
                                      exn_end = exn_end_t[exn_ind_snv],
                                      num_exn = length(exn_sta_t),
				                      txname  = txname,
                                      css_prox_dist = css_prox_dist,
                                      penultimate_prox_dist = penultimate_prox_dist),
                                apply_nmd_escape_rules)

    if(detailed){
        tbl_snv <- tbl_snv |> dplyr::mutate(#exon_index = ,exn_ind_snv
            ptc_position_alt_p = cdn_ind,
            tx_id = exn$tx_id[exn_ind_snv],
            exn_id = exn$exon_id[exn_ind_snv])
    }

    return(list(vars_snv = vars[evr_ind_snv], tbl_snv = tbl_snv))

}


annotate_vars_by_tx_idl <- function(txname, vars, exn, exn_x_vrs, css_prox_dist = 150L, 
                                    penultimate_prox_dist = 50L, detailed = FALSE){
#===================================================================================

    #-exon starts and ends (genome, "transcript", protein) for reference
    #- note: we could pre-caclulate that if it turns out slow, but unlikely.
    exn_sta_g <- GenomicRanges::start(exn)
    exn_end_g <- GenomicRanges::end(exn)
    exn_sta_t <- c(1, cumsum(GenomicRanges::width(exn))[-length(exn)]+1)
    exn_end_t <- cumsum(GenomicRanges::width(exn))
   # exn_sta_p <- (exn_sta_t -1) %/% 3 + 1
   # exn_end_p <- (exn_end_t -1) %/% 3 + 1
    #- which nucleotide in codon: 1, 2 or 3
   # exn_sta_p_nc <- (exn_sta_t -1) %% 3 + 1
   #exn_end_p_nc <- (exn_end_t -1) %% 3 + 1

    #- Here we need to make the new protein to find the PTC position
    #  (that's only strictly true for insertions, but okay)
    evr_ind_idl <- (vars$type != 'snv') |> which()
    if(length(evr_ind_idl) == 0){	
    	return(NULL)#- we are done here.
    }
    exn_ind_idl     <- exn_x_vrs[evr_ind_idl |> as.character()] #- LIST; possible multi-exon variants.
    exn_ind_idl_5p  <- exn_ind_idl |> lapply(min) |> unlist()   #- VECTOR; holds 5p-most exons for multi-exon variants

    #- get the alternative version of the DNA sequence for each variant
    #------------------------------------------------------------------
    if( is.null( seq_ref <- get0(txname, future::value(._EA_cds_env)) )){
        stop(paste0("Cannot find sequence for ", txname))
    }

    #- need to map the genomic variants to the reference protein (CDS/nuc) coordinates
    if(all(GenomicRanges::strand(exn)=="-")){
        ref_nuc_sta <- exn_end_g[exn_ind_idl_5p] - GenomicRanges::end(vars[evr_ind_idl])     + exn_sta_t[exn_ind_idl_5p]
        ref_nuc_end <- exn_end_g[exn_ind_idl_5p] - GenomicRanges::start(vars[evr_ind_idl])   + exn_sta_t[exn_ind_idl_5p]

    } else{
        ref_nuc_sta <- GenomicRanges::start(vars[evr_ind_idl]) - exn_sta_g[exn_ind_idl_5p] + exn_sta_t[exn_ind_idl_5p]
        ref_nuc_end <- GenomicRanges::end(vars[evr_ind_idl])   - exn_sta_g[exn_ind_idl_5p] + exn_sta_t[exn_ind_idl_5p]
    }

    #--------------------------------------------------------------------------------------
    #- deal with OVERHANGING variants, i.e. variants that extend beyond the coding sequence
    #--------------------------------------------------------------------------------------
    #- example:
    #
    #     XX|X                   Variant
    #  =====|ATG=...==TER|=====  Transcript / CDS
    #     AA|ATG                 Reference AAA
    #     A-|-TG                 Alternative A
    #
    #                   X|XXX    Variant 
    #  =====|ATG=...==TER|=====  Transcript / CDS
    #                 TAG|ATT    Reference GATT
    #                   T|A      Alternative TA
    #
    #  HOW WE DEAL WITH IT:
    #
    #  1. Append the reference sequence to include overhaning sequence
    #  2. Adjust ref_nuc_sta and ref_nuc_end accordingly.
    #  (note): there is only one ref_seq but potentially many variants.
    #          so we do this with the maximum "overhangs" of all the variants.
    #-----------------------------------------------------------------------------------------

    #- find out variants that need overhanging treatment
    over_5p_lgl <- ( ( ref_nuc_sta < min(exn_sta_t) ) & (ref_nuc_end >= min(exn_sta_t)) )
    over_5p_ind <- over_5p_lgl |> which()

    over_3p_lgl <-  ( ( ref_nuc_end > max(exn_end_t) ) & (ref_nuc_sta <= max(exn_end_t)) )
    over_3p_ind <-  over_3p_lgl |> which()

    #- find maximum overhangs
    over_5p <- ( (-1)*(ref_nuc_sta - 1) )      |> pmax(0)
    #over_3p <- ref_nuc_end |> (\(x) max(x) - max(exn_end_t))() |> (\(x) pmax(x,0))()
    over_3p <- ( ref_nuc_end - max(exn_end_t) ) |> pmax(0)


    max_over_5p_val <- max(over_5p)
    if(max_over_5p_val > 0){
        max_over_5p_ind <- which.max(over_5p) #- it is okay if there are more than one, can take any. (reference is the same)
    } else {
        max_over_5p_ind <- vector(length = 0, mode='integer')
    }
    max_over_3p_val <- max(over_3p)
    if(max_over_3p_val > 0){
        max_over_3p_ind <- which.max(over_3p) #- gain, if there are many, any of them will do. (reference is the same)
    } else {
        max_over_3p_ind <- vector(length = 0, mode='integer')
    }
    #- Append sequence on the 5' and fix ref_nuc_sta (to be positive), and the exon starts
    if(max_over_5p_val > 0){
        over_seq    <- vars[evr_ind_idl][max_over_5p_ind]$ref |> Biostrings::subseq(1L, max_over_5p_val)
        if( (GenomicRanges::strand(exn) == "-") |> all()) over_seq <- over_seq |> Biostrings::reverseComplement()
        seq_ref     <- Biostrings::xscat(over_seq,seq_ref)[[1]]
        ref_nuc_sta <- ref_nuc_sta + max_over_5p_val 
        ref_nuc_end <- ref_nuc_end + max_over_5p_val 
        exn_sta_t   <- exn_sta_t + max_over_5p_val
        exn_end_t   <- exn_end_t + max_over_5p_val
    }
    #- Append sequence on the 3' end, ref_nuc_end us already correct
    if(max_over_3p_val > 0){
        refwidth    <- vars[evr_ind_idl][max_over_3p_ind]$ref |> Biostrings::width()
        over_seq    <- vars[evr_ind_idl][max_over_3p_ind]$ref |> Biostrings::subseq(refwidth - max_over_3p_val +1, refwidth)
        if( (GenomicRanges::strand(exn) == "-") |> all()) over_seq <- over_seq |> Biostrings::reverseComplement()
        seq_ref     <- Biostrings::xscat(seq_ref, over_seq)[[1]]
    }


    #- GET alternative allele sequence
    #---------------------------------
    alt_vrs <- vars$alt[evr_ind_idl]
    if(all(GenomicRanges::strand(exn)=="-")) alt_vrs <- Biostrings::reverseComplement(alt_vrs)
    make_alt <- function(ref_nuc_sta, ref_nuc_end, alt, seq_ref){
        if(ref_nuc_sta < length(seq_ref)) {
            Biostrings::xscat(Biostrings::subseq(seq_ref, 1, ref_nuc_sta -1), alt,
                          Biostrings::subseq(seq_ref, ref_nuc_end +1, length(seq_ref)))[[1]]
        } else if (ref_nuc_sta == length(seq_ref)){
            Biostrings::xscat(Biostrings::subseq(seq_ref, 1, ref_nuc_sta -1), alt)[[1]]

        } else {
            stop('Variant starts out of (extended) CDS')
        }
    }

    seq_alt <- sapply(seq_len(length(evr_ind_idl)),
                      function(ind) make_alt(ref_nuc_sta = ref_nuc_sta[ind],
                                             ref_nuc_end[ind],
                                             alt_vrs[ind],
                                             seq_ref) )

    #- Divvy out variants overlapping the start codon and treat
    #----------------------------------------------------------
    over_1_lgl    <- (ref_nuc_sta <= 3 + max_over_5p_val) #- those overlap the start codon
    over_1_ind    <- which(over_1_lgl) |> sort()
    n_over_1_lgl  <- !over_1_lgl
    n_over_1_ind  <- which(n_over_1_lgl) |> sort()
    num_over_1    <- sum(over_1_lgl)
    num_n_over_1  <- sum(n_over_1_lgl)

    if(num_n_over_1>0 ){
 	    #- translate right away
   	    seq_alt_p_n_over_1 <- seq_alt[n_over_1_ind] |> Biostrings::DNAStringSet() |> Biostrings::translate() |> suppressWarnings()
        seq_alt_n_over_1   <- seq_alt[n_over_1_ind]
    } 
    if(num_over_1 > 0 ){
   	    #- find the start codons (if there are any)
	    tmp_atg <- Biostrings::vmatchPattern(pattern='atg', subject=Biostrings::DNAStringSet(seq_alt[over_1_ind])) |> Biostrings::start() 
	    tmp_ttg <- Biostrings::vmatchPattern(pattern='ttg', subject=Biostrings::DNAStringSet(seq_alt[over_1_ind])) |> Biostrings::start() 
	    tmp_ctg <- Biostrings::vmatchPattern(pattern='ctg', subject=Biostrings::DNAStringSet(seq_alt[over_1_ind])) |> Biostrings::start() 

	    #- get the first one of each type
	    tmp_atg <- tmp_atg |> lapply(min) |> suppressWarnings()
	    tmp_ttg <- tmp_ttg |> lapply(min) |> suppressWarnings()
	    tmp_ctg <- tmp_ctg |> lapply(min) |> suppressWarnings()

	    #- get the first one for each sequence (sci = start codon index)
	    sci <- sapply(seq_len(tmp_atg |> length()), function(ind) min(c(tmp_atg[[ind]], tmp_ttg[[ind]], tmp_ctg[[ind]])))
	 
	    #- do we have sequences we need to drop? such a pain...
	    inds_in <- (! is.infinite(sci)) |> which()
	    if((inds_in |> length()) == 0){
            seq_alt_over_1 = NULL
		    seq_alt_p_over_1 = NULL
	    } else {
		    sci <- sci[inds_in]
            seq_alt_over_1   <- seq_alt[over_1_ind][inds_in] |> Biostrings::DNAStringSet()  |>
                                    Biostrings::subseq(sci)                                 |>
                                    suppressWarnings()
		    seq_alt_p_over_1 <- seq_alt_over_1 |> Biostrings::translate()                   |> 
								    suppressWarnings()
	     }
    }

    #- seq_alt_p: conatenate and restore order, 
    #- drop sequences without start codon 
    ind_in <- seq_len(evr_ind_idl |> length())
    if( (num_over_1 > 0) && ( num_n_over_1 > 0) ){ #- both types
	    #-concatenating the NULL works
  	    seq_alt_p <- c(seq_alt_p_n_over_1, seq_alt_p_over_1)
        seq_alt   <- c(seq_alt_n_over_1, seq_alt_over_1)
	    ind_in    <- c(n_over_1_ind, over_1_ind[inds_in])
        sci       <- c(rep(1+max_over_5p_val, num_n_over_1),sci) #- need to add the 5'overhang to sci here!
    } else if(num_over_1 > 0){ #- only start-overlapping
	    seq_alt_p <- seq_alt_p_over_1
        seq_alt   <- seq_alt_over_1
    	ind_in    <- over_1_ind[inds_in]
        sci       <- sci
    } else if(num_n_over_1 > 0){ #- no start-overlapping
	    seq_alt_p <- seq_alt_p_n_over_1
        seq_alt   <- seq_alt_n_over_1 
	    ind_in    <- n_over_1_ind
        sci       <-  rep(1, num_n_over_1)
    }

    
    #- restore order
    ord       <- order(ind_in)
    ind_in    <- ind_in[ord]
    seq_alt   <- seq_alt[ord]
    seq_alt_p <- seq_alt_p[ord]
    sci       <- sci[ord]

    #-subset all relevant quantities to reflect omission of non-start alternatives
    evr_ind_idl <- evr_ind_idl[ind_in]
    exn_ind_idl <- exn_ind_idl[ind_in] #- can still be a LIST

    #- now we can look for stop codons
    #---------------------------------
    stop_pos  <- Biostrings::vmatchPattern("*",seq_alt_p)
    stop_pos  <- lapply(stop_pos, function(x) {
        if(length(x) == 0) {
            fst_stop <- NA
        } else {
            fst_stop <- x[1] |> Biostrings::start() |> as.integer()
        }
    }) |> unlist()

    #- ptc does not contain "original" stop codon at the end of the seuence
    #  the corner cases here that we miss contain indels affecting the
    #  last reference codon and wrongly throw them out...
    ptc_pos <- stop_pos
    ptc_pos[stop_pos == Biostrings::width(seq_alt_p)] <- NA

    #- update exon boundaries and apply rules
    #----------------------------------------
    #- now we need the updated exon boundaries in CDS/codon coordinates and TX coordinates.
    #  since all variants are contained in one exon, we just need to update that.
    #  dw = delta width

    d_w <- vars[evr_ind_idl]$alt |> Biostrings::width() -
        vars[evr_ind_idl]$ref |> Biostrings::width()

    #- exn_inds: exons (in the reference! that get ovdrlapped by variant.)
    afu <- function(txname, exn_sta_t, exn_end_t, exn_inds, num_exn, d_w, sci, ptc_pos, alt_seq,
                    css_prox_dist, penultimate_prox_dist){
        #-----------------------------------------------
        #- if ptc_pos is NA we return FALSE
        if(is.na(ptc_pos)){
            res <- rep(FALSE, 6)
            names(res) <- c("is_ptc","is_last", "is_penultimate", "is_cssProximal", "is_single", "is_407plus")
            return(res)
        }
       
        #argg <- c(as.list(environment()))
        #print(argg)
        #- Get exon boundaries in alternative sequence (note, we start at the new CSS)
        tmp <- mev_alt_exn_bnd(exn_sta_t, exn_end_t, exn_inds, d_w, sci)
        exn_sta_t_alt <- tmp$exn_sta_t_alt
        exn_end_t_alt <- tmp$exn_end_t_alt

        #- find the exon that contains the PTC (not necessarily the variant)
        #- note, if the PTC overlaps multiple exons, we use the 3'-most.
        #- note, since the reference sequence can have 5' and 3' overhang, we need to 
        #  check that the ptc is contained in the alternative exons of this variant.
        ptc_coords_alt_t <- (ptc_pos-1)*3 + 1:3
        cnd1 <- ptc_coords_alt_t[3] <= max(exn_end_t_alt)
        cnd2 <- ptc_coords_alt_t[1] >= min(exn_sta_t_alt)
        if( cnd1 && cnd2 ){
            #- have a "real" PTC
            exn_ind_ptc <- (exn_end_t_alt >= ptc_coords_alt_t[3]) |> which() |> min()
        } else {
            #- e.g., PTC was found in sequence related to another variant that is not present in this variant
            ptc_pos <- NA
            res <- rep(FALSE, 6)
            names(res) <- c("is_ptc","is_last", "is_penultimate", "is_cssProximal", "is_single", "is_407plus")
            return(res)
        }
        #- location of ptc in that exon
        res <- apply_nmd_escape_rules(ptc_pos, exn_ind_ptc, exn_sta_t_alt[exn_ind_ptc],
                                exn_end_t_alt[exn_ind_ptc], num_exn, txname,
                                css_prox_dist = css_prox_dist, penultimate_prox_dist = penultimate_prox_dist)
        #- FIXME: add these 
        if(FALSE){
            res                 <- as.list(res)
            res$seq_alt_p       <- alt_seq_p[1:ptc_pos] #- alternative coding sequence 
            res$seq_alt_p_start <- sci #- position of new CSS in alternative sequence
            res$seq_alt_p_end   <- ptc_pos #- AA position of TER in alternative sequence, starting at (potentially new) CSS.
            res$exn_sta_alt     <- exn_sta_t_alt #- with "alternative" CSS
            res$exn_end_alt     <- exn_end_t_alt #- with "alternative" CSS
        }

        return(res)
    }

    tmp <- lapply(seq_len( vars[evr_ind_idl] |> length()), function(var_ind) afu(
                                        txname    = txname,
                                        exn_sta_t = exn_sta_t,
                                        exn_end_t = exn_end_t,
                                        exn_inds  = exn_ind_idl[[var_ind]], #- can a vector
                                        num_exn   = length(exn_sta_t),
                                        d_w       = d_w[var_ind],
                                        sci       = sci[var_ind],
                                        ptc_pos   = ptc_pos[var_ind],
                                        alt_seq   = seq_alt_p[[var_ind]],
                                        css_prox_dist = css_prox_dist,
                                        penultimate_prox_dist = penultimate_prox_dist))

    #conversion; FIXME: does not work for tmp where tmp is a list of lists that cannnot be coerced to a matrix...
    mat <- tmp |> unlist() |> matrix(ncol = 6, byrow = TRUE) 
    colnames(mat) <- names(tmp[[1]])

    tbl_idl <- tibble::as_tibble(mat)

    if(detailed){
        tbl_idl <- tbl_idl |> dplyr::mutate(
            ptc_position_alt_p = stop_pos,
            tx_id = exn$tx_id[exn_ind_idl],
            exn_id = exn$exon_id[exn_ind_idl]) #- this is not the ptc exon but the variant exon
    }
    return(list(vars_idl = vars[evr_ind_idl], tbl_idl = tbl_idl))
}



#' Annotate variants overlapping a transcript
#'
#' @param txname Character. Transcript name.
#' @param vars GRanges. Overlapping variants.
#' @param css_prox_dist Integer. Distance to the CSS defining NMD escape regions.
#' @param penultimate_prox_dist Integer. Distance to the penultimate exon 3'end defining NMD escape regions.
#' @param detailed Logical. Should additional information be returned.
#' @return GRanges.
annotate_variants_by_tx <- function( txname, vars, css_prox_dist = 150L, 
                                     penultimate_prox_dist = 50L, detailed = FALSE){
#===================================================================================

    #- get vars overlapping the tx, and exons
    #----------------------------------------
    res_tmp <- collect_vars_by_exons(txname, vars)

    #- get the SNVs
    #--------------
    res_snv <- annotate_vars_by_tx_snv(txname, vars, res_tmp$exn, res_tmp$exn_x_vrs, 
                                       css_prox_dist = css_prox_dist, penultimate_prox_dist = penultimate_prox_dist, detailed = detailed)

    #- get the indels
    #----------------
    res_idl <- annotate_vars_by_tx_idl(txname, vars, res_tmp$exn, res_tmp$exn_x_vrs,
                                       css_prox_dist = css_prox_dist, penultimate_prox_dist = penultimate_prox_dist, detailed = detailed)

    #- collate and return
    if(is.null(res_snv)){
        if(is.null(res_idl)){
            #- NO variants, return empty range
            empty_range <- GenomicRanges::GRanges('chr0',IRanges::IRanges(start=0,width=0))
            return(empty_range)
        } else {
            #- return INDELs only
            res <- res_idl$vars_idl
            dfr <- res_idl$tbl_idl |> S4Vectors::DataFrame()
            S4Vectors::mcols(res)$res_aenmd <- dfr
            return(res)
        }
    } else if(is.null(res_idl)){
        #- return SNVs only
        res <- res_snv$vars_snv
        dfr <- res_snv$tbl_snv |> S4Vectors::DataFrame()
        S4Vectors::mcols(res)$res_aenmd <- dfr
        return(res)
    } else {
        #- return SNVs and INDELs
        res <- c(res_snv$vars_snv, res_idl$vars_idl)
        dfr <- dplyr::bind_rows(res_snv$tbl_snv, res_idl$tbl_idl) |> S4Vectors::DataFrame()
        S4Vectors::mcols(res)$res_aenmd <- dfr
        return(res)
    }
}


