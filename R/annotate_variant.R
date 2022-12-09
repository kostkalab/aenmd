
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
    #- we can now look up for each variant its exon
    sHm        <- sH
    names(sHm) <- qH

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
    exn_ind_snv <- exn_x_vrs[evr_ind_snv |> as.character()]

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
    exn_sta_p <- (exn_sta_t -1) %/% 3 + 1
    exn_end_p <- (exn_end_t -1) %/% 3 + 1
    #- which nucleotide in codon: 1, 2 or 3
    exn_sta_p_nc <- (exn_sta_t -1) %% 3 + 1
    exn_end_p_nc <- (exn_end_t -1) %% 3 + 1

    #- Here we need to make the new protein to find the PTC position
    #  (that's only strictly true for insertions, but okay)
    evr_ind_idl <- (vars$type != 'snv') |> which()
    if(length(evr_ind_idl) == 0){	
    	return(NULL)#- we are done here.
    }
    exn_ind_idl <- exn_x_vrs[evr_ind_idl |> as.character()]

    #- get the alternative version of the DNA sequence for each variant
    #------------------------------------------------------------------
    if( is.null( seq_ref <- get0(txname, future::value(._EA_cds_env)) )){
        stop(paste0("Cannot find sequence for ", txname))
    }

    #- need to map the genomic variants to the reference protein (CDS/nuc) coordinates
    if(all(GenomicRanges::strand(exn)=="-")){
        ref_nuc_sta <- exn_end_g[exn_ind_idl] - GenomicRanges::end(vars[evr_ind_idl])     + exn_sta_t[exn_ind_idl]
        ref_nuc_end <- exn_end_g[exn_ind_idl] - GenomicRanges::start(vars[evr_ind_idl])   + exn_sta_t[exn_ind_idl]

    } else{
        ref_nuc_sta <- GenomicRanges::start(vars[evr_ind_idl]) - exn_sta_g[exn_ind_idl] + exn_sta_t[exn_ind_idl]
        ref_nuc_end <- GenomicRanges::end(vars[evr_ind_idl])   - exn_sta_g[exn_ind_idl] + exn_sta_t[exn_ind_idl]
    }

    #- cut "ovrhanging" parts of variants
    ohi <- ( ( ref_nuc_end > max(exn_end_t) ) & (ref_nuc_sta <= max(exn_end_t)) ) |> which()
    ref_nuc_end[ohi] <- max(exn_end_t)

    #- Divvy out variants overlapping the start codon
    log_over_1  <- (ref_nuc_sta <= 3) #- those overlap the start codon (but not splice regions)
    ind_over_1  <- which(log_over_1) |> sort()
    log_nover1  <- !log_over_1
    ind_nover1  <- which(log_nover1) |> sort()
    n_nover1    <- sum(log_nover1)
    n_over_1    <- sum(log_over_1)

    alt_vrs <- vars$alt[evr_ind_idl]
    if(all(GenomicRanges::strand(exn)=="-")) alt_vrs <- Biostrings::reverseComplement(alt_vrs)
    make_alt <- function(ref_nuc_sta, ref_nuc_end, alt, seq_ref){
        if(ref_nuc_sta < length(seq_ref)) {
            Biostrings::xscat(Biostrings::subseq(seq_ref, 1, ref_nuc_sta -1), alt,
                          Biostrings::subseq(seq_ref, ref_nuc_end +1, length(seq_ref)))[[1]]
        } else if (ref_nuc_sta == length(seq_ref)){
            Biostrings::xscat(Biostrings::subseq(seq_ref, 1, ref_nuc_sta -1), alt)[[1]]

        } else {
            stop('Variant starts out of CDS')
        }
    }

    seq_alt <- sapply(seq_len(length(evr_ind_idl)),
                      function(ind) make_alt(ref_nuc_sta = ref_nuc_sta[ind],
                                             ref_nuc_end[ind],
                                             alt_vrs[ind],
                                             seq_ref) )

    if(n_nover1 >0 ){
 	    #- translate right away
   	    seq_alt_p_nover1 <- seq_alt[ind_nover1] |> Biostrings::DNAStringSet() |> Biostrings::translate() |> suppressWarnings()
    } 
    if(n_over_1 > 0 ){
   	    #- find the start codons (if there are any)
	    tmp_atg <- Biostrings::vmatchPattern(pattern='atg', subject=Biostrings::DNAStringSet(seq_alt[ind_over_1])) |> Biostrings::start() 
	    tmp_ttg <- Biostrings::vmatchPattern(pattern='ttg', subject=Biostrings::DNAStringSet(seq_alt[ind_over_1])) |> Biostrings::start() 
	    tmp_ctg <- Biostrings::vmatchPattern(pattern='ctg', subject=Biostrings::DNAStringSet(seq_alt[ind_over_1])) |> Biostrings::start() 

	    #- get the first one of each type
	    tmp_atg <- tmp_atg |> lapply(min) |> suppressWarnings()
	    tmp_ttg <- tmp_ttg |> lapply(min) |> suppressWarnings()
	    tmp_ctg <- tmp_ctg |> lapply(min) |> suppressWarnings()

	    #- get the first one for each sequence
	    sci <- sapply(seq_len(tmp_atg |> length()), function(ind) min(c(tmp_atg[[ind]], tmp_ttg[[ind]], tmp_ctg[[ind]])))
	 
	    #- do we have sequences we need to drop? such a pain...
	    inds_in <- (! is.infinite(sci)) |> which()
	    if((inds_in |> length()) == 0){
		    seq_alt_p_over1 = NULL
	    } else {
		    sci <- sci[inds_in]
		    seq_alt_p_over1 <- seq_alt[ind_over_1][inds_in] |> Biostrings::DNAStringSet() |>
		 						    Biostrings::subseq(sci)    |> 
		 						    Biostrings::translate()    |> 
								    suppressWarnings()
	     }
    }

    #- seq_alt_p: conatenate and restore order, drop sequences without start codon 
    ind_in <- seq_len(evr_ind_idl |> length())
    if( (n_over_1 > 0) && ( n_nover1 > 0) ){ #- both types
	#-concatenating the NULL works
  	    seq_alt_p <- c(seq_alt_p_nover1, seq_alt_p_over1)
	    ind_in    <- c(ind_nover1, ind_over_1[inds_in])
    } else if(n_over_1 > 0){ #- only start-overlapping
	    seq_alt_p <- seq_alt_p_over1
    	ind_in    <-  ind_over_1[inds_in]
    } else if(n_nover1 > 0){ #- no start-overlapping
	    seq_alt_p <- seq_alt_p_nover1
	    ind_in    <- ind_nover1
    }

    #- restory order
    ord       <- order(ind_in)
    ind_in    <- ind_in[ord]
    seq_alt_p <- seq_alt_p[ord]

    #-subset all relevant quantities to reflect omission of non-start alternatives
    evr_ind_idl <- evr_ind_idl[ind_in]
    exn_ind_idl <- exn_ind_idl[ind_in]
    seq_alt     <- seq_alt[ind_in] #- don't think we need to do that, actually

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

    #- ptc does not contain "original" stop codon.
    #  the corner cases here that we miss contain indels affecting the
    #  last reference codon
    ptc_pos <- stop_pos
    ptc_pos[stop_pos == Biostrings::width(seq_alt_p)] <- NA

    #- update exon boundaries and apply rules
    #----------------------------------------
    #- now we need the updated exon boundaries in CDS/codon coordinates and TX coordinates.
    #  since all variants are contained in one exon, we just need to update that.

    d_w <- vars[evr_ind_idl]$alt |> Biostrings::width() -
        vars[evr_ind_idl]$ref |> Biostrings::width()

    afu <- function(ptc_pos, exn_ind, d_w, num_exn, txname, 
                    css_prox_dist, penultimate_prox_dist){
        #-----------------------------------------------
        #- if ptc_pos is NA we return FALSE
        if(is.na(ptc_pos)){
            res <- rep(FALSE, 6)
            names(res) <- c("is_ptc","is_last", "is_penultimate", "is_cssProximal", "is_single", "is_407plus")
            return(res)
        }
        #- new exon boundaries (in cds/codon space), don't currently need it though
        exn_sta_p_alt <- exn_sta_p
        exn_end_p_alt <- exn_end_p
        inds          <- exn_ind : length(exn)
        delt          <- (d_w + (exn_end_p_nc[inds] %% 3) ) %/% 3

        exn_sta_p_alt[inds[-1]] <- exn_sta_p_alt[inds[-1]] + delt[-1]
        exn_end_p_alt[inds]     <- exn_end_p_alt[inds]     + delt

        #- new exon boundaries (in tx/nuc space), easier
        exn_sta_t_alt <- exn_sta_t
        exn_end_t_alt <- exn_end_t
        inds          <- exn_ind : length(exn)
        delt          <- rep(d_w, length(inds))

        exn_sta_t_alt[inds[-1]] <- exn_sta_t_alt[inds[-1]] + delt[-1]
        exn_end_t_alt[inds]     <- exn_end_t_alt[inds]     + delt

        #- find the exon that contains the PTC (not necessarily the variant)
        exn_ind_ptc <- (exn_sta_p_alt <= ptc_pos) |> which() |> max()
        #- location of ptc in that exon
        #ptc_loc <- ptc_pos - exn_sta_p[exn_ind_ptc] + 1
        apply_nmd_escape_rules(ptc_pos, exn_ind_ptc, exn_sta_t_alt[exn_ind_ptc],
                                exn_end_t_alt[exn_ind_ptc], num_exn, txname,
                                css_prox_dist = css_prox_dist, penultimate_prox_dist = penultimate_prox_dist)
    }

    tbl_idl <- purrr::pmap_dfr( list( ptc_pos = ptc_pos,
                                      exn_ind = exn_ind_idl,
                                      d_w     = d_w,
                                      num_exn = length(exn_sta_p),
				                      txname  = txname,
                                      css_prox_dist = css_prox_dist,
                                      penultimate_prox_dist = penultimate_prox_dist),
                                afu)

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
        dfr <- rbind(res_snv$tbl_snv, res_idl$tbl_idl) |> S4Vectors::DataFrame()
        S4Vectors::mcols(res)$res_aenmd <- dfr
        return(res)
    }
}


