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
    exn     <- future::value(._EA_exn_env)[[tx]] #- exons of the transcript;
    seq_ref <- future::value(._EA_cds_env)[[tx]] #- reference allele sequence;

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
    seq_ref_p <- suppressWarnings(Biostrings::translate(seq_ref))
    seq_alt_p <- suppressWarnings(Biostrings::translate(seq_alt))
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

        if(cnd1 || cnd2){
            is_ptc <- TRUE
        }
    }

    #- single exon rule
    is_single <- FALSE
    if(length(S4Vectors::metadata(seq_alt_p)$exon_starts) == 1){
        is_single <- TRUE
    }


    #- last exon rule (TRUE if first PTC overlaps last exon)
    is_last <- FALSE
    if(!is.na(fst_stop) && (!is_single)){
        cnd1 <- fst_stop >=  ( S4Vectors::metadata(seq_alt_p)$exon_starts |> tail(n=1) )
        if(cnd1){
            is_last <- TRUE
        }
    }

    #- penultimate exon rule (true if PTS overlaps last 50bp of penultimate exon)
    #- note that we use 17 AA, which is 51 bp and not 50
    is_penultimate <- FALSE
    if(!is.na(fst_stop) && (!is_single)){

        cnd1 <- fst_stop >=  (S4Vectors::metadata(seq_alt_p)$exon_starts |> tail(n=2))[1]
        cnd2 <- fst_stop <   (S4Vectors::metadata(seq_alt_p)$exon_starts |> tail(n=2))[2]
        cnd3 <- (( (S4Vectors::metadata(seq_alt_p)$exon_starts |> tail(n=2))[2]) - fst_stop) <= 17

        if(cnd1 && cnd2 && cnd3){
            is_penultimate <- TRUE
        }
    }

    #- first exon rule
    is_first <- FALSE
    if(!is.na(fst_stop) && (!is_single)){
        #- We are in the first exon
        cnd1 <- fst_stop <  (S4Vectors::metadata(seq_alt_p)$exon_starts |> head(n=2))[2]
        #-We are in the first 150 nucleotides (=50 AA) of the first exon
        cnd2 <- fst_stop <= 50
        if(cnd1 && cnd2){
            is_first <- TRUE
        }
    }

    #- 407 bp rule: if PTC is in large exon: escape.
    #  we do 408bp or bigger
    is_407plus <- FALSE
    if(!is.na(fst_stop)){
        #-exon starts
        tmp_es <- S4Vectors::metadata(seq_alt_p)$exon_starts
        if(is_single){
            if(length(seq_alt_p) >= 136) is_407plus <- TRUE
        } else {
            #- start of overlapping exon and the next exon
            strt_over    <- tmp_es[(tmp_es <= fst_stop) |> which() |> max()]
            strt_next    <- NA
            if(strt_over == max(tmp_es)){
                strt_next <- length(seq_alt_p) + 1
            } else {
                strt_next <- tmp_es[(tmp_es >= fst_stop) |> which() |> min()]
            }

            sze    <- strt_next - strt_over
            if(sze >= 136) is_407plus <- TRUE
        }
    }

    rule_res        <- c(is_ptc, is_last, is_penultimate, is_first, is_single, is_407plus)
    names(rule_res) <- c("is_ptc","is_last", "is_penultimate", "is_first", "is_single", "is_407plus")

    res <- list( rules   = rule_res,
                 seq_ref = seq_ref_p,
                 seq_alt = seq_alt_p)

    return(res)
}

#' Determine NMD escape rules
#'
#' @param ptc_loc Integer. Location of PTC. In CDS codon coordinates, 1-based.
#' @param exn_ind Integer. Index of the exon the variant overlaps, 1-based.
#' @param exn_sta Integer. Start of the exon. In CDS codon coordinates, 1-based.
#' @param exn_end Integer. End of the exon. In CDS codon coordinates, 1-based.
#' @param num_exn Integer. Number of exons in the CDS of the transcript.
#' @param txname  Character. Ensembl transcript id. 
#' @return Logical. Vogical vector with six named elements.
#' @details 
#' Return value (Logical vector) has six named elements.
#' - \code{is_ptc = TRUE} if there is a PTC in the alternative transcript. 
#' - \code{is_last = TRUE} if the PTC is in the last (3'-most) coding exon.
#' - \code{is_penultimate = TRUE} if the PTC is within the last 50bp of the penultimate exon.
#' - \code{is_cssProximal  = TRUE} if the PTC is within the first (5'-most) 150bp of the translation/coding start site. 
#' - \code{is_single = TRUE} if the PTC is in a single exon transcript.
#' - \code{is_407plus = TRUE} if the PTC is in an exon that is longer than 407bp.

get_rules <- function(ptc_loc, exn_ind, exn_sta, exn_end, num_exn, txname){
    #----------------------------------------------------------------------

    res <- rep(FALSE, 6)
    names(res) <- c("is_ptc","is_last", "is_penultimate", "is_cssProximal", "is_single", "is_407plus")

    #- didnt' get a PTC
    #------------------
    if(is.na(ptc_loc) || is.null(ptc_loc)) return(res)
    res["is_ptc"] <- TRUE

    #- are we in the first 150 bp proximal to the css?
    if( ptc_loc <= 50 )  res["is_cssProximal"] <- TRUE

    #- are we in the last (coding!) exon?
    if(exn_ind == num_exn) res["is_last"] <- TRUE

    #- are we in a single exon transcrtipt
    #- need to use ._EA_set_env b/c num_exn only counts the coding exons, here we need all
    if( exists(txname, where = future::value(._EA_set_env)) ) res["is_single"] <- TRUE

    #- the exon longer than 407bp?
    if( (exn_end - exn_sta) >= 136 ) res["is_407plus"] <- TRUE

    #- are we in (the last 51bp of) the penultimate exon?
    cnd1 <- exn_ind == (num_exn-1)
    cnd2 <- (exn_end - ptc_loc) <= 17
    if( cnd1 && cnd2) res["is_penultimate"] <- TRUE

    return(res)
}

#' Annotate variants overlapping a transcript
#'
#' @param txname Character. Transcript name.
#' @param vars GRanges. Overlapping variants.
#' @param detailed Logical. Should additional information be returned.
#' @return GRanges.
annotate_variants_by_tx <- function( txname, vars, detailed = FALSE){
#====================================================================

    #- get exons for the transcript
    if( is.null( exn <- get0(txname, future::value(._EA_exn_env)) )){
        stop(past0("Cannot find exons for ", txname))
    }

    #- find exon for each variant
    #- up till now: only one exon per variant
    #. i.e., multiple exon spanning variants are discarded.
    ov <- GenomicRanges::findOverlaps(vars,exn)
    qH <- S4Vectors::queryHits(ov)
    sH <- S4Vectors::subjectHits(ov)
    if( (qH |> max() |> table()) > 1) stop("multiple-exon variants are not supported.")
    #- we can now look up for each variant its exon
    sHm        <- sH
    names(sHm) <- qH

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
    #- SNVs
    #------
    #- All SNVs are PTC-generating (others have been filtered out)
    #- Exon boundaries don't change.
    evr_ind_snvs <- (vars$type == "snv") |> which()
    exn_ind_snvs <- sHm[evr_ind_snvs |> as.character()]
    #- get the codons
    if(all(GenomicRanges::strand(exn)=="-")){
        #- index of codon in "current" exon
        cdn_ind_curr <- ( exn_end_g[exn_ind_snvs] - GenomicRanges::start(vars[evr_ind_snvs]) +
                              exn_sta_p_nc[exn_ind_snvs] -1 ) %/% 3 + 1
        cdn_ind      <- cdn_ind_curr + exn_sta_p[exn_ind_snvs] - 1
    } else{
        cdn_ind_curr <- ( GenomicRanges::start(vars[evr_ind_snvs]) - exn_sta_g[exn_ind_snvs] +
                              exn_sta_p_nc[exn_ind_snvs] -1 ) %/% 3 + 1
        cdn_ind      <- cdn_ind_curr + exn_sta_p[exn_ind_snvs] - 1
    }

    tbl_snv <- purrr::pmap_dfr( list( ptc_loc = cdn_ind,
                                      exn_ind = exn_ind_snvs,
                                      exn_sta = exn_sta_p[exn_ind_snvs],
                                      exn_end = exn_end_p[exn_ind_snvs],
                                      num_exn = length(exn_sta_p),
				      txname  = txname),
                                get_rules)

    if(detailed){
        tbl_snv <- tbl_snv |> dplyr::mutate(#exon_index = exn_ind_snvs,
            ptc_position_alt_p = cdn_ind,
            tx_id = exn$tx_id[exn_ind_snvs],
            exn_id = exn$exon_id[exn_ind_snvs])
    }
 
    #------------------------
    #- INSERTIONS / DELETIONS
    #------------------------
    #- Here we need to make the new protein to find the PTC position
    #  (that's only strictly true for insertions, but okay)

    evr_ind_idl <- (vars$type != 'snv') |> which()
    exn_ind_idl <- sHm[evr_ind_idl |> as.character()]


    #- we are done here.
    if(length(evr_ind_idl) == 0){	
    	res <- vars[evr_ind_snvs]
    	dfr <- tbl_snv |> S4Vectors::DataFrame()
    	S4Vectors::mcols(res)$res_aenmd <- dfr
    	return(res)
    }



    #- get the alternative version of the DNA sequence for each variant
    #------------------------------------------------------------------
    if( is.null( seq_ref <- get0(txname, future::value(._EA_cds_env)) )){
        stop(past0("Cannot find sequence for ", txname))
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
    #- now we need the updated exon boundaries in CDS/codon coordinates.
    #  since all variants are contained in one exon, we just need to update that.

    d_w <- vars[evr_ind_idl]$alt |> Biostrings::width() -
        vars[evr_ind_idl]$ref |> Biostrings::width()

    afu <- function(ptc_pos, exn_ind, d_w, num_exn, txname){
        #-----------------------------------------------
        #- if ptc_pos is NA we return FALSE
        if(is.na(ptc_pos)){
            res <- rep(FALSE, 6)
            names(res) <- c("is_ptc","is_last", "is_penultimate", "is_cssProximal", "is_single", "is_407plus")
            return(res)
        }
        #- need to make new exon boundaries (in cds/codon space)
        exn_sta_p_alt <- exn_sta_p
        exn_end_p_alt <- exn_end_p
        inds          <- exn_ind : length(exn)
        delt          <- (d_w + (exn_end_p_nc[inds] %% 3) ) %/% 3

        exn_sta_p_alt[inds[-1]] <- exn_sta_p_alt[inds[-1]] + delt[-1]
        exn_end_p_alt[inds]     <- exn_end_p_alt[inds]     + delt

        #- find the exon that contains the PTC (not necessarily the variant)
        exn_ind_ptc <- (exn_sta_p_alt <= ptc_pos) |> which() |> max()
        #- location of ptc in that exon
        #ptc_loc <- ptc_pos - exn_sta_p[exn_ind_ptc] + 1
        get_rules(ptc_pos, exn_ind_ptc, exn_sta_p_alt[exn_ind_ptc],
                  exn_end_p_alt[exn_ind_ptc], num_exn, txname)
    }

    tbl_idl <- purrr::pmap_dfr( list( ptc_pos = ptc_pos,
                                      exn_ind = exn_ind_idl,
                                      d_w     = d_w,
                                      num_exn = length(exn_sta_p),
				      txname  = txname),
                                afu)

    if(detailed){
        tbl_idl <- tbl_idl |> dplyr::mutate(
            ptc_position_alt_p = stop_pos,
            tx_id = exn$tx_id[exn_ind_idl],
            exn_id = exn$exon_id[exn_ind_idl]) #- this is not the ptc exon but the variant exon
    }
    #- return empty range if we don't have any variants (example: all variants in tx overlap start, none has a stop)
    if( (c(evr_ind_snvs, evr_ind_idl) |> length()) == 0){
        empty_range <- GenomicRanges::GRanges('chr0',IRanges::IRanges(start=0,width=0))
        return(empty_range)
    }

    res <- vars[c(evr_ind_snvs, evr_ind_idl)]
    dfr <- rbind(tbl_snv,tbl_idl) |> S4Vectors::DataFrame()
    S4Vectors::mcols(res)$res_aenmd <- dfr
    return(res)
}
