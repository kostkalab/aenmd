
#' Read in vcf file, make GRanges object.
#'
#' Uses vcfR package. Keep only bare-bone info from vcf.
#' @param vcf_filename String. Filename of vcf file.
#' @param pass_only Logical. Keep only entries where FILTER = PASS
parse_vcf <- function(vcf_filename, pass_only = TRUE){
#=====================================================

    vcf      <- vcfR::read.vcfR(vcf_filename)
    vcf_rng  <- GenomicRanges::GRanges(vcfR::getCHROM(vcf),
                                       IRanges::IRanges(vcfR::getPOS(vcf),
                                                        vcfR::getPOS(vcf)))
    vcf_rng$ref     <- vcfR::getREF(vcf) |> Biostrings::DNAStringSet()
    vcf_rng$alt     <- vcfR::getALT(vcf) |> Biostrings::DNAStringSet()
    vcf_rng$id      <- vcfR::getID(vcf)
    vcf_rng$filter  <- vcfR::getFILTER(vcf)

    GenomeInfoDb::seqlevelsStyle(vcf_rng) <- 'NCBI'

    if(pass_only){
        vcf_rng <- vcf_rng[vcf_rng$filter == 'PASS']
    }

    return(vcf_rng)
}

#' Filters / processes variants
#' @param vcf_rng GRanges object of variants. Start is first base of variant.
#'            Each variant is of length 1.
#'            Each variant has an ``ref`` and ``alt`` metadata column with DNAStringSets of the sequences.
#' @param check_ref Logical. Should the reference alleles be checked against assembly.
#' @param verbose Logical. Report progress.
#' @return List.
#' @importFrom utils head tail
#' @examples
#' @export
process_variants <- function(vcf_rng, check_ref = FALSE, verbose = TRUE){

    if(verbose) message("Processing variants.")

    #- using GRCh38 only for now
    hsap <- BSgenome.Hsapiens.NCBI.GRCh38::Hsapiens

    #- expand the vcf ranges to contain first and last base of refrence allele
    vcf_rng <-  GenomicRanges::resize(vcf_rng,
                                      width = Biostrings::width(vcf_rng$ref),
                                      fix="start")

    #- check if reference variants are legal
    if(check_ref){
        if(verbose) message("Confirming reference alleles match assembly.")
        seqs <- BSgenome::getSeq(hsap, vcf_rng)
        if( ! all(as.character(seqs) == vcf_rng$ref) ){
            stop("Refrence variants annoations don't all match assembly.")
        }
    }

    #- filter out variants that overlap splice sites
    #  we filter out ranges where either start or end overlap a splice region.
    if(verbose) message("Filtering out splice variants.")
    ov_starts <- GenomicRanges::findOverlaps(GenomicRanges::resize(vcf_rng, width=1L, fix='start'),
                                             future::value(._EA_spl_grl))
    ov_ends   <- GenomicRanges::findOverlaps(GenomicRanges::resize(vcf_rng, width=1L, fix='end'),
                                             future::value(._EA_spl_grl))

    out_idx  <- sort(c(unique(S4Vectors::queryHits(ov_starts)),
                       unique(S4Vectors::queryHits(ov_ends))))
    if(length(out_idx)>0){
        vcf_rng  <- vcf_rng[-out_idx]
    }

    #- make sure asll variants are *contained* within *one* of our exons
    if(verbose) message("Filtering out splice variants that are not contained in at least one exon.")
    ovc     <- GenomicRanges::countOverlaps(vcf_rng, future::value(._EA_exn_grl), type="within")
    out_lgl <- ovc == 0
    if(length(out_idx)>0){
        vcf_rng  <- vcf_rng[!out_lgl]
    }

    #- Classify ins, del, sbs, assign key; snv for
    vcf_rng$type                                                                   <- NA
    vcf_rng$type[Biostrings::width(vcf_rng$ref) >  Biostrings::width(vcf_rng$alt)] <- 'del'
    vcf_rng$type[Biostrings::width(vcf_rng$ref) <  Biostrings::width(vcf_rng$alt)] <- 'ins'
    vcf_rng$type[Biostrings::width(vcf_rng$ref) == Biostrings::width(vcf_rng$alt)] <- 'sbs'
    vcf_rng$type[(vcf_rng$type == 'sbs') & (Biostrings::width(vcf_rng$alt) == 1)]  <- 'snv'

    starts <- GenomicRanges::start(vcf_rng) |> stringr::str_pad(9L, pad="0")
    vcf_rng$key <- paste0(GenomicRanges::seqnames(vcf_rng), ":", starts,"|" ,vcf_rng$ref, "|", vcf_rng$alt)

    #- for snvs, only use "stop-making" snvs
    if(verbose) message("Filtering out snvs that don't create stop codons.")
    ind <- rep(TRUE, length(vcf_rng))
    tps <- vcf_rng$type
    kys <- vcf_rng$key
    mtc <- triebeard::longest_match(future::value(._EA_snv_tri), kys)
    #- exclude SNVs that do not generate stop codons
    ind[ (is.na(mtc)) & (vcf_rng$type == 'snv') ] <- FALSE
    vcf_rng <- vcf_rng[ind]

    #- check that we only have unique variants
    if( !(length(vcf_rng$key) == length(unique(vcf_rng$key))) ){
        stop("duplicated variants")
    }

    return(vcf_rng)
}

#' Wrapper function for annotating a set of variants with NMD escape
#' @param vcf_rng GRanges object of variants. Start is first base of variant.
#'            Each variant is of length 1.
#'            Each variant has an ``ref`` and ``alt`` metadata column with DNAStrings of the sequences.
#' @param check_ref Logical. Should the reference alleles be checked against assembly.
#' @param verbose Logical. Report progress.
#' @return List.
#' @importFrom utils head tail
#' @examples
#' @export
annotate_nmd <- function(vcf_rng, check_ref = TRUE, verbose = FALSE){
#====================================================================

    if(verbose) message("Processing variants.")

    #vcf_rng <- process_variants(vcf_rng, check_ref, verbose)

    #- explode variants into variant/transcript pairs
    #  FIXME: we could make the transcript set more flexible...
    #         maybe pass it as a function arg, provide more than just tsl=1
    if(verbose) message("Creating variant-transcript pairs.")
    ov                 <- GenomicRanges::findOverlaps(vcf_rng,
                                                      future::value(._EA_exn_grl))
    vcf_rng_by_tx      <- vcf_rng[S4Vectors::queryHits(ov)] #- 405,094
    vcf_rng_by_tx$enst <- names(future::value(._EA_exn_grl))[S4Vectors::subjectHits(ov)]

    #- TODO: why is this still getting vars? Should be in process_variant already...
    if(length(vcf_rng_by_tx) == 0){
        message("No transcript-overlapping variants.")
        return(vcf_rng_by_tx)
    }

    #- actually annotate variants
    #- TODO: add parallel capability.
    if(verbose) message("Annotating variant-trainscript pairs for NMD escape.")
    rr <- lapply(seq_len(length(vcf_rng_by_tx)), function(ind) annotate_variant(vcf_rng_by_tx[ind],vcf_rng_by_tx[ind]$enst ))
    rmat <- unlist(lapply(rr,function(x) x[[1]])) |> matrix(byrow=TRUE, ncol=6)
    colnames(rmat) <- rr[[1]]$rules |> names()
    vcf_rng_by_tx$rules <- rmat

    if(verbose) message("Done.")
    return(vcf_rng_by_tx)
}

#' Wrapper function for annotating a set of variants with NMD escape
#' @param vcf_rng GRanges object of variants. Start is first base of variant.
#'            Each variant is of length 1.
#'            Each variant has an ``ref`` and ``alt`` metadata column with DNAStrings of the sequences.
#' @param check_ref Logical. Should variants be verified against refernce sequence.
#' @param verbose Logical. Report progress.
#' @return List.
#' @importFrom utils head tail
#' @examples
#' @export
annotate_nmd_v2 <- function(vcf_rng, check_ref = FALSE, verbose = FALSE){
#=======================================================================

    #- connect variants to exons
    ov      <- GenomicRanges::findOverlaps(vcf_rng,
                                           future::value(._EA_exn_grl))
    sH  <- S4Vectors::subjectHits(ov)
    qH  <- S4Vectors::queryHits(ov)
    sHu <- sort(unique(sH))

    #- stratify variants by overlapped transcripts
    afu <- function(sHi){
        txname <- names(future::value(._EA_exn_grl))[sHi]
        res    <- vcf_rng[qH[sH == sHi]]
        S4Vectors::mcols(res)$tx_id <- txname
        return(res)
    }
    rlst        <- lapply(sHu, afu)
    names(rlst) <- names(future::value(._EA_exn_grl))[sHu]

    #- annotate variants for each transcript
    res <- lapply(seq_len( rlst |> length()),
                  \(i) annotate_variants_by_tx(names(rlst)[i], rlst[[i]]))

    res <- res |> GenomicRanges::GRangesList() |> unlist()

    return(res)
}


if(FALSE){
#- why so slow.. better way?

#- collect all hits by transcripts
#--------------------------------

vcf_rng <- parse_vcf("../../dat/gnomad.exomes.r2.1.1.sites.liftover_grch38_chr21.vcf")
vcf_rng <- process_variants(vcf_rng)
}
#'
#' txname <- names(rlst)[1]
#' evars  <- rlst[[1]]
#'
#'
#' if( is.null( seq <- get0(txname, future::value(._EA_cds_env)) )){
#'     return(NULL)
#' }
#'
#' if( is.null( exn <- get0(txname, future::value(._EA_exn_env)) )){
#'     return(NULL)
#' }
#'
#' #cdn_bsv <- Biostrings::codons(seq)
#' #cdn_vec <- cdn_bsv |> as.character()
#'
#' ov <- findOverlaps(evars,exn)
#' qH <- queryHits(ov)
#' sH <- subjectHits(ov)
#'
#' #- up till now: only one exon per variant
#' #. i.e., multiple exon spanning variants are discarded.
#'
#' if( (qH |> max() |> table()) > 1) stop("multiple-exon variants are not supported.")
#' #- we can now look up for each variant its exon
#' sHm        <- subjectHits(ov)
#' names(sHm) <- queryHits(ov)
#'
#' #-exon starts and ends (genome, "transcript", protein)
#' #-----------------------------------------------------
#' #- note: we could pre-caclulate that if it turns out slow, but unlikely.
#' exn_sta_g <- start(exn)
#' exn_end_g <- end(exn)
#' exn_sta_t <- c(1, cumsum(width(exn))[-length(exn)]+1)
#' exn_end_t <- cumsum(width(exn))
#' exn_sta_p <- (exn_sta_t -1) %/% 3 + 1
#' exn_end_p <- (exn_end_t -1) %/% 3 + 1
#' #- which nucleotide in codon: 1, 2 or 3
#' exn_sta_p_nc <- (exn_sta_t -1) %% 3 + 1
#' exn_end_p_nc <- (exn_end_t -1) %% 3 + 1
#'
#' #------
#' #- SNVs
#' #------
#' #- All SNVs are PTC-generating (others have been filtered out)
#' #- Exon boundaries don't change.
#' evr_ind_snvs <- (evars$type == "snv") |> which()
#' exn_ind_snvs <- sHm[evr_ind_snvs |> as.character()]
#' #- get the codons
#' if(all(strand(exn)=="-")){
#'     #- index of codon in "current" exon
#'     cdn_ind_curr <- ( exn_end_g[exn_ind_snvs] - start(evars[evr_ind_snvs]) +
#'                       exn_sta_p_nc[exn_ind_snvs] -1 ) %/% 3 + 1
#'     cdn_ind      <- cdn_ind_curr + exn_sta_p[exn_ind_snvs] - 1
#' } else{
#'     cdn_ind_curr <- ( start(evars[evr_ind_snvs]) - exn_sta_g[exn_ind_snvs] +
#'                       exn_sta_p_nc[exn_ind_snvs] -1 ) %/% 3 + 1
#'     cdn_ind      <- cdn_ind_curr + exn_sta_p[exn_ind_snvs] - 1
#' }
#'
#' tbl_snv <- purrr::pmap_dfr( list( ptc_loc = cdn_ind_curr,
#'                        exn_ind = exn_ind_snvs,
#'                        exn_sta = exn_sta_p[exn_ind_snvs],
#'                        exn_end = exn_end_p[exn_ind_snvs],
#'                        num_exn = length(exn_sta_p)),
#'                  get_rules)
#'
#' if(detailed){
#'     tbl_snv <- tbl_snv |> dplyr::mutate(#exon_index = exn_ind_snvs,
#'                                 ptc_position_alt_p = cdn_ind,
#'                                 tx_id = exn$tx_id[exn_ind_snvs],
#'                                 exn_id = exn$exon_id[exn_ind_snvs])
#' }
#'
#' #------------------------
#' #- INSERTIONS / DELETIONS
#' #------------------------
#' #- Here we need to make the new protein to find the PTC position
#' #  (that's only strictly true for insertions, but okay)
#'
#' #
#' evr_ind_idl <- (evars$type != 'snv') |> which()
#' exn_ind_idl <- sHm[evr_ind_idl |> as.character()]
#'
#' #- get the alternative version of the DNA sequence
#'
#' seq_ref <- future::value(._EA_cds_env)[[txname]]
#'
#' #- need to map the genomic variants to the reference protein (CDS/nuc) coordinates
#'
#' if(all(strand(exn)=="-")){
#'     ref_nuc_sta <- exn_end_g[exn_ind_idl] - end(evars[evr_ind_idl])     + exn_sta_t[exn_ind_idl]
#'     ref_nuc_end <- exn_end_g[exn_ind_idl] - start(evars[evr_ind_idl])   + exn_sta_t[exn_ind_idl]
#'
#' } else{
#'     ref_nuc_sta <- start(evars[evr_ind_idl]) - exn_sta_g[evr_ind_idl] + exn_sta_t[exn_ind_idl]
#'     ref_nuc_end <- end(evars[evr_ind_idl])   - exn_sta_g[evr_ind_idl] + exn_sta_t[exn_ind_idl]
#' }
#'
#' ind_out     <- (ref_nuc_sta <= 3) #- those overlap the start codon (but not splice regions)
#' evr_ind_idl <- evr_ind_idl[!ind_out]
#' exn_ind_idl <- sHm[evr_ind_idl |> as.character()]
#' ref_nuc_sta <- ref_nuc_sta[!ind_out]
#' ref_nuc_end <- ref_nuc_end[!ind_out]
#'
#' alt_vrs <- evars$alt[evr_ind_idl]
#' if(all(strand(exn)=="-")) alt_vrs <- Biostrings::reverseComplement(alt_vrs)
#' make_alt <- function(ref_nuc_sta, ref_nuc_end, alt, seq_ref){
#'     Biostrings::xscat(Biostrings::subseq(seq_ref, 1, ref_nuc_sta -1), alt,
#'                       Biostrings::subseq(seq_ref, ref_nuc_end +1, length(seq_ref)))[[1]]
#' }
#'
#' seq_alt <- sapply(seq_len(length(evr_ind_idl)),
#'                   function(ind) make_alt(ref_nuc_sta = ref_nuc_sta[ind],
#'                                          ref_nuc_end[ind],
#'                                          alt_vrs[ind],
#'                                          seq_ref) )
#'
#' seq_alt_p <- seq_alt |> DNAStringSet() |> Biostrings::translate() |> suppressWarnings()
#' stop_pos  <- Biostrings::vmatchPattern("*",seq_alt_p)
#' stop_pos  <- lapply(stop_pos, function(x) {
#'     if(length(x) == 0) {
#'         fst_stop <- NA
#'     } else {
#'         fst_stop <- x[1] |> Biostrings::start() |> as.integer()
#'     }
#'     }) |> unlist()
#'
#' #- ptc does not contain "original" stop codon.
#' #  the corner cases here that we miss contain indels affecting the
#' #  last reference codon
#'
#' ptc_pos <- stop_pos
#' ptc_pos[stop_pos == Biostrings::width(seq_alt_p)] <- NA
#'
#' #- now we need the updated exon boundaries in CDS/codon coordinates.
#' #  since all variants are contained in one exon, we just need to update that.
#'
#' d_w <- evars[evr_ind_idl]$alt |> Biostrings::width() -
#'        evars[evr_ind_idl]$ref |> Biostrings::width()
#'
#' afu <- function(ptc_pos, exn_ind, d_w, num_exn){
#' #-----------------------------------------------
#'     #- if ptc_pos is NA we return FALSE
#'     if(is.na(ptc_pos)){
#'         res <- rep(FALSE, 6)
#'         names(res) <- c("is_ptc","is_last", "is_penultimate", "is_first", "is_single", "is_407plus")
#'         return(res)
#'     }
#'     #- need to make new exon boundaries (in cds/codon space)
#'     exn_sta_p_alt <- exn_sta_p
#'     exn_end_p_alt <- exn_end_p
#'     inds          <- exn_ind : length(exn)
#'     delt          <- (d_w + (exn_end_p_nc[inds] %% 3) ) %/% 3
#'
#'     exn_sta_p_alt[inds[-1]] <- exn_sta_p_alt[inds[-1]] + delt[-1]
#'     exn_end_p_alt[inds]     <- exn_end_p_alt[inds]     + delt
#'
#'     #- find the exon that contains the PTC (not necessarily the variant)
#'     exn_ind_ptc <- (exn_sta_p_alt <= ptc_pos) |> which() |> max()
#'     #- location of ptc in that exon
#'     ptc_loc <- ptc_pos - exn_sta_p[exn_ind_ptc] + 1
#'     get_rules(ptc_loc, exn_ind_ptc, exn_sta_p_alt[exn_ind_ptc],
#'               exn_end_p_alt[exn_ind_ptc], num_exn)
#' }
#'
#' tbl_idl <- purrr::pmap_dfr( list( ptc_pos = ptc_pos,
#'                                   exn_ind = exn_ind_idl,
#'                                   d_w     = d_w,
#'                                   num_exn = length(exn_sta_p)),
#'                             afu)
#'
#' if(detailed){
#'     tbl_idl <- tbl_idl |> dplyr::mutate(
#'         ptc_position_alt_p = stop_pos,
#'         tx_id = exn$tx_id[exn_ind_idl],
#'         exn_id = exn$exon_id[exn_ind_idl])
#' }
#'
#' #- now string both results together
#'
#'
#'
#'
#'
#' #- details: both sequences
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' se <- start(exn)
#' ee <- end(exn)
#' if(all(strand(exn)=='+')) bp_map <- lapply(seq_len(length(se)),function(ind) seq(se[ind],ee[ind]))
#' if(all(strand(exn)=='-')) bp_map <- lapply(seq_len(length(se)),function(ind) seq(se[ind],ee[ind]) |> rev())
#' bp_map               <- bp_map |> unlist()
#' bp_map_cdn           <- matrix(bp_map, nrow=3)
#' colnames(bp_map_cdn) <- cdn_vec
#'
#' #- codon coordinates in transcript (CDS in our case)
#' bp_tra               <- seq_len(length(seq))
#' bp_tra_cdn           <- matrix(bp_tra, nrow=3)
#' colnames(bp_tra_cdn) <- cdn_vec
#'
#'
#'
#'
#' evr_ind_idl <- (evars$type != 'snv') |> which()
#'
#' alt_var <- evars$alt[evr_ind_idl]
#'
#' if(all(GenomicRanges::strand(exn) == '-')) {
#'     alt_var <- Biostrings::reverseComplement(alt_var)
#' }
#' seq_alt <- Biostrings::xscat(Biostrings::subseq(seq_ref, 1, sae[1]-1), alt_var,
#'                              Biostrings::subseq(seq_ref, sae[2]+1, length(seq_ref)))[[1]]
#'
#'
#'
#'
#' #' Determine NMD escape rules
#' #'
#' #' @param ptc_loc Integer. Location of PTC. In CDS codon coordinates, 1-based.
#' #' @param exn_ind Integer. Index of the exon the variant overlaps, 1-based.
#' #' @param exn_sta Integer. Start of the exon. In CDS codon coordinates, 1-based.
#' #' @param exn_end Integer. End of the exon. In CDS codon coordinates, 1-based.
#' #' @param num_exn Integer. Number of exons in transcript.
#' get_rules <- function(ptc_loc, exn_ind, exn_sta, exn_end, num_exn){
#' #------------------------------------------------------------------
#'
#'     res <- rep(FALSE, 6)
#'     names(res) <- c("is_ptc","is_last", "is_penultimate", "is_first", "is_single", "is_407plus")
#'
#'     #- didnt' get a PTC
#'     #------------------
#'     if(is.na(ptc_loc) | is.null(ptc_loc)) return(res)
#'     res["is_ptc"] <- TRUE
#'
#'     #- are we in the first 150 bp of the first exon?
#'     if( (exn_ind ==1) & (ptc_loc <= 50) )  res["is_first"] <- TRUE
#'
#'     #- are we in the last exon?
#'     if(exn_ind == num_exn) res["is_last"] <- TRUE
#'
#'     #- are we in a single exon transcrtipt
#'     if(num_exn == 1) res["is_single"] <- TRUE
#'
#'     #- the exon longer than 407bp?
#'     if( (exn_end - exn_sta) >= 136 ) res["is_407plus"] <- TRUE
#'
#'     #- are we in (the last 51bp of) the penultimate exon?
#'     cnd1 <- exn_ind == (num_exn-1)
#'     cnd2 <- (exn_end - ptc_loc) <= 17
#'     if( cnd1 & cnd2) res["is_penultimate"] <- TRUE
#'
#'     return(res)
#' }
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' #- codon cooordinates in genome
#' se <- start(exn)
#' ee <- end(exn)
#' if(all(strand(exn)=='+')) bp_map <- lapply(seq_len(length(se)),function(ind) seq(se[ind],ee[ind]))
#' if(all(strand(exn)=='-')) bp_map <- lapply(seq_len(length(se)),function(ind) seq(se[ind],ee[ind]) |> rev())
#' bp_map               <- bp_map |> unlist()
#' bp_map_cdn           <- matrix(bp_map, nrow=3)
#' colnames(bp_map_cdn) <- cdn_vec
#'
#' #- codon coordinates in transcript (CDS in our case)
#' bp_tra               <- seq_len(length(seq))
#' bp_tra_cdn           <- matrix(bp_tra, nrow=3)
#' colnames(bp_tra_cdn) <- cdn_vec
#'
#'
#'
#'
#'
#' #- SNVs [exon boundaries for all the same]; We know they generate PTCs.
#'
#' #- INSERTIONS (more sequence in alt) [individual exon boundaries]
#'
#' #- DELETIONS (less sequence in alt) [individual exon boundaries]
#'
#'
#' }
