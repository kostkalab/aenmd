
#' Filters / processes variants
#' @param vcf_rng GRanges object of variants. Start is first base of variant.
#'            Each variant is of length 1.
#'            Each variant has an ``ref`` and ``alt`` metadata column with DNAStringSets of the sequences.
#' @param verbose Logical. Report progress.
#' @return List.
#' @importFrom utils head tail
#' @examples
#' @export
process_variants <- function(vcf_rng, 
                             verbose   = TRUE){

    if(verbose) message("Processing variants.")

    #- expand the vcf ranges to contain first and last base of reference allele
    if( ! all(GenomicRanges::width(vcf_rng) == Biostrings::width(vcf_rng$ref)) ){
        vcf_rng <-  GenomicRanges::resize(vcf_rng,
                                      width = Biostrings::width(vcf_rng$ref),
                                      fix="start")
    }

    #- filter out variants that overlap splice sites
    #  we filter out ranges where either start or end overlap a splice region.
    if(verbose) message("Filtering out splice variants.")
    vcf_rng <- filter_splice_overlap(vcf_rng)

    #- make sure asll variants are *contained* within *one* of our exons
    #- step 1: they need to *overlap* coding sequence
    if(verbose) message("Filtering out variants that are not in coding sequence.")
    vcf_rng <- filter_in_coding_exon(vcf_rng)

    #- Classify ins, del, sbs, snv; assign key;
    vcf_rng$type                                                                   <- NULL
    vcf_rng$type[Biostrings::width(vcf_rng$ref) >  Biostrings::width(vcf_rng$alt)] <- 'del'
    vcf_rng$type[Biostrings::width(vcf_rng$ref) <  Biostrings::width(vcf_rng$alt)] <- 'ins'
    vcf_rng$type[Biostrings::width(vcf_rng$ref) == Biostrings::width(vcf_rng$alt)] <- 'sbs'
    vcf_rng$type[(vcf_rng$type == 'sbs') & (Biostrings::width(vcf_rng$alt) == 1)]  <- 'snv'
    vcf_rng$key <- make_keys(vcf_rng)
    

    #- for snvs, only use "stop-making" snvs
    if(verbose) message("Filtering out snvs that don't create stop codons.")
    vcf_rng <- filter_ptc_generating_snvs(vcf_rng)

    #- check for duplicated variants
    #- Happens for example in the GRCh38-lifted-over version of gnomad v 2.2.1
    #- when GRCh37 SNVs map to GRCh38 SNVs with different RSIDS (e.g, rs782209862 and rs61732915)
    if( !(length(vcf_rng$key) == length(unique(vcf_rng$key))) ){
        stop("duplicated variants")
    }

    return(vcf_rng)
}

#' Wrapper function for annotating a set of variants with NMD escape
#' @param vcf_rng GRanges object of variants. Start is first base of variant.
#'            Each variant is of length 1.
#'            Each variant has an ``ref`` and ``alt`` metadata column with DNAStrings of the sequences.
#' @param css_prox_dist  Distance to the CSS defining NMD escape regions.
#' @param penultimate_prox_dist Integer. Distance to the penultimate exon 3'end defining NMD escape regions.
#' @param verbose Logical. Report progress.
#' @param multicore Logical. Should multiple cores be used in the computations (via the parallel pacakge).
#' @param num_cores Integer. The number of cores to use. Only when \code{multicore = TRUE}.
#' @param rettype Character. Should a `GRangesList` be returned (default, `rettype = 'grl'`) or shold results be collated into a `GRanges` object (`rettype = 'gr'`).
#' @return GRanges Object, containing ranges in \code{vcf_rng} that overlap transcripts in \code{aenmd}'s transcript set.
#' The metadata column of this object contains a column named \code{res_aenmd} that contains annotation results.
#' @details Coming up.
#' @importFrom utils head tail
#' @examples
#' @export
annotate_nmd <- function(vcf_rng, css_prox_dist = 150L, penultimate_prox_dist = 50L,
                             verbose = FALSE , multicore = FALSE, num_cores = 2, rettype = 'grl'){
#============================================================================================================================

    #- stratify variants by tx
    vri_by_tx <- stratify_by_tx(vcf_rng)

    #- annotate variants for each transcript        
    if(multicore == FALSE){
        res <- pbapply::pblapply(seq_len( vri_by_tx |> length()),
                                 function(i) {
                                    res <- annotate_variants_by_tx(names(vri_by_tx)[i], vcf_rng[vri_by_tx[[i]]],
                                                                         css_prox_dist = css_prox_dist, 
                                                                         penultimate_prox_dist = penultimate_prox_dist)
                                    res})
    } else {
        res <- parallel::mclapply(seq_len( vri_by_tx |> length()),
                                   function(i) {
                                    res <- annotate_variants_by_tx(names(vri_by_tx)[i], vcf_rng[vri_by_tx[[i]]],
                                                                         css_prox_dist = css_prox_dist, 
                                                                         penultimate_prox_dist = penultimate_prox_dist)
                                    res },
				                  mc.cores = num_cores)
    }

    names(res) <- names(vri_by_tx)
    if(rettype == "raw") return(res)
    if (rettype == "grl") {
	    minw    <- lapply(res, function(x) min(GenomicRanges::width(x))) |> unlist()
        ind_out <- which(minw == 0) |> sort()
	if(length(ind_out) >0){
        	res <- res[-ind_out]  |> GenomicRanges::GRangesList()
	} else {
		res <- res |>  GenomicRanges::GRangesList()
	}
    } else if (rettype == "gr") {
        #res <- res |> GenomicRanges::GRangesList() |> unlist()
        #- remove empty ranges
        minw    <- lapply(res, function(x) min(GenomicRanges::width(x)))
        ind_out <- which(minw == 0) |> sort()
        #- return results as GRanges
        if (length(ind_out) > 0) {
            res <- res[-ind_out] |> GenomicRanges::GRangesList() |> unlist()
        } else {
            res <- res |> GenomicRanges::GRangesList() |> unlist()
        }
        names(res) <- paste0(res$res_aenmd$transcript, '|', res$key)
        res <- sort(res)
    }

    return(res)
}
