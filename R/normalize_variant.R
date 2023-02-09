#' Left-normalize genetic variants
#'
#' Convert variants into left-normalized version (PMID: 25701572)
#' @param rng GRanges object representing a variant. Range here is *first* and *last* base of reference variant.
#' @param ref DNAString. The reference allele of the variant
#' @param alt DNAString. The alternative allele of the variant
#' @param prefix DNAString. Nucleotides preceeding the variant, used for normalization.
#' @return The normalized variant (as GRanges object)
#' @examples
#' prefix     <- Biostrings::DNAString("GATAT")
#' ref        <- Biostrings::DNAString("GAGATAT")
#' alt        <- Biostrings::DNAString("GAGAT")
#' start      <- 17
#' end        <- 17+6
#' rng        <- GenomicRanges::GRanges("1",IRanges::IRanges(17,17+6))
#' rng$ref    <- Biostrings::DNAStringSet(ref)
#' rng$alt    <- Biostrings::DNAStringSet(alt)
#' rng$prefix <- Biostrings::DNAStringSet(prefix)
#' normalize_var(rng, rng$ref, rng$alt, rng$prefix)
#' @export
normalize_var <- function(rng, ref = NULL, alt = NULL , prefix = NULL){
#======================================================================

    if(is.null(ref))    ref    <- rng$ref[[1]]
    if(is.null(alt))    alt    <- rng$alt[[1]]
    if(is.null(prefix)) prefix <- rng$prefix[[1]]

    start   <- GenomicRanges::start(rng)
    end     <- GenomicRanges::end(rng)
    ref_p   <- c(prefix[-1], ref) #- keep the first base in reserve
    alt_p   <- c(prefix[-1], alt)
    start   <- start - length(prefix) +1
    min_len <- min(length(ref_p), length(alt_p))

    #- shorten right side of alt and ref if the same
    for(i in seq_len(min_len)){
        if( utils::tail(alt_p,1) == utils::tail(ref_p,1)){
            alt_p <- utils::head(alt_p,-1)
            ref_p <- utils::head(ref_p,-1)
            end   <- end - 1
        } else {
            break
        }
        if( (length(ref_p) == 0 ) || (length(alt_p) ==0 ) ){
            ref_p <- c(prefix[1],ref_p);
            alt_p <- c(prefix[1],alt_p);
            end   <- end - 1
            break
        }
    }

    #- shorten the left side
    for(i in seq_len(min_len)){
        if( utils::head(alt_p,1) == utils::head(ref_p,1)){
            cut_base <- utils::head(alt_p,1)
            alt_p    <- utils::tail(alt_p,-1)
            ref_p    <- utils::tail(ref_p,-1)
            start    <- start + 1
        } else {
            break
        }
        if( (length(ref_p) == 0 ) || (length(alt_p) ==0 ) ){
            ref_p <- c(cut_base,ref_p);
            alt_p <- c(cut_base,alt_p);
            start <- start - 1
            break
        }
    }

    GenomicRanges::start(rng) <- start
    GenomicRanges::end(rng)   <- end
    rng$ref    <- Biostrings::DNAStringSet(ref_p)
    rng$alt    <- Biostrings::DNAStringSet(alt_p)
    if('prefix' %in% colnames(GenomicRanges::mcols(rng))) rng$prefix <- NULL

    return(rng)
}

#' Left-normalize genetic variants
#'
#' Convert variants into left-normalized version (PMID: 25701572);
#' Naive (i.e., ,slow) "vectorization" of the ``normalize_var`` function.
#' @param rng GRanges object representing a variant. Range here is *first* and *last* base of refernce variant.
#' Needs to have ``ref``, ``alt`` and ``prefix`` in the metadata.
#' @return The normalized variant (as GRanges object)
#' @examples
#' @export
normalize_vars <- function(rng){
  res <- sapply( seq_len(length(rng)) ,
                 function(ind) normalize_var(rng[ind])) |>
                               GenomicRanges::GRangesList() |>
                               unlist()
}
