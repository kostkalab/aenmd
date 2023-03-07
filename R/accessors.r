#' Genome information for transcript set
#' @param details Logical. If TRUE, return seqinfo object. Default: FALSE.
#' @return String (if details == FALSE). Seqinfo object otherwise
#' 
#' @export
ad_get_genome <- function() NULL

#' Transcripts considered
#' @return GRanges, one range per transcript 
#' 
#' @export
ad_get_txs <-  function() NULL

#' Transcript mask
#' @return GRanges, one range for exon - to - transcript mapping 
#' 
#' @export
ad_get_txs_mask <-  function() NULL

#' Splice mask
#' @return GRanges splice regions
#' 
#' @export
ad_get_spl_mask <-  function() NULL

#' Exons, stratified by transcript
#' @param txname String. Name of transcript
#' @return GRanges. Exons, sorted 5' to 3'
#'
#' @export
ad_get_exns_by_tx <-  function() NULL

#' Coding sequence, stratified by transcript
#' @param txname String. Name of transcript
#' @return DNAString. Coding sequence
#'
#' @export
ad_get_cds_by_tx <-  function() NULL

#' Query for single exon transcripts
#' @param txname String. Name of transcript 
#' @return Logical. TRUE if txname is a single exon transcript.
#' @details Caveat: This function returns "FALSE" for single exon transcripts that are not in the transcript set.
#' 
#' @export
ad_is_single_exn_tx <-  function() NULL

#' Query stop-making SNVs by transcript
#' @param keys Character vector. SNVs to be checked.
#' @param txname String. Transcript name. 
#' @return Logical. For each key, whether it is PTC-generating in txname.
#' 
#' @export
ad_is_ptc_snv <-  function() NULL





