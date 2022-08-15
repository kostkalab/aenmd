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

        gr1     <- GenomicRanges::intersect(gr1, gr2, ignore.strand = TRUE)
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


