
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

#' Wrapper function for annotating a set of variants with NMD escape
#' @param vcf_rng GRanges object of variants. Start is first base of variant.
#'            Each variant is of length 1.
#'            Each variant has an ``ref`` and ``alt`` metadata column with DNAStrings of the sequences.
#' @param check_ref Logical. Should the reference alleles be checked against assembly.
#' @return List.
#' @importFrom utils head tail
#' @examples
#' @export
annotate_nmd <- function(vcf_rng, check_ref = TRUE){
#===================================================

    #- using GRCh38 only for now
    hsap <- BSgenome.Hsapiens.NCBI.GRCh38::Hsapiens

    #- expand the vcf ranges to contain first and last base of refrence allele
    vcf_rng <-  GenomicRanges::resize(vcf_rng,
                                      width = Biostrings::width(vcf_rng$ref),
                                      fix="start")

    #- check if reference variants are legal
    if(check_ref){
        seqs <- BSgenome::getSeq(hsap, vcf_rng)
        if( ! all(as.character(seqs) == vcf_rng$ref) ){
            stop("Refrence variants annoations don't all match assembly")
        }
    }

    #- filter out variants that overlap splice sites
    # FIXME: this filters too much, e.g. when (sets of) exons get deleted.
    #        -> filter only when fist/last nucleotide overlap splice region?
    ov       <- GenomicRanges::findOverlaps(vcf_rng, ._EA_spl_grl)
    out_idx  <- sort(unique(S4Vectors::queryHits(ov)))
    if(length(out_idx)>0){
        vcf_rng  <- vcf_rng[-out_idx]
    }

    #- Classify ins, del, sbs, assign key
    vcf_rng$type                                                                   <- NA
    vcf_rng$type[Biostrings::width(vcf_rng$ref) >  Biostrings::width(vcf_rng$alt)] <- 'del'
    vcf_rng$type[Biostrings::width(vcf_rng$ref) <  Biostrings::width(vcf_rng$alt)] <- 'ins'
    vcf_rng$type[Biostrings::width(vcf_rng$ref) == Biostrings::width(vcf_rng$alt)] <- 'sbs'
    vcf_rng$key <- paste(as.character(vcf_rng),vcf_rng$ref,vcf_rng$alt,sep="|")

    #- check that we only have unique variants
    if( !(length(vcf_rng$key) == length(unique(vcf_rng$key))) ){
        stop("duplicated variants")
    }

    #- explode variants into variant/transcript pairs
    #  FIXME: we could make the transcript set more flexible...
    #         maybe pass it as a function arg, provide more than just tsl=1
    ov                 <- GenomicRanges::findOverlaps(vcf_rng, ._EA_txs_grl)
    vcf_rng_by_tx      <- vcf_rng[S4Vectors::queryHits(ov)] #- 405,094
    vcf_rng_by_tx$enst <- names(._EA_txs_grl)[S4Vectors::subjectHits(ov)]

    #- actually annotate variants
    ##res <- annotate_vars(vcf_rng_by_tx)
    rr <- lapply(seq_len(length(vcf_rng_by_tx)), function(ind) annotate_variant(vcf_rng_by_tx[ind],vcf_rng_by_tx[ind]$enst ))
    rmat <- unlist(lapply(rr,function(x) x[[1]])) |> matrix(byrow=TRUE, ncol=6)
    colnames(rmat) <- rr[[1]]$rules |> names()
    vcf_rng_by_tx$rules <- rmat

    return(vcf_rng_by_tx)
}
