#' Simple wrapper to read in in vcf file, make GRanges object.
#'
#' Uses vcfR package.
#' @param vcf_filename String. Filename of vcf file.
#' @param pass_only Logical. Keep only entries where FILTER = PASS
#' @return List. \code{vcf_rng} contains the fixed part of the vcf-file as GRanges. \code{vcf} is the complete \code{vcfR::vcfR} object.
parse_vcf_vcfR <- function(vcf_filename, pass_only = TRUE){
#==========================================================

    vcf      <- vcfR::read.vcfR(vcf_filename)
    vcf_rng  <- GenomicRanges::GRanges(vcfR::getCHROM(vcf),
                                       IRanges::IRanges(vcfR::getPOS(vcf),
                                                        vcfR::getPOS(vcf)))
    vcf_rng$ref     <- vcfR::getREF(vcf) |> Biostrings::DNAStringSet()
    vcf_rng$alt     <- vcfR::getALT(vcf) |> Biostrings::DNAStringSet()
    vcf_rng$id      <- vcfR::getID(vcf)
    vcf_rng$filter  <- vcfR::getFILTER(vcf)
    vcf_rng$qual    <- vcfR::getQUAL(vcf)
    GenomeInfoDb::seqlevelsStyle(vcf_rng) <- 'NCBI'

    if(pass_only){
        vcf_rng <- vcf_rng[vcf_rng$FILTER == 'PASS']
    }

    return(list(vcf_rng = vcf_rng, vcf = vcf))
}

#' Simple wrapper to read in in vcf file, make GRanges object.
#'
#' Uses VariantAnnotation package.
#' @param vcf_filename String. Filename of vcf file.
#' @param pass_only Logical. Keep only entries where FILTER = PASS
#' @param verbose Logicale. Report progress.
#' @return List. \code{vcf_rng} contains the fixed part of the vcf-file as GRanges. \code{vcf} is the complete \code{vcfR::vcfR} object.
parse_vcf_VariantAnnotation<- function(vcf_filename, pass_only = TRUE, verbose = TRUE){
#==========================================================

    if(verbose) message('Reading in vcf file ...', appendLF = FALSE)
    vcf        <- VariantAnnotation::readVcf(vcf_filename)
    vcf        <- VariantAnnotation::expand(vcf)
    vcf_rng    <- rowRanges(vcf) #- FIXME: how do I use rowRanges from teh VariantAnnotation::class:VCF here?
    vcf_rng$ID <- names(vcf_rng)
    if(verbose) message(' done.')

    #- don't really know why this is necessary, but sometimes does not work otherwise

    if(verbose) message('Creating GRanges ...', appendLF = FALSE)
    gn <- GenomeInfoDb::genome(vcf_rng) |> unique()
    GenomeInfoDb::genome(vcf_rng) <- NA
    GenomeInfoDb::seqlevelsStyle(vcf_rng) <- 'NCBI'
    GenomeInfoDb::genome(vcf_rng) <- gn

    if(pass_only){
        ind     <- vcf_rng$FILTER == 'PASS'
        vcf_rng <- vcf_rng[ind]
        vcf     <- vcf[ind]
    }

    #- we lower-case, clean up 
    colnames(mcols(vcf_rng)) <- colnames(mcols(vcf_rng)) |> janitor::make_clean_names()
    
if(verbose) message(' done.')

    return(list(vcf_rng = vcf_rng, vcf = vcf))
}