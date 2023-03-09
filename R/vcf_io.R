#' Simple wrapper to read in in vcf file, make GRanges object.
#'
#' Uses vcfR package.
#' @param vcf_filename String. Filename of vcf file.
#' @param verbose Logical. Report progress.
#' @return List. \code{vcf_rng} contains the fixed part of the vcf-file as GRanges. \code{vcf_obj} is the complete \code{vcfR::vcfR} object.
parse_vcf_vcfR <- function(vcf_filename, verbose = TRUE){
#========================================================

    if(verbose) message('Reading in vcf file using VariantAnnotation::readVcf...', appendLF = FALSE)
    vcf      <- vcfR::read.vcfR(vcf_filename)
    if(verbose) message(' done.')

    if(verbose) message('Creating GRanges object...', appendLF = FALSE)
    vcf_rng  <- GenomicRanges::GRanges(vcfR::getCHROM(vcf),
                                       IRanges::IRanges(vcfR::getPOS(vcf),
                                                        vcfR::getPOS(vcf)))

    ref <- vcfR::getREF(vcf) 
    alt <- vcfR::getALT(vcf) 
    out <- (is.na(alt) | is.na(ref)) |> which()  #- yes, this happens...

    vcf_rng <- vcf_rng[-out]
    ref     <- ref[-out] |> Biostrings::DNAStringSet()
    alt     <- alt[-out] |> Biostrings::DNAStringSet()
    vcf     <- vcf[-out]

    vcf_rng$id      <- vcfR::getID(vcf)
    vcf_rng$filter  <- vcfR::getFILTER(vcf)
    vcf_rng$qual    <- vcfR::getQUAL(vcf)
    GenomeInfoDb::seqlevelsStyle(vcf_rng) <- 'NCBI'
    if(verbose) message(' done.')

    return(list(vcf_rng = vcf_rng, vcf_obj = vcf))
}

#' Simple wrapper to read in in vcf file, make GRanges object.
#'
#' Uses VariantAnnotation package.
#' @param vcf_filename String. Filename of vcf file.
#' @param verbose Logical. Report progress.
#' @param ... dotdotdot. Passed to \code{VariantAnnotation::readVcf}.
#' @return List. \code{vcf_rng} contains the fixed part of the vcf-file as GRanges. \code{vcf_obj} is the complete \code{vcfR::vcfR} object.
#' @importClassesFrom VariantAnnotation VCF
parse_vcf_VariantAnnotation<- function(vcf_filename, verbose = TRUE, ...){
#=========================================================================

    if(verbose) message('Reading in vcf file using VariantAnnotation::readVcf...', appendLF = FALSE)
    vcf        <- VariantAnnotation::readVcf(vcf_filename, ...)
    vcf        <- VariantAnnotation::expand(vcf)
    vcf_rng    <- SummarizedExperiment::rowRanges(vcf) #- FIXME: how do I use rowRanges from VariantAnnotation::class:VCF here?
    vcf_rng$ID <- names(vcf_rng)
    if(verbose) message(' done.')

    #- don't really know why this is necessary, but sometimes does not work otherwise

   # if(verbose) message('Creating GRanges ...', appendLF = FALSE)
   # gn <- GenomeInfoDb::genome(vcf_rng) |> unique()
   # GenomeInfoDb::genome(vcf_rng) <- NA
   # GenomeInfoDb::seqlevelsStyle(vcf_rng) <- 'NCBI'
   # GenomeInfoDb::genome(vcf_rng) <- gn

   #if(pass_only){
    #    ind     <- vcf_rng$FILTER == 'PASS'
    #    vcf_rng <- vcf_rng[ind]
    ##    vcf     <- vcf[ind]
    #}

    #- we lower-case, clean up 
    colnames(S4Vectors::mcols(vcf_rng)) <- colnames(S4Vectors::mcols(vcf_rng)) |> janitor::make_clean_names()
    
if(verbose) message(' done.')

    return(list(vcf_rng = vcf_rng, vcf_obj = vcf))
}