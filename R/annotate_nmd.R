
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
#' @param check_ref Logical. Should variants be verified against refernce sequence.
#' @param verbose Logical. Report progress.
#' @param multicore Logical. Should multiple cores be used in the computations (via the parallel pacakge).
#' @param num_cores Integer. The number of cores to use. Only when \code{multicore = TRUE}.
#' @param rettype Character. Should a `GRangesList` be returned (default, `rettype = 'grl'`) or shold results be collated into a `GRanges` object (`rettype = 'gr'`).
#' @return List.
#' @details For multicore, the `BiocParallel` default backand returned by `BiocParallel::bpparam()` is used.
#' @importFrom utils head tail
#' @examples
#' @export
annotate_nmd_v2 <- function(vcf_rng, check_ref = FALSE, verbose = FALSE , multicore = FALSE, num_cores = 2, rettype = 'grl'){
#============================================================================================================================

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
    if(multicore == FALSE){
        rlst  <- lapply(sHu, afu)
    } else {
        rlst  <- parallel::mclapply(sHu, afu, mc.cores = num_cores)
    }
    names(rlst) <- names(future::value(._EA_exn_grl))[sHu]

    #- annotate variants for each transcript
    if(multicore == FALSE){
        res <- pbapply::pblapply(seq_len( rlst |> length()),
                                 \(i) annotate_variants_by_tx(names(rlst)[i], rlst[[i]]))
    } else {
        res <- parallel::mclapply(seq_len( rlst |> length()),
                                      \(i) annotate_variants_by_tx(names(rlst)[i], rlst[[i]]),
				      mc.cores = num_cores)
    }

    names(res) <- names(rlst)
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
	    names(res) <- NULL #- transcript names do not make sense here.
        }
    }

    return(res)
}
