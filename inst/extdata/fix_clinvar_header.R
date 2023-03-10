#- there is an issue vcfR and VariantAnnoation have with the ClinVar header.
#  the offending line is 
#
#  ##ID=<Description="ClinVar Variation ID">
#
#  ...there is no ID= in "<...>", as in ##ID=<ID=...,Description="ClinVar Variation ID">;
#  we'll change it to
#
#  ##ID="ClinVar Variation ID"
#
#  ...which does not make any trouble. 

library(VariantAnnotation)
fix_clinvar_header <- function(vcf_in, vcf_out){

	#- extract vcf header
	vcf2 <- VariantAnnotation::readVcf(vcf_in)
	hd   <- VariantAnnotation::header(vcf2)

	#- fix ID entry
	colnames(meta(hd)[[5]]) <- 'Value'

	#- write out vcf with fixed header
	VariantAnnotation::header(vcf2) <- hd
	VariantAnnotation::writeVcf(vcf2, filename = vcf_out, index = FALSE)
}


