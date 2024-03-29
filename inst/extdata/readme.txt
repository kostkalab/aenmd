#- make a small vcf from clinvar, here it is documented

#- strip info
#============
bcftools annotate -x INFO clinvar_20221211.vcf.gz > clinvar_20221211_noinfo.vcf
bgzip clinvar_20221211_noinfo.vcf

#-get subsample
#==============

R -e '  library(vcfR);
        vcf <- read.vcfR("clinvar_20221211_noinfo.vcf.gz");
        set.seed(1223214);
        vcf_s <- vcf[sample(seq_len(dim(vcf)[1]), 1000) |> sort() ];
        write.vcf(vcf_s, file="clinvar_20221211_noinfo_sample1k.vcf.gz")'
gunzip clinvar_20221211_noinfo_sample1k.vcf.gz
bgzip clinvar_20221211_noinfo_sample1k.vcf

#- get grch37 version
#====================

wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2022/clinvar_20221211.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2022/clinvar_20221211.vcf.gz.tbi
mv clinvar_20221211.vcf.gz clinvar_grch37_20221211.vcf.gz
mv clinvar_20221211.vcf.gz.tbi clinvar_grch37_20221211.vcf.gz.tbi

#- remove info 
#=============

bcftools annotate -x INFO clinvar_grch37_20221211.vcf.gz > clinvar_grch37_20221211_noinfo.vcf
bgzip clinvar_grch37_20221211_noinfo.vcf

#-get subsample
#==============

R -e '  library(vcfR);
        vcf <- read.vcfR("clinvar_grch37_20221211_noinfo.vcf.gz");
        set.seed(1223214);
        vcf_s <- vcf[sample(seq_len(dim(vcf)[1]), 1000) |> sort() ];
        write.vcf(vcf_s, file="clinvar_grch37_20221211_noinfo_sample1k.vcf.gz")'
gunzip clinvar_grch37_20221211_noinfo_sample1k.vcf.gz
bgzip clinvar_grch37_20221211_noinfo_sample1k.vcf

#- fix ClinVar headers (for both files)
#======================================

R -e '  library(vcfR);
	source("./fix_clinvar_header.R")
	fix_clinvar_header(vcf_in  = "clinvar_grch37_20221211_noinfo_sample1k.vcf.gz",
			   vcf_out = "clinvar_grch37_20221211_noinfo_sample1k_fxd.vcf");
	fix_clinvar_header(vcf_in  = "clinvar_20221211_noinfo_sample1k.vcf.gz",
			   vcf_out = "clinvar_20221211_noinfo_sample1k_fxd.vcf")'
bgzip clinvar_grch37_20221211_noinfo_sample1k_fxd.vcf
bgzip clinvar_20221211_noinfo_sample1k_fxd.vcf
mv clinvar_20221211_noinfo_sample1k_fxd.vcf.gz clinvar_20221211_noinfo_sample1k.vcf.gz
mv clinvar_grch37_20221211_noinfo_sample1k_fxd.vcf.gz clinvar_grch37_20221211_noinfo_sample1k.vcf.gz


