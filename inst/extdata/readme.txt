#- make a small vcf from clinvar, here it is documented

#- strip info
bcftools annotate -x INFO clinvar_20221211.vcf.gz > clinvar_20221211_noinfo.vcf
bgzip clinvar_20221211_noinfo.vcf

#-get subsample
R -e '  library(vcfR);
        vcf <- read.vcfR("clinvar_20221211_noinfo.vcf.gz");
        set.seed(1223214);
        vcf_s <- vcf[sample(seq_len(dim(vcf)[1]), 1000) |> sort() ];
        write.vcf(vcf_s, file="clinvar_20221211_noinfo_sample1k.vcf")'
bgzip clinvar_20221211_noinfo_sample1k.vcf

