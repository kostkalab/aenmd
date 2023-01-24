## `aenmd` &mdash; annotating escape from nonsense-mediated decay

This package provides functionality to annotate variant-transcript pairs with predicted escape from nonsense-mediated decay (NMD).
Details can be found in the package documentation, and in our preprint [here](www.here.com). The `aenmd` package depends on transcript annotations that are provided via data packages that can be found in the [aenmd_data](https://github.com/kostkalab/aenmd_data) repository on GitHub. 


#### Installing `aenmd`
```
#- need the remotes package
if(!require(remotes){
    install.packages("remotes")
}

#- install data package (current default is ENSEMBL v 105)
remotes::install_github("kostkalab/aenmd_data/aenmd.data.ensdb.v105")

#- install aenmd package
remotes::install_github("kostkalab/aenmd_data/aenmd")
```

#### Getting started
Below is a minimal example:
```
library(aenmd) #- automatically loads annotations also

#- annotate a SNV
#----------------

#- make a GRanges object
library(GenomicRanges)
library(Biostrings)
snv     <- GRanges("1:100210750-100210750:")
snv$ref <- Biostrings::DNAStringSet("G")
snv$alt <- Biostrings::DNAStringSet("A")
snv$type <- "snv"
snv$key <- aenmd:::make_keys(snv)

#- annotate
snv_annotated <- annotate_nmd(snv, rettype="gr")

#- inspect results: snv causes a PTC in ENST00000370132,
#- but it does not trigger any rules annotating it as 
#- being predicted to escape NMD in this transcript
res <- snv_annotated$res_aenmd
res
#      is_ptc   is_last is_penultimate is_cssProximal is_single is_407plus
#   <logical> <logical>      <logical>      <logical> <logical>  <logical>
# 1      TRUE     FALSE          FALSE          FALSE     FALSE      FALSE
#        transcript
#       <character>
# 1 ENST00000370132

```

This works also with >1 variants; variants should be normalized:

```
snv2     <- GRanges("4:169406633-169406633:")
snv2$ref <- Biostrings::DNAStringSet("C")
snv2$alt <- Biostrings::DNAStringSet("A")
snv2$type <- "snv"
snv2$key <- aenmd:::make_keys(snv2)

snvs <- c(snv,snv2) |> suppressWarnings()
snvs_annotated <- annotate_nmd(snvs, rettype="gr")
snvs_annotated$res_aenmd
#- snv2 overlaps multiple transcripts, reflected in an expanded result
# DataFrame with 6 rows and 7 columns
#      is_ptc   is_last is_penultimate is_cssProximal is_single is_407plus
#   <logical> <logical>      <logical>      <logical> <logical>  <logical>
#1       TRUE     FALSE          FALSE          FALSE     FALSE      FALSE
#2       TRUE     FALSE          FALSE          FALSE     FALSE      FALSE
#3       TRUE     FALSE          FALSE          FALSE     FALSE      FALSE
#4       TRUE     FALSE          FALSE          FALSE     FALSE      FALSE
#5       TRUE     FALSE          FALSE          FALSE     FALSE      FALSE
#6       TRUE     FALSE          FALSE          FALSE     FALSE      FALSE
#        transcript
#       <character>
# 1 ENST00000370132
# 2 ENST00000439128
# 3 ENST00000507142
# 4 ENST00000510533
# 5 ENST00000511633
# 6 ENST00000512193


```

#### Known Issues

With the current version the following known issues exist:

* Mitochondrial chromosomes have [alternative start/stop codons](https://en.wikipedia.org/wiki/Stop_codon#Alternative_stop_codons), also in vertebrates.\
  This is not currently implemented in `aenmd`, so that PTC calls for locations on the MT chromosome are only with respect to "normal" stop codons.

* `aenmd` filters out variants that overlap splice regions in *any* annotated transcript. Therefore, variants that overlap splice regions in some transcripts but generate PTCs in others are currently not evaluated.


