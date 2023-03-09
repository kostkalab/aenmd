## aenmd 

### R package for annotating predicted escape from nonsense-mediated decay (NMD)

- [Introduction](#introduction)
- [Installing aenmd](#installing-aenmd)
- [Getting started](#getting-started)
- [Known issues](#known-issues)

#### Introduction

This `R` package provides functionality to annotate transcripts with genetic variants that introduce premature termination codons (stop codons) with predicted escape from NMD. More details can be found in the package documentation, and in our preprint [here](www.here.com). 

The `aenmd` package depends on transcript annotations that are provided via data packages that can be found in the [aenmd_data](https://github.com/kostkalab/aenmd_data) repository on GitHub. Further on, a command line interface/utility ([aenmd_cli](https://github.com/kostkalab/aenmd_cli)) is available on GitHub as well.


#### Installing `aenmd`
```R
#- need the remotes package
if(!require(remotes){
    install.packages("remotes")
}

#- install data package (current default is ENSEMBL v 105)
remotes::install_github("kostkalab/aenmd_data/aenmd.data.ensdb.v105")

#- install aenmd package
remotes::install_github("kostkalab/aenmd")
```

#### Getting started
Below is a short example, where we use a small vcf file provided with the `aenmd` package. 

##### 1. Import variants from vcf file
First, we load variants we'd like to annotate. `aenmd` uses `GenomicRanges` with `ref` and `alt` metadata columns to represent variants.
```R
#- load aenmd package
library(aenmd) #- automatically loads annotation data as well

#- load variants to annotate (1,000 random ClinVar variants)
vcf_file <- system.file('extdata/clinvar_20221211_noinfo_sample1k.vcf.gz', package = 'aenmd')
vcf      <- aenmd:::parse_vcf_VariantAnnotation(vcf_file)
vcf_rng  <- vcf$vcf_rng
vcf_rng
#GRanges object with 1000 ranges and 5 metadata columns:
#         seqnames              ranges strand | param_range_id                     ref            alt      qual      filter
#            <Rle>           <IRanges>  <Rle> |       <factor>          <DNAStringSet> <DNAStringSet> <numeric> <character>
#     [1]        1       940501-941150      * |             NA GGAGCCTGCA...CAGATCTCCT              G        NA           .
#     [2]        1       942504-942505      * |             NA                      CG              C        NA           .
#     [3]        1             1041417      * |             NA                       C              T        NA           .
#...
```

##### 2. Filter out variants we don't annotate
Next, we filter out variants that:
  - overlap splice regions
  - don't overlap coding sequence
  - are SNVs but don't generate stop codons

We also generate a unique key for each variant.

```R
#- Variant filtering
vcf_rng_fil <- process_variants(vcf_rng)
vcf_rng_fil
#GRanges object with 70 ranges and 7 metadata columns:
#       seqnames              ranges strand | param_range_id                     ref            alt      qual      filter        type                    key
#          <Rle>           <IRanges>  <Rle> |       <factor>          <DNAStringSet> <DNAStringSet> <numeric> <character> <character>            <character>
#   [1]        1            12011551      * |             NA                       C              T        NA           .         snv        1:012011551|C|T
#   [2]        1            99902707      * |             NA                       C              T        NA           .         snv        1:099902707|C|T
#   [3]        1           115727014      * |             NA                       C              A        NA           .         snv        1:115727014|C|A
#...
```

##### 3. Annotate remaining variant/transcript pairs
For the remaining variants we annotate each overlapping transcript:

```R
#- annotate
vcf_rng_ann <- annotate_nmd(vcf_rng_fil, rettype="gr")
vcf_rng_ann
#GRanges object with 142 ranges and 8 metadata columns:
#                                  seqnames    ranges strand | param_range_id            ref            alt      qual      filter        type             key            res_aenmd
#                                     <Rle> <IRanges>  <Rle> |       <factor> <DNAStringSet> <DNAStringSet> <numeric> <character> <character>     <character>          <DataFrame>
#  ENST00000235329|1:012011551|C|T        1  12011551      * |             NA              C              T        NA           .         snv 1:012011551|C|T  TRUE:TRUE:FALSE:...
#  ENST00000294724|1:099902707|C|T        1  99902707      * |             NA              C              T        NA           .         snv 1:099902707|C|T TRUE:FALSE:FALSE:...
#  ENST00000361915|1:099902707|C|T        1  99902707      * |             NA              C              T        NA           .         snv 1:099902707|C|T TRUE:FALSE:FALSE:...
#...
```

Note the `res_aenmd` column in the metadata of `vcf_rng_ann.` This contains `aenmd`'s results. 
Also, each range is named by the transcript-variant combination analyzed.
Below, we explore `aenmd`'s results

```R
vcf_rng_ann$res_aenmd
#DataFrame with 142 rows and 7 columns
#       is_ptc   is_last is_penultimate is_cssProximal is_single is_407plus      transcript
#    <logical> <logical>      <logical>      <logical> <logical>  <logical>     <character>
#1        TRUE      TRUE          FALSE          FALSE     FALSE      FALSE ENST00000235329
#2        TRUE     FALSE          FALSE          FALSE     FALSE      FALSE ENST00000294724
#3        TRUE     FALSE          FALSE          FALSE     FALSE      FALSE ENST00000361915
#4        TRUE     FALSE          FALSE          FALSE     FALSE      FALSE ENST00000370163
#5        TRUE     FALSE          FALSE          FALSE     FALSE      FALSE ENST00000370165
#...       ...       ...            ...            ...       ...        ...             ...
#138      TRUE     FALSE          FALSE          FALSE     FALSE      FALSE ENST00000288447
#139      TRUE     FALSE          FALSE          FALSE     FALSE      FALSE ENST00000357033
#140      TRUE     FALSE          FALSE          FALSE     FALSE      FALSE ENST00000447523
#141     FALSE     FALSE          FALSE          FALSE     FALSE      FALSE ENST00000303391
#142     FALSE     FALSE          FALSE          FALSE     FALSE      FALSE ENST00000453960
```
The first column indicates whether a variant causes a premature termination codon (PTC) in a transcript. The next five columns indicate rules as to whether the variant/transcript combination is expected to escape NMD. The last column gives the transcript name.

Above we see that the last two variant/transcript combinations do not have a PTC; the first combination is predicted to escape NMD, because the PTC is located in the last exon.



#### Known Issues

With the current version the following known issues exist:

* Mitochondrial chromosomes have [alternative start/stop codons](https://en.wikipedia.org/wiki/Stop_codon#Alternative_stop_codons), also in vertebrates.\
  This is not currently implemented in `aenmd`, so that PTC calls for locations on the MT chromosome are only with respect to "normal" stop codons.

* `aenmd` filters out variants that overlap splice regions in *any* annotated transcript. Therefore, variants that overlap splice regions in some transcripts but generate PTCs in others are currently not evaluated.


