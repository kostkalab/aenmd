
## `aenmd` &mdash; Annotating escape from nonsense-mediated decay


#### Overview
Some genetic variants in coding sequence introduce premature termination (i.e., stop) codons, called PTCs.
Often, mRNA from such altered transcripts is degraded by a process called [nonsense-mediated decay (NMD)](https://en.wikipedia.org/wiki/Nonsense-mediated_decay). However, depending on the location of a PTC with respect to the transcript model (i.e., transcription start site, coding start site, exon boundaries, etc.) a set of rules exist to predict escape from NMD for transcripts with PTCs [CITE]. The `aenmd` package implements such rules and provides functionality to annotate human genetic variants that generate PTCs with (predicted) escape from NMD.

#### Quickstart 

##### Installing `aenmd`

```
> token1 <- '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'
> token2 <- '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'

> aenmd_data_gh <- 'kostkalab/aenmd.data.ensdb.v105'
> aenmd_gh      <- 'kostkalab/aenmd'

> devtools::install_github(aenmd_data_gh, auth_token = token1)
> devtools::install_github(aenmd_gh, auth_token = token2)
```

##### Annotating variants for escape from NMD

```
```

#### Details

##### Transcript sets and assembly

##### Rules for predicting escape from NMD

##### Input and output

##### Command line usage

#### Known Issues

With the current version the following known issues exist:

* Mitochondrial chromosomes have [alternative start/stop codons](https://en.wikipedia.org/wiki/Stop_codon#Alternative_stop_codons), also in vertebrates.\
  This is not currently implemented in `aenmd`, so that PTC calls for locations on the MT chromosome are only with respect to "normal" stop codons.

* `aenmd` filters out variants that overlap splice regions in *any* annotated transcript. Therefore, variants that overlap splice regions in some transcripts but generate PTCs in others are currently not evaluated.

* `aenmd` filters out variants that overlap splice regions, even when the splice region proper is not altered.
