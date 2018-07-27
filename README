# ICGC PCAWG-11 consensus copy number pipeline

This repository contains the pipeline used to create the consensus copy number profiles of the ICGC PCAWG data set. The procedure has been described in detail in https://www.biorxiv.org/content/early/2018/05/07/312041

## General procedure

The pipeline takes copy number profiles from six individual callers (ABSOLUTE, Aceseq, Battenberg, cloneHD, JaBbA and Sclust) and combines them into a single consensus profile using a multi-tiered approach. Prerequisite is that a consensus segmentation has been established and that the six callers have been run using that segmentation.

The procedure then works as follows

* Subclonal copy number is rounded, both up and down: consensus_subclonal_copynumber.R
* The profiles are combined into a consensus: get_agreement.R
* Finally, the release profile is created: create_release_profile.R

## Dependencies

R libraries
```
readr
GenomicRanges
gtools
ggplot2
grid
gridExtra
parallel
stringr
reshape2
```

## Data bundle

This pipeline takes input data that has been bundled and is available TODO

* Inferred sex of the patient (processed output of the Sanger PCAWG pipeline)
* Consensus breakpoints (established by the consensus breakpoints procedure available TODO)
* Copy number profiles from the six methods with various annotations

## How to run the pipeline

```
TODO
```

## TODO setup instructions

Steps

* download the data bundle (TODO refer to a synapse id)
* run the setup_environment.sh script that:
    * extract bundle
    * mkdir output output/broad_rounded_clonal output/broad_rounded_alt_clonal output/consensus_profile output/consensus_profile_final output/dkfz_rounded_alt_clonal output/dkfz_rounded_clonal output/figures output/mustonen_rounded_alt_clonal output/mustonen_rounded_clonal output/pcawg11_consensus_profile output/peifer_rounded_alt_clonal output/peifer_rounded_clonal output/saves output/status_inventory output/summary_stats output/vanloowedge_rounded_alt_clonal output/vanloowedge_rounded_clonal output/raw_ploidy
* run the pipeline (TODO make a single bash script that runs all steps per sample)