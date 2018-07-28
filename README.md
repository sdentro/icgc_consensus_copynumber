# ICGC PCAWG-11 consensus copy number pipeline

This repository contains the pipeline used to create the consensus copy number profiles of the ICGC PCAWG data set. The procedure has been described in detail in https://www.biorxiv.org/content/early/2018/05/07/312041

## General procedure

The pipeline takes copy number profiles from six individual callers (ABSOLUTE, Aceseq, Battenberg, cloneHD, JaBbA and Sclust) and combines them into a single consensus profile using a multi-tiered approach. Prerequisite is that a consensus segmentation has been established and that the six callers have been run using that segmentation.

The procedure then works as follows

* Subclonal copy number is rounded, both up and down: `consensus_subclonal_copynumber.R`
* The profiles are combined into a consensus: `get_agreement.R`
* Finally, the release profile is created: `create_release_profile.R`

The runtime of the procedure is dependent upon the number of segments, but typically takes about 30 minutes on a single core.

## Dependencies

Software packages used to develop the code and run the pipeline on the PCAWG dataset. Installation of these packages should normally take a few minutes via Bioconductor.

```
R (version 3.1.0)
```

R libraries (all installed via Bioconductor)
```
Bioconductor (version 3.0)
BiocInstaller (version 1.16.5)
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

This pipeline takes input data that has been bundled and is available [here](https://www.synapse.org/#!Synapse:syn15426870) (access restricted)

* Inferred sex of the patient (processed output of the Sanger PCAWG pipeline)
* Consensus breakpoints (established by the consensus breakpoints procedure available TODO)
* Copy number profiles from the six methods with various annotations

## Setup instructions - this is work in progress

Steps

* download the data [bundle](https://www.synapse.org/#!Synapse:syn15426870)
* run the setup_environment.sh script that:
    * extract bundle
    * mkdir output output/broad_rounded_clonal output/broad_rounded_alt_clonal output/consensus_profile output/consensus_profile_final output/dkfz_rounded_alt_clonal output/dkfz_rounded_clonal output/figures output/mustonen_rounded_alt_clonal output/mustonen_rounded_clonal output/pcawg11_consensus_profile output/peifer_rounded_alt_clonal output/peifer_rounded_clonal output/saves output/status_inventory output/summary_stats output/vanloowedge_rounded_alt_clonal output/vanloowedge_rounded_clonal output/raw_ploidy
* run the pipeline (TODO make a single bash script that runs all steps per sample)

## How to run the pipeline - work in progress

```
TODO - integrate the different steps into a single pipeline
```

## Produced output
The pipeline produces two files for each sample:

A file with the publically released PCAWG consensus copy number profile. It contains a row per copy number segment and with the following columns

| Column | Description |
| --- | --- |
| chromosome | Chromosome the segment is based |
| start and end | Start and end coordinates of the segment |
| total_cn, major_cn and minor_cn | The total copy number and the copy number of the major and minor allele respectively |
| star | Quality rating where 3 represents the highest quality, 2 medium and 1 the lowest |

A file with the internally released PCAWG-11 consensus copy number profile. The copy number is exactly the same as in the above file, but more columns have been added that provide information required for other analysis within PCAWG-11. It contains a row per copy number segment and with the following columns

| Column | Description |
| --- | --- |
| chromosome | Chromosome the segment is based |
| start and end | Start and end coordinates of the segment |
| total_cn, major_cn and minor_cn | The total copy number and the copy number of the major and minor allele respectively |
| star | Quality rating where 3 represents the highest quality, 2 medium and 1 the lowest |
| level | More fine grained detail assignment of quality, please see our manuscript for more details |
| methods_agree | The number of methods that have agreed upon the copy number of this segment |
| absolute* | Additional columns provided by ABSOLUTE |
| aceseq* | Additional columns provided by Aceseq |
| battenberg* | Additional columns provided by Battenberg |
| clonehd* | Additional columns provided by cloneHD |
| jabba* | Additional columns provided by JaBbA |



