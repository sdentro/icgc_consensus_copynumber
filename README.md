# ICGC PCAWG-11 consensus copy number pipeline

This repository contains the pipeline used to create the consensus copy number profiles of the ICGC PCAWG data set. The procedure has been described in detail in https://www.biorxiv.org/content/early/2018/05/07/312041

## General procedure

The pipeline takes copy number profiles from six individual callers (ABSOLUTE, Aceseq, Battenberg, cloneHD, JaBbA and Sclust) and combines them into a single consensus profile using a multi-tiered approach. Prerequisite is that a consensus segmentation has been established and that the six callers have been run using that segmentation.

The procedure then works as follows

* Subclonal copy number is rounded, both up and down: `round_profile.R`
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

## Setup instructions

Steps to perform in the directory where the pipeline is to be run

* download the data [bundle](https://www.synapse.org/#!Synapse:syn15426870)
* extract the bundle
* run `setup.sh` to create all the output directories
* run the pipeline using `run.sh`

## How to run the pipeline

```
/path/to/code/run.sh /path/to/code/ [samplename]
```

## Produced output
The pipeline produces two files for each sample:

A file with the publically released PCAWG consensus copy number profile in `output/consensus_profile_final`. It contains a row per copy number segment and with the following columns

| Column | Description |
| --- | --- |
| chromosome | Chromosome the segment is based |
| start and end | Start and end coordinates of the segment |
| total_cn, major_cn and minor_cn | The total copy number and the copy number of the major and minor allele respectively |
| star | Quality rating where 3 represents the highest quality, 2 medium and 1 the lowest |

A file with the internally released PCAWG-11 consensus copy number profile in `output/pcawg11_consensus_profile`. The copy number is exactly the same as in the above file, but more columns have been added that provide information required for other analysis within PCAWG-11. It contains a row per copy number segment and with the following columns

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



