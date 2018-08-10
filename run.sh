#!/bin/bash
libpath=$1
item=$2
sex=`grep ${item} data_bundle/sex.txt | cut -f 7`
Rscript ${libpath}/consensus_subclonal_copynumber.R ${libpath} ${item} output/ TRUE FALSE ${sex};
Rscript ${libpath}/consensus_subclonal_copynumber.R ${libpath} ${item} output/ FALSE FALSE ${sex};
Rscript ${libpath}/get_agreement.R ${libpath} ${item} output/ ${sex};
Rscript ${libpath}/create_release_profile.R ${libpath} ${item} output/ ${sex};
