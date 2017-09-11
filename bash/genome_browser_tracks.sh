#!/bin/bash

# Webserver URL: http://cbdm-01.zdv.uni-mainz.de/~jibnsale

# example track:
# http://cbdm-01.zdv.uni-mainz.de/~jibnsale/tracks/hg19/factorbookMotifs/factorbookMotifPos.bed.CTCF.sort.bed.gz

# http://cbdm-01.zdv.uni-mainz.de/~jibnsale/tracks/hg19/chromloop/v04_selected_models.motifSig6_w1000_b1.motifs.bed.sorted.bed.gz

#=======================================================================
# set some variables here:
#=======================================================================

export WEBSERVER="jibnsale@cbdm-01.zdv.uni-mainz.de"
export WEBSERVER_PORT="jibnsale@cbdm-01.zdv.uni-mainz.de:22022"
export WEB_DIR="/home/jibnsale/public_html/tracks/hg19/chromloop"
export WEB_URL="http://cbdm-01.zdv.uni-mainz.de/~jibnsale/tracks/hg19/chromloop"
export BEDTOOLS=data/bin/bedtools2/bin/bedtools
export TABIX_DIR="data/bin/tabix-0.2.6"

#=======================================================================
# Copy motifs to genome browser.
#=======================================================================

# create folder on webserver
ssh -p 22022 ${WEBSERVER} mkdir ${WEB_DIR}

# function to prepare .bed file (using tabix and upload to webserver)
# see: http://wiki.wubrowse.org/Simple_bed
function upload_bed {
  
  FILE=$1
  
  # sort bed file
  ${BEDTOOLS} sort -i ${FILE} \
    > ${FILE}.sort.bed
  
  # gunzip file
  ${TABIX_DIR}/bgzip -f ${FILE}.sort.bed
  
  # create index
  ${TABIX_DIR}/tabix -p bed ${FILE}.sort.bed.gz
  
  # copy motifs to webserver
  scp -P 22022 ${FILE}.sort.bed.gz \
               ${FILE}.sort.bed.gz.tbi \
    ${WEBSERVER}:${WEB_DIR}/
  
  # report link
  echo ${WEB_URL}/$(basename ${FILE}).sort.bed.gz
}
export -f upload_bed


upload_bed results/v04_selected_models.motifSig6_w1000_b1.motifs.bed

#=======================================================================
# copy ChIP-seq tracks to browser
#=======================================================================

# get all .bigWig file paths
cut -f 3 results/v04_selected_models.motifSig6_w1000_b1.meta.tsv \
  |tail -n +2 \
  |while read FILE ; do
    scp -P 22022 ${FILE} ${WEBSERVER}:${WEB_DIR}/
  done
  

cut -f 2,52 results/v04_selected_models.motifSig6_w1000_b1.meta.tsv
