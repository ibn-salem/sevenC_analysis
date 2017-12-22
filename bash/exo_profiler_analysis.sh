#!/bin/bash

#=======================================================================
# Analyse ChIP-nexus data using Exo-Profiler tool at loops and non looping sites
#=======================================================================
SAMTOOLS="data/bin/bin/samtools"

LOOP_BED="results/CTCF_JASPAR.v01.pval_2.5e-06.loop_anchors.bed"
NON_LOOP_BED="results/CTCF_JASPAR.v01.pval_2.5e-06.non-loop_anchors.bed"

EXOPROFILER_DIR="../ExoProfiler"

mkdir -p results/exo_profiler

SAMPLE="SRR2312570"
for SAMPLE in "SRR2312570" "SRR2312571" ; do
  
  BAM=data/Tang2015/${SAMPLE}.flexcat_matched_barcode.fastq.bowtie2.hg19_filtered.bam
  
  # sort and index bam file
  $SAMTOOLS sort $BAM > $BAM.sorted.bam
  $SAMTOOLS index $BAM.sorted.bam
  
  python ${EXOPROFILER_DIR}/python/5primeCounter.py \
    --input_sites $LOOP_BED \
    --input_format bed \
    --bam_file $BAM.sorted.bam  \
    --output_prefix results/exo_profiler/${SAMPLE}_loop_ExoProfile

  python ${EXOPROFILER_DIR}/python/5primeCounter.py \
    --input_sites $NON_LOOP_BED \
    --input_format bed \
    --bam_file $BAM.sorted.bam  \
    --output_prefix results/exo_profiler/${SAMPLE}_non-loop_ExoProfile

  Rscript ${EXOPROFILER_DIR}/R/exoPlotter.R results/exo_profiler/${SAMPLE}_loop_ExoProfile
  Rscript ${EXOPROFILER_DIR}/R/exoPlotter.R results/exo_profiler/${SAMPLE}_non-loop_ExoProfile

  
done
