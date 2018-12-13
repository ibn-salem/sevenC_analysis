#!/bin/bash

#===============================================================================
# This script documents all downlaods and raw data processing in the data folder.
#
# It is supposed to be run from within the 'data' directory
#
# >cd data
# >sh download.sh
#
#===============================================================================
#module load tools/parallel/20170622

# set some variables here:
export BIN=bin
export Q_DIR=../../Q
export Q=${Q_DIR}/bin/Q
export BEDTOOLS=${BIN}/bedtools2/bin/bedtools
export FASTQ_DUMP="bin/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump"
export SAMTOOLS=${BIN}"/bin/samtools"
export FASTQC=${BIN}/FastQC/fastqc
mkdir -p ${BIN}

#-----------------------------------------------------------------------
# General genome assembly based data and tools from UCSC:
#-----------------------------------------------------------------------

# UCSC liftover chains
mkdir -p UCSC
wget -P UCSC http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz
wget -P UCSC http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
gunzip UCSC/*.gz

# get chromosome sizes
wget -P UCSC http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
head -n 24 UCSC/hg19.chrom.sizes > UCSC/hg19.chrom.sizes.real_chroms


# download liftOver tool from UCSC:
wget -P ${BIN} http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod u+x ${BIN}/liftOver

# download bedGraphToBigWig tool from UCSC
wget -P ${BIN} http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
chmod u+x ${BIN}/bedGraphToBigWig

# download twoBitToFa tool
wget -P ${BIN} http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
chmod u+x ${BIN}/twoBitToFa

#-----------------------------------------------------------------------
# download and compile BEDtools:
#-----------------------------------------------------------------------
mkdir -p ${BIN}
wget -P ${BIN} wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz
tar xvfz  ${BIN}/bedtools-2.26.0.tar.gz -C ${BIN}
make -C ${BIN}/bedtools2


#=======================================================================
# Download specific data sets
#=======================================================================

#=======================================================================
# JASPAR 2018 CTCF motifs from UCSC GenomeBrowser track
#=======================================================================
mkdir -p JASPAR2018

wget -P JASPAR2018 http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/hg19/tsv/MA0139.1.tsv.gz
gunzip JASPAR2018/MA0139.1.tsv.gz

#=======================================================================
# Hi-C data from Rao et al 2014 Cell
#=======================================================================
mkdir -p Rao2014

RAO_CELLS="GM12878_primary+replicate HeLa"

for CELL in ${RAO_CELLS} ; do
    
    # download
    wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${CELL}_HiCCUPS_looplist_with_motifs.txt.gz

    # unzip 
    gunzip Rao2014/GSE63525_${CELL}_HiCCUPS_looplist_with_motifs.txt.gz
done

#=======================================================================
# ChIA-pet data from Tang et al 2015 Cell
#=======================================================================
mkdir -p Tang2015 

wget -P Tang2015 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1872nnn/GSM1872886/suppl/GSM1872886%5FGM12878%5FCTCF%5FPET%5Fclusters.txt.gz
wget -P Tang2015 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1872nnn/GSM1872887/suppl/GSM1872887%5FGM12878%5FRNAPII%5FPET%5Fclusters.txt.gz
wget -P Tang2015 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1872nnn/GSM1872888/suppl/GSM1872888%5FHeLa%5FCTCF%5FPET%5Fclusters.txt.gz
wget -P Tang2015 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1872nnn/GSM1872889/suppl/GSM1872889%5FHeLa%5FRNAPII%5FPET%5Fclusters.txt.gz

# unzip all files
gunzip Tang2015/*.gz

#=======================================================================
# ENCODE batch download 
#=======================================================================
mkdir -p ENCODE
# manully selct for ChIP-seq, TF, bigWig in hg19
# on this site https://www.encodeproject.org/search
# result in this URL: https://www.encodeproject.org/report/?type=Experiment&assay_title=ChIP-seq&assembly=hg19&target.investigated_as=transcription+factor&files.file_type=bigWig&limit=all
# column description of metadata.tsv can be found here: https://www.encodeproject.org/help/batch-download/


# download files.txt with download URLs to all files:
wget -O ENCODE/files.txt https://www.encodeproject.org/batch_download/type%3DExperiment%26assay_title%3DChIP-seq%26assembly%3Dhg19%26target.investigated_as%3Dtranscription%2Bfactor%26files.file_type%3DbigWig%26files.file_type%3Dbam


# download report.tsv files
wget -O ENCODE/report.tsv "https://www.encodeproject.org/report.tsv?type=Experiment&assay_title=ChIP-seq&assembly=hg19&target.investigated_as=transcription+factor&files.file_type=bigWig&files.file_type=bam"

# download only metadata.tsv file (with first link in files.txt)
head -n 1 ENCODE/files.txt | wget -O ENCODE/metadata.tsv -i -

# filter by metadata using R script filter_ENCODE.R
cd ..
Rscript R/filter_ENCODE.R
cd data

# define bash function for download update
function update_download {
  URL=$1
  FILE=ENCODE/Experiments/$(basename $URL)

  if test -e "$FILE"
  then zflag="-z $FILE"
  else zflag=
  fi
  
  # echo $zflag
  curl -o "$FILE" $zflag -L "$URL"  
}
export -f update_download

# download files in parallel
cat ENCODE/URLs.fcDF.txt | parallel -j 10 update_download
cat ENCODE/URLs.fcDF_selectedTF.txt | parallel -j 3 --no-notice update_download
cat ENCODE/URLs.fc_HELA_selected.txt | parallel -j 1 --no-notice update_download

mkdir -p ENCODE/Experiments
cd ENCODE/Experiments 

xargs -P 10 -n 1 curl -O -L < ../URLs.fltOuttype.txt

cd ../..

#=======================================================================
# Download ENCODE bigWig files from UCSC 
#=======================================================================

SELECTED_FIELS="
wgEncodeSydhTfbsGm12878Stat1StdSig.bigWig
wgEncodeSydhTfbsGm12878Stat3IggmusSig.bigWig
wgEncodeSydhTfbsGm12878Yy1StdSig.bigWig
wgEncodeSydhTfbsGm12878Znf143166181apStdSig.bigWig
wgEncodeSydhTfbsGm12878Znf274StdSig.bigWig
wgEncodeSydhTfbsGm12878Znf384hpa004051IggmusSig.bigWig
wgEncodeSydhTfbsGm12878Ctcfsc15914c20StdSig.bigWig
wgEncodeSydhTfbsGm12878Rad21IggrabSig.bigWig
wgEncodeSydhTfbsGm12878Pol2IggmusSig.bigWig
wgEncodeSydhTfbsGm12878Pol2StdSig.bigWig
wgEncodeSydhTfbsGm12878NfkbTnfaIggrabSig.bigWig
"

mkdir -p ENCODE/UCSC
for F in $SELECTED_FIELS; do
    wget -P ENCODE/UCSC http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/${F} 
done

#=======================================================================
# Download other data types 
#=======================================================================

#-------------------------------------------------------------------------------
# Open chromatin (DNase-seq) from ENCODE
#-------------------------------------------------------------------------------
# See https://www.encodeproject.org/experiments/ENCSR000EJD/
# First file is "base overlap signal", second file is "signal"
mkdir -p ENCODE/open_chromatin
wget -P ENCODE/open_chromatin https://www.encodeproject.org/files/ENCFF000SLA/@@download/ENCFF000SLA.bigWig
wget -P ENCODE/open_chromatin https://www.encodeproject.org/files/ENCFF000SLH/@@download/ENCFF000SLH.bigWig

#-------------------------------------------------------------------------------
# ChIP-seq input data in GM12878
#-------------------------------------------------------------------------------
# See: https://www.encodeproject.org/experiments/ENCSR398JTO/

mkdir p ENCODE/input
wget -P ENCODE/input http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878InputStdSig.bigWig

#-------------------------------------------------------------------------------
# ChIP-seq BAM file for Rad21 in GM12878 from ENCODE / UCSC
#-------------------------------------------------------------------------------
# See http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/
mkdir -p ENCODE/bam
wget -P ENCODE/bam http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Rad21IggrabAlnRep1.bam
wget -P ENCODE/bam http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Rad21IggrabAlnRep1.bam.bai

# and STAT1
wget -P ENCODE/bam http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Stat1StdAlnRep1.bam
wget -P ENCODE/bam http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Stat1StdAlnRep1.bam.bai

# and CTCF
wget -P ENCODE/bam http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Ctcfsc15914c20StdAlnRep1.bam
wget -P ENCODE/bam http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Ctcfsc15914c20StdAlnRep1.bam.bai


#=======================================================================
# Run Q pipeline on BAM files
#=======================================================================

BAM_FILES="
ENCODE/bam/wgEncodeSydhTfbsGm12878Rad21IggrabAlnRep1.bam
ENCODE/bam/wgEncodeSydhTfbsGm12878Stat1StdAlnRep1.bam
ENCODE/bam/wgEncodeSydhTfbsGm12878Ctcfsc15914c20StdAlnRep1.bam
"

for BAM in $BAM_FILES ; do

  echo INFO: Working on file $BAM

  # iterate over each chromosome
  cut -f 1 UCSC/hg19.chrom.sizes.real_chroms  | while read CHR
  do
    echo INFO: Working on $CHR
    
    # get qfrags
    ${Q} \
      --treatment-sample ${BAM} \
      --out-prefix ${BAM} \
      -w ${CHR}
    
  done
  
  for OUT_TYPE in "qfrags" "shifted-reads" ; do
    
    echo "INFO output type:" $OUT_TYPE 
    #ls -lht ${BAM}-${OUT_TYPE}-*-chip.bed
  
    # combine all chromosome files
    cat ${BAM}-${OUT_TYPE}-*-chip.bed > ${BAM}-${OUT_TYPE}_allChr_chip.bed
    
    # get bedgraph file
    ${BEDTOOLS} genomecov -bg \
      -i ${BAM}-${OUT_TYPE}_allChr_chip.bed \
      -g UCSC/hg19.chrom.sizes.real_chroms \
      > ${BAM}-${OUT_TYPE}_allChr_chip.bed.bedGraph
  
    # is not case-sensitive sorted at line 2906741.  Please use "sort -k1,1 -k2,2n" with LC_COLLATE=C,  or bedSort and try again.
    export LC_COLLATE=C
    cat ${BAM}-${OUT_TYPE}_allChr_chip.bed.bedGraph \
      | sort -k1,1 -k2,2n \
      > ${BAM}-${OUT_TYPE}_allChr_chip.bed.sorted.bedGraph
    
    # convert bedGraph into BigWig format
    ${BIN}/bedGraphToBigWig \
      ${BAM}-${OUT_TYPE}_allChr_chip.bed.sorted.bedGraph \
      UCSC/hg19.chrom.sizes.real_chroms \
      ${BAM}-${OUT_TYPE}_allChr_chip.bed.sorted.bedGraph.bw

  done  # OUT_TYPE
  
done # BAM

#=======================================================================
# ChIP-nexus data from Tang2015 
#=======================================================================
# GEO https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1872890

# download SRA tool
wget -P ${BIN} https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2-1/sratoolkit.2.8.2-1-centos_linux64.tar.gz
tar xvfz  ${BIN}/sratoolkit.2.8.2-1-centos_linux64.tar.gz -C ${BIN}
SRA_DIR=sratoolkit.2.8.2-1-centos_linux64

#=======================================================================
# get Bowtie
#=======================================================================

wget -P ${BIN} https://kent.dl.sourceforge.net/project/bowtie-bio/bowtie2/2.3.2/bowtie2-2.3.2-linux-x86_64.zip
unzip ${BIN}/bowtie2-2.3.2-linux-x86_64.zip -d ${BIN}

#=======================================================================
# get Samtools
#=======================================================================
wget -P ${BIN} https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2
tar xvfj ${BIN}/samtools-1.5.tar.bz2 -C ${BIN}

cd ${BIN}
BIN_ABS=$(pwd)
cd samtools-1.5
./configure --prefix=${BIN_ABS}
make
make install

#=======================================================================
# get tabix
#=======================================================================
wget -P ${BIN} https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2
tar jxvf ${BIN}/tabix-0.2.6.tar.bz2 -C ${BIN}
cd ${BIN}/tabix-0.2.6
make

#=======================================================================
# get fastqc
#=======================================================================
wget -P ${BIN} https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip ${BIN}/fastqc_v0.11.5.zip -d ${BIN}
chmod ug+rx ${BIN}/FastQC/fastqc

#-------------------------------------------------------------------------------
# get human reference genome:
#-------------------------------------------------------------------------------
mkdir -p hg19
wget -P hg19 ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip
unzip hg19/hg19.zip -d hg19

#-------------------------------------------------------------------------------
# get SRA of RAD21 ChIP-nexus
#-------------------------------------------------------------------------------
./bin/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump -O Tang2015 --gzip SRR2312570
./bin/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump -O Tang2015 --gzip SRR2312571

#=======================================================================
# Process ChIP-nexus data with the Q-nexus pipeline
#=======================================================================

for SAMPLE in "SRR2312570" "SRR2312571" ; do
  
  # 1) flexcat: Adapter Trimming and filtering for fixed Barcode
  ${Q_DIR}/bin/flexcat_noavx2 Tang2015/${SAMPLE}.fastq.gz \
    -tl 5 -tt -t -ml 0 -er 0.2 -ol 4 -app -ss \
    -a ${Q_DIR}/data/adapters.fa \
    -b ${Q_DIR}/data/barcodes.fa \
    -tnum 20 \
    -o Tang2015/${SAMPLE}.flexcat.fastq

  # Mapping:
  ${BIN}/bowtie2-2.3.2/bowtie2 \
    -p 20 \
    -x hg19/hg19 \
    -U Tang2015/${SAMPLE}.flexcat_matched_barcode.fastq \
    -S Tang2015/${SAMPLE}.flexcat_matched_barcode.fastq.bowtie2.hg19.sam

  # 2) remove duplicates
  ${Q_DIR}/bin/nexcat \
    Tang2015/${SAMPLE}.flexcat_matched_barcode.fastq.bowtie2.hg19.sam \
    -fc "(.*)[H|U|M|_]+(.*)"
  
  BAM=Tang2015/${SAMPLE}.flexcat_matched_barcode.fastq.bowtie2.hg19_filtered.bam
  
  # iterate over each chromosome
  cut -f 1 UCSC/hg19.chrom.sizes.real_chroms  | while read CHR
  do
    echo INFO: Working on $CHR
    
    # get qfrags
    ${Q} --nexus-mode -t ${BAM} -o ${BAM} -w ${CHR}
    
  done
  
  # iterate over output types
  for OUT_TYPE in "qfrags" "shifted-reads" ; do
    
    echo "INFO output type:" $OUT_TYPE 

    # combine all chromosome files
    cat ${BAM}-${OUT_TYPE}-*-chip.bed > ${BAM}-${OUT_TYPE}_allChr_chip.bed
    
    # get bedgraph file
    ${BEDTOOLS} genomecov -bg \
      -i ${BAM}-${OUT_TYPE}_allChr_chip.bed \
      -g UCSC/hg19.chrom.sizes.real_chroms \
      > ${BAM}-${OUT_TYPE}_allChr_chip.bed.bedGraph
  
    # is not case-sensitive sorted at line 2906741.  Please use "sort -k1,1 -k2,2n" with LC_COLLATE=C,  or bedSort and try again.
    export LC_COLLATE=C
    cat ${BAM}-${OUT_TYPE}_allChr_chip.bed.bedGraph \
      | sort -k1,1 -k2,2n \
      > ${BAM}-${OUT_TYPE}_allChr_chip.bed.sorted.bedGraph
    
    # convert bedGraph into BigWig format
    ${BIN}/bedGraphToBigWig \
      ${BAM}-${OUT_TYPE}_allChr_chip.bed.sorted.bedGraph \
      UCSC/hg19.chrom.sizes.real_chroms \
      ${BAM}-${OUT_TYPE}_allChr_chip.bed.sorted.bedGraph.bw

  done  # OUT_TYPE

done # SAMPLE

#=======================================================================
# Download predicted loops from Oti et al. 2016 
#=======================================================================
mkdir -p Oti2016

wget -P Oti2016 https://zenodo.org/record/29423/files/ctcf_predictedloops_ENCODE_chipseq_datasets.tar.gz