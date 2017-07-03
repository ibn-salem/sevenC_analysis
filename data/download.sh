#!/bin/bash

#=======================================================================
# this script is suppost do document all downlaods in the data folder
#=======================================================================

# set some variables here:
BIN=bin
Q_DIR=../../Q
Q=${Q_DIR}/bin/Q
BEDTOOLS=${BIN}/bedtools2/bin/bedtools

mkdir -p ${BIN}

#-----------------------------------------------------------------------
# General genome assembly based data and tools from UCSC:
#-----------------------------------------------------------------------

# UCSC liftover chains
mkdir -p UCSC
wget -P UCSC http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz
wget -P UCSC http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
gunzip UCSC/*.gz

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


#-----------------------------------------------------------------------
# specific data sets
#-----------------------------------------------------------------------

#=======================================================================
# TF motif sites from ENCODE factor book
#=======================================================================
# manually follow this guide: 
#~ You can also access all these tables from the Table Browser:
#~ http://genome.ucsc.edu/cgi-bin/hgTables

#~ 1. Select the hg19 human assembly.
#~ 2. Set "group:" to "All Tables"
#~ 3. From "table:" select the factorbookMotifPos table.
#~ *At this point you could use the filter or intersection tools to limit the
#~ output to factors or locations of interest (via a bed file of coordinates,
#~ see more about the Table Browser here:
#~ http://www.openhelix.com/cgi/tutorialInfo.cgi?id=28)
#~ 4. Set "output format:" to either "custom track" or "BED" and click "get
#~ output".
# SOURCE: http://redmine.soe.ucsc.edu/forum/index.php?t=msg&goto=18147&S=216396035f6ccb52545e9b1627f21599



mkdir -p factorbook

wget -P factorbook http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/factorbookMotifPos.txt.gz

gunzip factorbook/factorbookMotifPos.txt.gz

# download as factorbook/factorbookMotifPos.bed
cat  factorbook/factorbookMotifPos.txt \
	| awk -F"\t" '{if ($5 == "CTCF") print }' OFS="\t" \
	| cut -f2-7 \
	> factorbook/factorbookMotifPos.txt.CTCF.bed


#=======================================================================
# Hi-C data from Rao et al 2014 Cell
#=======================================================================
mkdir -p Rao2014


RAO_CELLS="GM12878_primary+replicate HMEC HUVEC HeLa IMR90 K562 KBM7 NHEK"

for CELL in ${RAO_CELLS} ; do
    
    # download
    #wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${CELL}_Arrowhead_domainlist.txt.gz
    #wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${CELL}_HiCCUPS_looplist.txt.gz
    wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${CELL}_HiCCUPS_looplist_with_motifs.txt.gz

    # unzip 
    # gunzip Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt.gz
    # gunzip Rao2014/GSE63525_${CELL}_HiCCUPS_looplist.txt.gz
    gunzip Rao2014/GSE63525_${CELL}_HiCCUPS_looplist_with_motifs.txt.gz
        
    # # re-format TADs into bed file
    # tail -n +2 Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt \
    #   |cut -f 1-3 \
    #   | sed -e 's/^/chr/' \
    #   > Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt.bed 
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
# Capture Hi-C data from Mifsud2015
#=======================================================================
mkdir -p Mifsud2015
wget -P Mifsud2015 http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2323/E-MTAB-2323.additional.1.zip
unzip Mifsud2015/E-MTAB-2323.additional.1.zip -d Mifsud2015

#=======================================================================
# ENCODE batch download 
#=======================================================================
mkdir -p ENCODE
# manully selct for ChIP-seq, TF, bigWig in hg19
# on this site https://www.encodeproject.org/search
# result in this URL: https://www.encodeproject.org/report/?type=Experiment&assay_title=ChIP-seq&assembly=hg19&target.investigated_as=transcription+factor&files.file_type=bigWig&limit=all
# column description of metadata.tsv can be found here: https://www.encodeproject.org/help/batch-download/


# download files.txt with download URLs to all files:
wget -O ENCODE/files.txt https://www.encodeproject.org/batch_download/type%3DExperiment%26assay_title%3DChIP-seq%26assembly%3Dhg19%26target.investigated_as%3Dtranscription%2Bfactor%26files.file_type%3DbigWig

# download report.tsv files
wget -O ENCODE/report.tsv "https://www.encodeproject.org/report.tsv?type=Experiment&assay_title=ChIP-seq&assembly=hg19&target.investigated_as=transcription+factor&files.file_type=bigWig"

# download only metadata.tsv file (with first link in files.txt)
head -n 1 ENCODE/files.txt \
| wget -O ENCODE/metadata.tsv -i -

# filter by metadata using R scirpt
# cd ..
# Rscript R/filter_ENCODE.R
# cd data

mkdir -p ENCODE/Experiments
cd ENCODE/Experiments 
xargs -n 1 curl -O -L < ../URLs.flt.txt
xargs -P 10 -n 1 curl -O -L < ../URLs.fltOuttype.txt
xargs -n 1 curl -O -L < ../URLs.fcDF.txt

cd ../..


#=======================================================================
# ENCODE from UCSC 
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
# Other data types 
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
# wget -P ENCODE/input https://www.encodeproject.org/files/ENCFF941TZZ/@@download/ENCFF941TZZ.bam
# samtools index ENCODE/input/ENCFF941TZZ.bam

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


# get chromosome sizes
wget -P UCSC http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
head -n 24 UCSC/hg19.chrom.sizes > UCSC/hg19.chrom.sizes.real_chroms

BAM_FILES="
ENCODE/bam/wgEncodeSydhTfbsGm12878Rad21IggrabAlnRep1.bam
ENCODE/bam/wgEncodeSydhTfbsGm12878Stat1StdAlnRep1.bam
ENCODE/bam/wgEncodeSydhTfbsGm12878Ctcfsc15914c20StdAlnRep1.bam
"

# BAM="ENCODE/bam/wgEncodeSydhTfbsGm12878Rad21IggrabAlnRep1.bam"

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


#-------------------------------------------------------------------------------
# get human reference genome:
#-------------------------------------------------------------------------------
# wget -P UCSC http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
# ${BIN}/twoBitToFa UCSC/hg19.2bit  UCSC/hg19.fa
# 
# # build index for mapping with bowtie
# ./${BIN}/bowtie-1.2.1.1/bowtie-build UCSC/hg19.fa UCSC/hg19.fa

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
  
# can be used as:
# ${BIN}/bowtie-1.2.1.1/bowtie

#=======================================================================
# DNA shape from GBshape
#=======================================================================
# mkdir -p GBshape
# wget -P GBshape ftp://rohslab.usc.edu/hg19/

#=======================================================================
# OLD PARTS (do not run)
#=======================================================================

#=======================================================================
# JASPAR vertebrate core data base (Downloaded 26.05.15)
#=======================================================================
#~ mkdir -p JASPAR
#~ wget -P JASPAR http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt
#~ # convert matrix to transfac format using RASAT convert-matrix
#~ convert-matrix -i JASPAR/pfm_vertebrates.txt -o JASPAR/pfm_vertebrates.txt.tf -from jaspar -to transfac
#~ 


#=======================================================================
# ENCODE data from UCSC Genome Browser
#=======================================================================
# see https://genome.ucsc.edu/ENCODE/downloads.html
# see https://www.encodeproject.org/search/

#~ ENCODE_LINKS="
#~ http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsGm12878Rad21V0416101RawRep1.bigWig
#~ http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Znf143166181apStdSig.bigWig
#~ http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Stat3IggmusSig.bigWig
#~ "
#~ ALL_TF_LINKS="
#~ wgEncodeSydhTfbsGm12878Bhlhe40cIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Brca1a300IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Cdpsc6327IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878CfosStdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Chd1a301218aIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Chd2ab68301IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Corestsc30189IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Ctcfsc15914c20StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878E2f4IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Ebf1sc137065StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Elk112771IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878ErraIggrabSig.bigWig
#~ wgEncodeSydhTfbsGm12878Gcn5StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Ikzf1iknuclaStdSig.bigWig
#~ wgEncodeSydhTfbsGm12878InputIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878InputIggrabSig.bigWig
#~ wgEncodeSydhTfbsGm12878InputStdSig.bigWig
#~ wgEncodeSydhTfbsGm12878InputTnfaIggrabSig.bigWig
#~ wgEncodeSydhTfbsGm12878Irf3IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878JundIggrabSig.bigWig
#~ wgEncodeSydhTfbsGm12878JundStdSig.bigWig
#~ wgEncodeSydhTfbsGm12878MafkIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878MaxIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878MaxStdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Mazab85725IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Mxi1IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Nfe2sc22827StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878NfkbTnfaIggrabSig.bigWig
#~ wgEncodeSydhTfbsGm12878NfyaIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878NfybIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Nrf1IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878P300IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878P300bStdSig.bigWig
#~ wgEncodeSydhTfbsGm12878P300sc584IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Pol2IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Pol2StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Pol2s2IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Pol3StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Rad21IggrabSig.bigWig
#~ wgEncodeSydhTfbsGm12878Rfx5200401194IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Sin3anb6001263IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Smc3ab9263IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Spt20StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Srebp1IggrabSig.bigWig
#~ wgEncodeSydhTfbsGm12878Srebp2IggrabSig.bigWig
#~ wgEncodeSydhTfbsGm12878Stat1StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Stat3IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Tblr1ab24550IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878TbpIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Tr4StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Usf2IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878WhipIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Yy1StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Znf143166181apStdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Znf274StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Znf384hpa004051IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Zzz3StdSig.bigWig
#~ "
#~ 
SELECTED_ENCODE="
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

mkdir -p ENCODE
for F in $SELECTED_ENCODE; do
    wget -P ENCODE http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/${F} 
done


#=======================================================================
# ENCODE TF ChIP-seq data:
#=======================================================================
# metadata file from ENCODE : metadata.csv
wget "https://www.encodeproject.org/metadata/type=Experiment&assay_term_name=ChIP-seq&status=released&assembly=hg19&target.investigated_as=transcription+factor&files.file_type=bed+narrowPeak&replication_type=isogenic/metadata.tsv"

# filter for CTCF ChIP-seq peaks
cat metadata.tsv \
    | grep "optimal idr thresholded peaks" \
    | grep "CTCF-human" \
    > metadata.tsv.idr.CTCF
# add header
head -n 1 metadata.tsv > metadata.tsv.header
cat  metadata.tsv.header  metadata.tsv.idr.CTCF > metadata.tsv.idr.CTCF.withHeader

cat metadata.tsv.idr.CTCF \
    |awk -F'\t' '{print $39}' \
    > metadata.tsv.idr.CTCF.links

# download all files
mkdir -p ENCODE
wget -P ENCODE -i metadata.tsv.idr.CTCF.links

# unzip all files
gunzip ENCODE/*.bed.gz

# 
# #=======================================================================
# # Tang2015 ChIP-nexus data for Rad21 and SMC3:
# #=======================================================================
mkdir -p Tang2015 
# wget -P Tang2015 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1872nnn/GSM1872890/suppl/GSM1872890%5FGM12878%5FRAD21%5Fnexus%5FMACE%5Fforward.bw
# wget -P Tang2015 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1872nnn/GSM1872890/suppl/GSM1872890%5FGM12878%5FRAD21%5Fnexus%5FMACE%5Freverse.bw
# wget -P Tang2015 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1872nnn/GSM1872891/suppl/GSM1872891%5FGM12878%5FSMC3%5Fnexus%5FMACE%5Fforward.bw
# wget -P Tang2015 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1872nnn/GSM1872891/suppl/GSM1872891%5FGM12878%5FSMC3%5Fnexus%5FMACE%5Freverse.bw


