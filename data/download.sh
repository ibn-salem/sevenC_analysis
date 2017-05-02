#!/bin/bash

#=======================================================================
# this script is suppost do document all downlaods in the data folder
#=======================================================================

# set some variables here:
BIN=bin
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
xargs -n 1 curl -O -L < ../URLs.fltOuttype.txt
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
# mkdir -p Tang2015 
# wget -P Tang2015 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1872nnn/GSM1872890/suppl/GSM1872890%5FGM12878%5FRAD21%5Fnexus%5FMACE%5Fforward.bw
# wget -P Tang2015 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1872nnn/GSM1872890/suppl/GSM1872890%5FGM12878%5FRAD21%5Fnexus%5FMACE%5Freverse.bw
# wget -P Tang2015 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1872nnn/GSM1872891/suppl/GSM1872891%5FGM12878%5FSMC3%5Fnexus%5FMACE%5Fforward.bw
# wget -P Tang2015 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1872nnn/GSM1872891/suppl/GSM1872891%5FGM12878%5FSMC3%5Fnexus%5FMACE%5Freverse.bw
