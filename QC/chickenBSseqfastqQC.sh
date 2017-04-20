#!/bin/bash
# wrapper for chicken data fastq QC

qcpy=~/Code/ilee/util/fastqQC.py
qcR=~/Code/ilee/util/fastqQC.R
fastqpath=/atium/Data/NGS/Aligned/170120_chicken/QC/fastqfiltered
rawpath=/mithril/Data/NGS/Raw/150415_HiSeqNorwayChicken/wtimp1_118512/FASTQ
qcpath=~/Dropbox/Data/Genetics/MethSeq/170120_chicken/QC
for lane in {1..8}
do
  samp=${lane}_CCGTCC
  read1=`ls ${fastqpath}/${samp}*val_1.fq.gz`
  read2=`ls ${fastqpath}/${samp}*val_2.fq.gz`
  base=$(basename "$read1" _sub_1_val_1.fq.gz)
  raw1=`ls ${rawpath}/*${base}_1*`
  raw2=`ls ${rawpath}/*${base}_2*`
  
  path1=~/Dropbox/Data/Genetics/MethSeq/170120_chicken/QC/${samp}_filtered_read1_qualcount.txt
  path2=~/Dropbox/Data/Genetics/MethSeq/170120_chicken/QC/${samp}_filtered_read2_qualcount.txt
  rawpath1=~/Dropbox/Data/Genetics/MethSeq/170120_chicken/QC/${samp}_read1_qualcount.txt
  rawpath2=~/Dropbox/Data/Genetics/MethSeq/170120_chicken/QC/${samp}_read2_qualcount.txt
  
  rawtag=${samp}_raw
  outtag=${samp}_filtered
  subnum=10000
  
#  $qcpy $read1 -s $subnum -o $path1
#  $qcpy $read2 -s $subnum -o $path2
#  $qcpy $raw1 -s $subnum -o $rawpath1
#  $qcpy $raw2 -s $subnum -o $rawpath2
#  
#  $qcR --input=${path1}.gz --input2=${path2}.gz --out=$outtag 
#  $qcR --input=${rawpath1}.gz --input2=${rawpath2}.gz --out=$rawtag 
done
echo all
$qcR --input=$qcpath --pattern=qualcount.txt.gz --out=QC --multiple=1
