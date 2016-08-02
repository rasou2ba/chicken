#!/bin/bash

fastqdir=/mithril/Data/NGS/Raw/150415_HiSeqNorwayChicken/wtimp1_118512/FASTQ
outdir=~/Dropbox/Data/Genetics/MethSeq/150415_chicken/fastqc
codedir=~/Code/ilee/util

for samp in `ls ${fastqdir}/*1.fastq.gz`; do
    lines=`gunzip -c $samp | wc -l`
    numread=`expr $lines / 4`
    fname=${samp##*/}
    echo ${fname},${numread} >> ${outdir}/readnum.csv
    
done    



#for samp in `ls $subsetdir`; do
#    echo $samp
#    fastqc -t 8 -o $outdir -f fastq ${subsetdir}$samp
#done    





