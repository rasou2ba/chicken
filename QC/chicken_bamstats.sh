#!/bin/bash

bamdir=/mithril/Data/NGS/Aligned/150415_HiSeqChick
outdir=~/Dropbox/Data/Genetics/MethSeq/150415_chicken/QC
codedir=~/Code/ilee/util

for bampath in `ls ${bamdir}/*/*full.bam`; do
    #echo $bampath
    bamname=${bampath##*/}
    #echo $bamname
    num=`samtools view $bampath| wc -l`
    echo ${bamname},${num} >> ${outdir}/alignnum.csv
    
done    



#for samp in `ls $subsetdir`; do
#    echo $samp
#    fastqc -t 8 -o $outdir -f fastq ${subsetdir}$samp
#done    





