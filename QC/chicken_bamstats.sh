#!/bin/bash

bamdir=/atium/Data/NGS/Aligned/170120_chicken
outdir=~/Dropbox/Data/Genetics/MethSeq/170120_chicken
codedir=~/Code/ilee/util

rm ${outdir}/alignnum.csv
for bampath in `ls ${bamdir}/*/*full.bam`; do
    echo $bampath
    bamname=${bampath##*/}
    echo $bamname
    num=`samtools view $bampath| wc -l`
    echo ${bamname},${num} >> ${outdir}/alignnum.csv
    
done    



#for samp in `ls $subsetdir`; do
#    echo $samp
#    fastqc -t 8 -o $outdir -f fastq ${subsetdir}$samp
#done    





