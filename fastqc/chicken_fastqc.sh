#!/bin/bash

subsetdir=/mithril/Data/NGS/Aligned/150415_HiSeqChick/subset
outdir=~/Dropbox/Data/Genetics/MethSeq/150415_chicken/fastqc
codedir=~/Code/ilee/util

fastqc -t 8 -o $outdir -f fastq `ls -d $subsetdir/*`

#for samp in `ls $subsetdir`; do
#    echo $samp
#    fastqc -t 8 -o $outdir -f fastq ${subsetdir}$samp
#done    





