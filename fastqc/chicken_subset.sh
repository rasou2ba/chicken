#!/bin/bash

fastqdir=/mithril/Data/NGS/Raw/150415_HiSeqNorwayChicken/wtimp1_118512/FASTQ
subsetdir=/mithril/Data/NGS/Aligned/150415_HiSeqChick/subset
codedir=~/Code/ilee/util

for samp in `ls $fastqdir`; do
    ${codedir}/subset.sh $fastqdir $samp $subsetdir    
done    





