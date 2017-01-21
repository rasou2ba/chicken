#!/bin/bash

wdir=$PWD
outdir=/atium/Data/NGS/Aligned/170120_chicken/${1}
bismarkpath=~/Code/Bismark
refpath=/atium/Data/Reference/chicken/galGal5

#${bismarkpath}/bismark_methylation_extractor -p --multicore 8 --gzip \
#    --genome_folder ${refpath} \
#    ${outdir}/${1}.full.bam -o ${outdir} --no_header

#${bismarkpath}/bismark2bedGraph --dir ${outdir}/ --scaffolds -o ${1}.full.bedGraph ${outdir}/CpG*${1}*.txt.gz 
${bismarkpath}/coverage2cytosine --gzip --dir ${outdir}/ --genome_folder ${refpath} -o ${1}.cyto.txt.gz ${outdir}/${1}.full.bismark.cov.gz

cd $outdir

${bismarkpath}/bismark2report 

cd $wdir
