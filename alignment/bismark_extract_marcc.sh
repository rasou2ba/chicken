#!/bin/bash

outdir=/home/sgeadmin/data/chickenaligned/${1}

~/bismark_v0.14.2/bismark_methylation_extractor -p --multicore 4 --gzip \
    --genome_folder /home/sgeadmin/ref/Reference/human/ \
    ${outdir}/${1}.full.bam -o ${outdir} --no_header

~/bismark_v0.14.2/bismark2bedGraph --dir ${outdir}/ --counts -o ${1}.full.bedGraph ${1}.CpG*.txt.gz 

cd $outdir

~/bismark_v0.14.2/bismark2report

cd /home/sgeadmin/data