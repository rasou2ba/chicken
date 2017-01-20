#!/bin/bash

rawdir=/scratch/groups/wtimp1/170119_chicken/fastq

outdir=/scratch/groups/wtimp1/170119_chicken/aligned/${2}

tmpdir=/scratch/tmp

lanesamp=${1}_${2}

mkdir ${tmpdir}/${lanesamp}

rm ${tmpdir}/${lanesamp}/*

~/trim_galore_zip/trim_galore --paired ${rawdir}/C6HRUANXX_${lanesamp}_1.fastq.gz ${rawdir}/C6HRUANXX_${lanesamp}_2.fastq.gz \
    -o ${tmpdir}/${lanesamp}

~/bismark_v0.14.2/bismark --bam --non_directional --bowtie2 \
    -p 4 \
    /home/sgeadmin/ref/Reference/chicken/ \
    -1 ${tmpdir}/${lanesamp}/*val_1.fq.gz \
    -2 ${tmpdir}/${lanesamp}/*val_2.fq.gz \
    --o ${outdir}

rm -R ${tmpdir}/${lanesamp}
