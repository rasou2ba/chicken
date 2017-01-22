#!/bin/bash

trimpath=/home/isac/Code/trim_galore_zip/trim_galore
bismarkpath=/home/isac/Code/Bismark

refpath=/atium/Data/Reference/chicken/galGal5cln
rawdir=/mithril/Data/NGS/Raw/150415_HiSeqNorwayChicken/wtimp1_118512/FASTQ
outdir=/scratch/groups/wtimp1/170119_chicken/aligned/${2}
tmpdir=/scratch/users/ilee29@jhu.edu/tmp

lanesamp=${1}_${2}

mkdir ${tmpdir}/${lanesamp}
rm ${tmpdir}/${lanesamp}/*

fastq1=`ls ${rawdir}/C6HRUANXX_${lanesamp}_1.fastq.gz`
fastq2=`ls ${rawdir}/C6HRUANXX_${lanesamp}_2.fastq.gz`

${trimpath} --paired ${fastq1} ${fastq2} \
    -o ${tmpdir}/${lanesamp}

trim1=`ls ${tmpdir}/${lanesamp}/*val_1.fq.gz`
trim2=`ls ${tmpdir}/${lanesamp}/*val_2.fq.gz`

echo ${trim1}
echo ${trim2}

${bismarkpath}/bismark --bam --non_directional --bowtie2 \
    -p 4 \
    --genome ${refpath} \
    -1 ${trim1} \
    -2 ${trim2} \
    --temp_dir ${tmpdir}/${lanesamp} \
    --output_dir ${outdir}

rm -R ${tmpdir}/${lanesamp}

