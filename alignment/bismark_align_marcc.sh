#!/bin/bash

#SBATCH
#SBATCH --job-name=bismark-align-ilee
#SBATCH --time=5:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#####SBATCH --cpus-per-task=3
#SBATCH --mail-type=end
#SBATCH --mail-user=ilee29@jhu.edu

trimpath=/home-2/ilee29@jhu.edu/Code/trim_galore_v0.4.2/trim_galore
bismarkpath=/home-2/ilee29@jhu.edu/Code/bismark_v0.17.0

refpath=/scratch/groups/wtimp1/Reference/chicken/galGal5
rawdir=/scratch/groups/wtimp1/170119_chicken/fastq
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

