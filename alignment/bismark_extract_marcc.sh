#!/bin/bash

#SBATCH
#SBATCH --job-name=bismark-align-ilee
#SBATCH --time=5:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#####SBATCH --cpus-per-task=3
#SBATCH --mail-type=end
#SBATCH --mail-user=ilee29@jhu.edu

outdir=/scratch/groups/wtimp1/170119_chicken/aligned/${1}
bismarkpath=/home-2/ilee29@jhu.edu/Code/bismark_v0.17.0
refpath=/scratch/groups/wtimp1/Reference/chicken/galGal5

${bismarkpath}/bismark_methylation_extractor -p --multicore 8 --gzip \
    --genome_folder ${refpath} \
    ${outdir}/${1}.full.bam -o ${outdir} --no_header

${bismarkpath}/bismark2bedGraph --dir ${outdir}/ --counts -o ${1}.full.bedGraph CpG*${1}*.txt.gz 

cd $outdir

${bismarkpath}/bismark2report

