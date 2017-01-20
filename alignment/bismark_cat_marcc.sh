#!/bin/bash

#SBATCH
#SBATCH --job-name=bismark-align-ilee
#SBATCH --time=0:5:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#####SBATCH --cpus-per-task=3
#SBATCH --mail-type=end
#SBATCH --mail-user=ilee29@jhu.edu

### loading modules
module load samtools

### execution
dir=/scratch/groups/wtimp1/170119_chicken/aligned/${1}
bamid=${dir}/*${1}*pe.bam
outpath=${dir}/${1}.full.bam

samtools cat -o ${outpath} ${bamid}

echo ${1} " done concatanating"



