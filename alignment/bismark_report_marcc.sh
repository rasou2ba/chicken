#!/bin/bash

#SBATCH
#SBATCH --job-name=bismarkReport
#SBATCH --time=5:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#####SBATCH --cpus-per-task=3
#SBATCH --mail-type=end
#SBATCH --mail-user=ilee29@jhu.edu

###load module
module load perl

###execute
wdir=$PWD
outdir=${2}
bismarkpath=/home-2/ilee29@jhu.edu/Code/Bismark
refpath=/scratch/groups/wtimp1/Reference/chicken/galGal5cln

${bismarkpath}/bismark2bedGraph --dir ${outdir}/ -o ${1}.full.bedGraph ${outdir}/CpG*${1}*.txt.gz

${bismarkpath}/coverage2cytosine --gzip --dir ${outdir}/ --genome_folder ${refpath} -o ${1}.cyto.txt.gz ${outdir}/${1}.full.bismark.cov.gz

cd $outdir

${bismarkpath}/bismark2report 

cd $wdir
