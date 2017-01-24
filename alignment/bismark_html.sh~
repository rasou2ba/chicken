#!/bin/bash

wdir=$PWD
outdir=${2}
bismarkpath=~/Code/Bismark
refpath=/atium/Data/Reference/chicken/galGal5cln

${bismarkpath}/bismark2bedGraph --dir ${outdir}/ -o ${1}.full.bedGraph ${outdir}/CpG*${1}*.txt.gz

${bismarkpath}/coverage2cytosine --gzip --dir ${outdir}/ --genome_folder ${refpath} -o ${1}.cyto.txt.gz ${outdir}/${1}.full.bismark.cov.gz

cd $outdir

${bismarkpath}/bismark2report 

cd $wdir
