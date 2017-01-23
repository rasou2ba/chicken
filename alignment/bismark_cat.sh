#!/bin/bash

dir=${2}
bamid=${dir}/*${1}*pe.bam
outpath=${dir}/${1}.full.bam

samtools cat -o ${outpath} ${bamid}

echo $1 " done concatanating"



