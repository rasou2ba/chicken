#!/bin/bash

outdir=/scratch/groups/wtimp1/170119_chicken/aligned

for samp in ACTTGA #AGTCAA AGTTCC ATGTCA CAGATC CCGTCC CTTGTA GATCAG GGCTAC GTCCGC GTGAAA TAGCTT
do
    mkdir ${outdir}/${samp}

    for lane in {1..8}
    do
	##bismark align
	if [ 0 -eq 1 ]; then	    
	    sbatch bismark_align_marcc.sh ${lane} ${samp}
	fi
    done
    
    if [ ! -e "${outdir}/${samp}/${1}.full.bam" ] && [ 0 -eq 0 ]; then
	sbatch bismark_cat_marcc.sh ${samp}
    else
	echo ${samp} "already concatanated"
    fi
    
    if [ 0 -eq 0 ]; then
	sbatch bismark_extract.sh $samp
    fi

done    





