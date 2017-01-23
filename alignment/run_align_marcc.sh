#!/bin/bash

outdir=/scratch/groups/wtimp1/170119_chicken/aligned/clnref

for samp in ACTTGA AGTCAA AGTTCC ATGTCA CAGATC CCGTCC CTTGTA GATCAG GGCTAC GTCCGC GTGAAA TAGCTT
do
    mkdir ${outdir}/${samp}

    for lane in {1..8}
    do
	##bismark align
	if [ 1 -eq 1 ]; then	    
	    sbatch bismark_align_marcc.sh ${lane} ${samp} ${outdir}/${samp}
	fi
    done
    
    if [[  1 -eq 0 ]]; then
	sbatch bismark_cat_marcc.sh ${samp}
    else
	echo ${samp} "already concatanated"
    fi
    
    if [ 1 -eq 0 ]; then
	sbatch bismark_extract_marcc.sh $samp
    fi

    if [ 1 -eq 0 ]; then
	sbatch bismark_cytoReport_marcc.sh $samp
    fi

done    





