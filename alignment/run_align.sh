#!/bin/bash

outdir=/scratch/groups/wtimp1/170119_chicken/aligned

for samp in ACTTGA AGTCAA AGTTCC ATGTCA CAGATC CCGTCC CTTGTA GATCAG GGCTAC GTCCGC GTGAAA TAGCTT
do
    mkdir ${outdir}/${samp}

    for lane in {1..8}
    do
	##bismark align
	if [ 0 -eq 1 ]; then	    
	    sbatch bismark_align_marcc.sh ${lane} ${samp}
	fi
    done
    
    if [ !-e ${outdir}/${samp}/${1}.full.bam ]; then
	echo "${1} is done"
	#sbatch bismark_cat.sh $samp
    fi
    
    if [ 0 -eq 1 ]; thenz
	    qsub -pe orte 16 -hold_jid bismarkcat_${samp} -N bismarkextract_${samp} -cwd bismark_extract.sh $samp
    fi

done    





