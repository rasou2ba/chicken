#!/bin/bash

outdir=/atium/Data/NGS/Aligned/170120_chicken

for samp in ACTTGAsubset
do
    mkdir ${outdir}/${samp}

    for lane in 1
    do
	##bismark align
	if [ 0 -eq 0 ]; then
	    ./bismark_align.sh ${lane} ${samp}
	fi
    done

    if [ 0 -eq 1 ]; then
	qsub -pe orte 16 -hold_jid bismarkalign_${samp} -N bismarkcat_${samp} -cwd bismark_cat.sh $samp
    fi
    
    if [ 0 -eq 1 ]; then
	    qsub -pe orte 16 -hold_jid bismarkcat_${samp} -N bismarkextract_${samp} -cwd bismark_extract.sh $samp
    fi

done    





