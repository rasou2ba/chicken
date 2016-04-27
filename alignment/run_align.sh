#!/bin/bash

for samp in ACTTGA AGTCAA AGTTCC ATGTCA CAGATC CCGTCC CTTGTA GATCAG GGCTAC GTCCGC GTGAAA TAGCTT
do
    mkdir chickenaligned/${samp}

    for lane in {1..8}
    do
	##bismark align
	if [ 0 -eq 0 ]; then	    
	    qsub -pe orte 16 -N bismarkalign_${samp} -cwd bismark_align.sh ${lane} ${samp}
	fi
    done

    if [ 0 -eq 0 ]; then
	qsub -pe orte 16 -hold_jid bismarkalign_${samp} -N bismarkcat_${samp} -cwd bismark_cat.sh $samp
    fi
    
    if [ 0 -eq 0 ]; then
	    qsub -pe orte 16 -hold_jid bismarkcat_${samp} -N bismarkextract_${samp} -cwd bismark_extract.sh $samp
    fi

done    





