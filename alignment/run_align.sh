#!/bin/bash

outdir=/atium/Data/NGS/Aligned/170120_chicken/clnref

for samp in ACTTGA #AGTCAA AGTTCC ATGTCA CAGATC CCGTCC CTTGTA GATCAG GGCTAC GTCCGC GTGAAA TAGCTT
do
    mkdir ${outdir}/${samp}

    for lane in 1 #.8}
    do
	##bismark align
	if [ 0 -eq 1 ]; then	    
	    ./bismark_align.sh ${lane} ${samp} ${outdir}/${samp}
	fi
    done
    
    if [[  1 -eq 0 ]]; then
	./bismark_cat.sh ${samp} ${outdir}/${samp}
    else
	echo $samp "already concatanated"
    fi
    
    if [ 0 -eq 0 ]; then
	./bismark_extract.sh ${samp} ${outdir}/${samp}
    fi

done    





