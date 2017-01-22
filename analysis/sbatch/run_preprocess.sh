#!/bin/bash

outdir=/atium/Data/NGS/Aligned/170120_chicken

for samp in ACTTGA #AGTCAA AGTTCC ATGTCA CAGATC CCGTCC CTTGTA GATCAG GGCTAC GTCCGC GTGAAA TAGCTT
do
    sbatch preprocess_batch.sh ${samp}
done    





