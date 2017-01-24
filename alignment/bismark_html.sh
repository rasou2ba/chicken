#!/bin/bash

wdir=$PWD
outdir=${2}
bismarkpath=~/Code/Bismark

cd $outdir

${bismarkpath}/bismark2report --splitting_report ${2}/${1}.full_splitting_report.txt \
	      --mbias_report ${2}/${1}.full.M-bias.txt

cd $wdir
