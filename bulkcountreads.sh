#!/bin/bash

setdir="run68mousernaseqtregsplnlung/"
basedir="/home/trevor/jen/"
alignedsub="aligned/"
outsub="counts/"
tailsub=".txt"
aligneddir=$basedir$alignedsub$setdir
outdir=$basedir$outsub$setdir

im1_files=($(ls "$aligneddir" | grep Aligned.out.sam))

for ((i=0;i<${#im1_files[@]};i++)); do
	alignedpath=$aligneddir${im1_files[i]}
	baseword=${im1_files[i]}
	outcutone=(${baseword//Aligned/ })
	baseout=${outcutone[0]}
	echo $baseout
	outpath=$outdir$baseout$tailsub
	htseq-count $alignedpath /home/trevor/jen/data/starindex/ENSEMBL.mus_musculus.release-75/Mus_musculus.GRCm38.75.gtf > $outpath
done
