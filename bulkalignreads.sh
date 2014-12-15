#!/bin/bash

setdir="run68mousernaseqtregsplnlung/"
basedir="/home/trevor/jen/"
readsub="data/"
outsub="aligned/"
readsdir=$basedir$readsub$setdir
outdir=$basedir$outsub$setdir

im1_files=($(ls "$readsdir" | grep R1))
im2_files=($(ls "$readsdir" | grep R2))

for ((i=0;i<${#im1_files[@]};i++)); do
	readpathone=$readsdir${im1_files[i]}
	readpathtwo=$readsdir${im2_files[i]}
	baseword=${im1_files[i]}
	outcutone=(${baseword//./ })
	basecut=${outcutone[0]}
	outcuttwo=(${basecut//_/ })
	baseout=${outcuttwo[0]}
	echo $baseout
	outpath=$outdir$baseout
	/home/trevor/Programs/STAR --runThreadN 10 --genomeDir /home/trevor/jen/data/starindex/ENSEMBL.mus_musculus.release-75 --readFilesIn $readpathone $readpathtwo --outFileNamePrefix $outpath --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000
done
