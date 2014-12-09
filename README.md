Basic RNA-Seq Analysis Pipeline
===========================

A straightforward but powerful RNA-seq analysis pipeline. Processes files in bulk based on a simple folder based organization system. 

Inputs: 

* .fastq files from sequencing machine
* STAR genome build for genome reads will be aligned to
* .gff file for genome reads will be aligned to
* .txt file with desired differential expression comparisons

Outputs: 

* .sam file of aligned reads
* .txt file of read counts
* .txt file of differentially expressed genes

Version 1.0 of this pipeline is written for mouse data, other organisms and a flag for choosing organism will be added as time permits or they are requested. 
