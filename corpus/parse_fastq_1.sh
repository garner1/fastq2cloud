#!/usr/bin/env bash

fastq=$1			# fastq.gz file as input
chunks=$2			# numb of chunks to split file for parallel processing
outPrefix=$3			# /home/garner1/Work/dataset/fastq2cloud/reads/reads_

zcat $fastq | paste - - - - | cut -f2 | split -l$chunks "$fastq"_readonly $outPrefix 
