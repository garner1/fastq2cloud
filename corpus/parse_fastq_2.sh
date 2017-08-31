#!/usr/bin/env bash

fastq=$1			# fastq.gz file as input
chunks=$2			# numb of lines per files
outPrefix=$3			# /home/garner1/Work/dataset/fastq2cloud/reads/reads_

cat $fastq | paste - - - - | cut -f2 | split -l$chunks - $outPrefix 
