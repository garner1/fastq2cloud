#!/usr/bin/env bash

fastq=$1			# fastq.gz file as input
lines=$2			# numb of liness to split file for parallel processing
outPrefix=$3			# /home/garner1/Work/dataset/fastq2cloud/reads/reads_

zcat $fastq | paste - - - - | cut -f2 | split -l$lines - $outPrefix 
