#!/usr/bin/env bash

fastq=$1			# fastq.gz file as input
numblines=$2
outPrefix=$3			# /home/garner1/Work/dataset/fastq2cloud/reads/reads_

zcat $fastq | paste - - - - | cut -f2 | split --lines=$numblines - $outPrefix 
