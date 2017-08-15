#!/usr/bin/env bash

fastq=$1			# fastq.gz file as input
zcat $fastq | paste - - - - | cut -f2 | split --lines=100000 - /home/garner1/Work/dataset/fastq2cloud/reads/reads_

# zcat $fastq | paste - - - - | cut -f2 | head -1000 > /home/garner1/Work/dataset/fastq2cloud/reads/reads_test
