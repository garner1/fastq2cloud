#!/usr/bin/env bash

fastq=$1			# fastq.gz file as input
zcat $fastq | paste - - - - |cut -f2| split --lines=1000000 - /home/garner1/Work/dataset/fastq2cloud
