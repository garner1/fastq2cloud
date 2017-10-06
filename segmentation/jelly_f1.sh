#!/usr/bin/env bash

kmerlen=$1
fastq=$2
output_file=$3

jellyfish count -C -m $kmerlen -o $output_file -s 100M -t 32 --out-counter-len 4 $fastq
