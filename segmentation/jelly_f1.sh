#!/usr/bin/env bash

kmerlen=$1
fastq=$2
output_file=$3

jellyfish count -C -m $kmerlen -o $output_file -s 5G -t 32 $fastq --disk
