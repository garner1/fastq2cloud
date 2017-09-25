#!/usr/bin/env bash

kmerlen=$1
fasta=$2
output_dir=$3

jellyfish count -m $kmerlen -o $output_dir -c 3 -s 10000000 -t 32 --both-strands $fasta
