#!/usr/bin/env bash

kmerlen=$1
fasta=$2
output_dir=$3

# consider both strands in the counting otherwise you will not find all matches in the model when searching for kmers and assigning information values
jellyfish count -m $kmerlen -o $output_dir -c 3 -s 10000000 -t 32 $fasta
