#!/usr/bin/env bash

fastqgz=$1			# full path to fastq file
kmer_length=$2			# kmer length
dim=$3				# dimension of the MC
output_dir=$4			# output directory

zcat $fastqgz | paste - - - - | cut -f-2 | tr '\t' '\n' > "$output_dir"reads.fasta
jellyfish count -m $kmer_length -o $output_dir -c 3 -s 10000000 -t 32 --both-strands $output_dir/reads.fasta

jellyfish dump "$output_dir"_0 | paste - - | tr -d '>' | 
awk -v dim="$dim" '{OFS="\t";print $1,substr($2,1,dim),substr($2,dim+1,dim)}' |
sort -k2,2 > "$output_dir"_count-kmer-kmer.tsv

jellyfish dump "$output_dir"_0 | paste - - | tr -d '>' | 
awk -v dim="$dim" '{OFS="\t";print $1,substr($2,1,dim)}' |
datamash -s -g 2 sum 1 | sort -k1,1 > "$output_dir"_kmer-sumCount.tsv

join -1 2 -2 1 -o 1.2,1.3,1.1,2.2 "$output_dir"_count-kmer-kmer.tsv "$output_dir"_kmer-sumCount.tsv | 
awk '{OFS="\t"; print $1$2,-log($3/$4)}' > "$output_dir"transitionMatrix_fromFastq_"$dim"to"$dim".csv

rm -f "$output_dir"_count-kmer-kmer.tsv "$output_dir"_kmer-sumCount.tsv "$output_dir"reads.fasta "$output_dir"_0
