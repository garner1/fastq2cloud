#!/usr/bin/env bash

export LC_ALL=C

jf_file=$1
dim=$2
alpha=$3			# alpha estimator in [0,1] 

# THE CANONICAL REPRESENTATION OF JELLYFISH NEEDS TO INTRODUCE THE REVERSE COMPLEMENT OF EACH COUNTED KMER!!!

# LEFT TO RIGHT PROCESSING
jellyfish dump "$jf_file" | paste - - | tr -d '>' | 
awk -v dim="$dim" '{OFS="\t";print $1,substr($2,1,dim-1),substr($2,dim,1)}' |
sort -k2,2 > "$jf_file"_count-kmer-kmer.tsv

jellyfish dump "$jf_file" | paste - - | tr -d '>' | 
awk -v dim="$dim" '{OFS="\t";print $1,substr($2,1,dim-1)}' |
datamash -s -g 2 sum 1 | sort -k1,1 > "$jf_file"_kmer-sumCount.tsv

join -1 2 -2 1 -o 1.2,1.3,1.1,2.2 "$jf_file"_count-kmer-kmer.tsv "$jf_file"_kmer-sumCount.tsv |
awk -v alpha=$alpha '{OFS="\t"; print $1$2,-log(($3+alpha)/($4+alpha*4))}' > "$jf_file"transitionMatrix_fromFastq_"$dim"to1.csv # use the alpha-estimator

## RIGHT TO LEFT PROCESSING##!!!!TBD!!!
