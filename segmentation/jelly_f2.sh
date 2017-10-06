#!/usr/bin/env bash

export LC_ALL=C

jf_file=$1
dim=$2
alpha=$3			# alpha estimator in [0,1] 

# THE CANONICAL REPRESENTATION OF JELLYFISH NEEDS TO INTRODUCE THE REVERSE COMPLEMENT OF EACH COUNTED KMER!!!

jellyfish dump "$jf_file" | paste - - | tr -d '>' > "$jf_file"_FW.fa
jellyfish dump "$jf_file" | perl -lpe '$_ = reverse $_ unless $.%2' - | tr ATCG TAGC | paste - - | tr -d '>' > "$jf_file"_RC.fa

cat "$jf_file"_FW.fa "$jf_file"_RC.fa | datamash -s -g 2 mean 1 > "$jf_file".fa

awk -v dim="$dim" '{OFS="\t";print $2,substr($1,1,dim-1),substr($1,dim,1)}' "$jf_file".fa |
sort -k2,2 > "$jf_file"_count-kmer-kmer.tsv 

awk -v dim="$dim" '{OFS="\t";print $2,substr($1,1,dim-1)}' "$jf_file".fa |
datamash -s -g 2 sum 1 > "$jf_file"_kmer-sumCount.tsv

join -1 2 -2 1 -o 1.2,1.3,1.1,2.2 "$jf_file"_count-kmer-kmer.tsv "$jf_file"_kmer-sumCount.tsv |
awk -v alpha=$alpha '{OFS="\t"; print $1$2,-log(($3+alpha)/($4+alpha*4))}' > "$jf_file".csv # use the alpha-estimator
