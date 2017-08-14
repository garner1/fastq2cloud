#!/usr/bin/env bash

echo "Create MC model ..." 	# approx 4m to run
Rscript /home/garner1/Dropbox/pipelines/fastq2cloud/module_1/MC_model_from_fasta.R

echo "Preapare reads ..."	# approx 14s to run
bash ./module_2/parse_fastq.sh /home/garner1/Work/dataset/bliss/BB23_R1.fastq.gz

echo "Prepare corpus ..."
parallel python module_2/create_sentences.py {} ~/Work/dataset/fastq2cloud/transitionMatrix_3to3.csv ::: ~/Work/dataset/fastq2cloud/reads_a?
