#!/usr/bin/env bash

echo `cat /home/garner1/Work/dataset/fastq2cloud/corpus/reads_test_sentences.txt_counter.txt | tr -d '(),' | cut -d' ' -f-2 | tr ' ' '\n' | LC_ALL=C sort | LC_ALL=C uniq | wc -l`
