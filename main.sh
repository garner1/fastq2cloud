#!/usr/bin/env bash

# echo "Create MC model ..." 	# approx 4m to run
# Rscript /home/garner1/Dropbox/pipelines/fastq2cloud/module_1/MC_model_from_fasta.R

# # GENERATE FASTQ FILES USING NEAT OR VARSIM_RUN                                                                                                                                              
# # !!! check directories and parameters !!!!                                                                                                                                                 # # parallel "python ~/tools/neat-genreads/genReads.py -r ~/igv/genomes/hg19.fasta -R 50 -o ~/Work/dataset/genome_segmentation/simulated_data -c 1 --job {} 32" ::: `seq 32`
# # python ~/tools/neat-genreads/mergeJobs.py -i  /home/garner1/tool_test/neat -o  /home/garner1/tool_test -s ~/tools/samtools-1.2

# # FASTQ.gz TO READS
# echo "Preapare reads ..."	# approx 14s to run
# bash ./module_2/parse_fastq.sh /home/garner1/Work/dataset/bliss/BB23_R1.fastq.gz
# echo "Done"

# # READS to CORPUS
# # These are pickle files TO BE USED IN GENSIM WORD2VEC sequentially, do not concatenate with cat
# echo "Prepare corpus ..."	# approx 2h with BB23_R1.fastq.gz
# parallel python module_2/create_corpus.py {} ::: ~/Work/dataset/fastq2cloud/reads/reads_*
# # python module_2/create_corpus.py ~/Work/dataset/fastq2cloud/reads/reads_test
# echo "Done"

# echo "Move corpus into specific directory"
# mkdir -p ~/Work/dataset/fastq2cloud/corpus
# mv ~/Work/dataset/fastq2cloud/reads/*_sentences.txt ~/Work/dataset/fastq2cloud/corpus
# echo "Done"

# # echo 'Run gensim word2vec implementation ...'
# # python ./module_3/word2vector.py ~/Work/dataset/fastq2cloud/corpus
# # echo 'Done!

echo 'Parse corpus ...'
cat ~/Work/dataset/fastq2cloud/corpus/reads_test_sentences.txt |
grep "'" | tr -d '^a' | tr -d '^S' | LC_ALL=C sort | LC_ALL=C uniq -c | cat -n | awk '{OFS="\t";print$1,$2,$3}'|LC_ALL=C sort -k3,3 > ~/Work/dataset/fastq2cloud/corpus/corpus__ind-counts-word.tsv
dim=`cat ~/Work/dataset/fastq2cloud/corpus/corpus__ind-counts-word.tsv | wc -l`
echo 'Done'

# echo 'Count pairs ...'
# python ./module_3/count_pairs.py ~/Work/dataset/fastq2cloud/corpus/reads_test_sentences.txt
# echo 'Done'

echo 'Parse cooccurrences ...'
cat ~/Work/dataset/fastq2cloud/corpus/reads_test_sentences.txt_counter.txt | tr -d '(),' | tr ' ' '\t'|LC_ALL=C sort -k1,1  > ~/Work/dataset/fastq2cloud/word1-word2-count
LC_ALL=C join -1 1 -2 3 -o 1.1 1.2 1.3 2.1 ~/Work/dataset/fastq2cloud/word1-word2-count ~/Work/dataset/fastq2cloud/corpus/corpus__ind-counts-word.tsv|
tr ' ' '\t'|LC_ALL=C sort -k2,2 > ~/Work/dataset/fastq2cloud/word1-word2-count-ind1
LC_ALL=C join -1 2 -2 3 -o 1.1 1.2 1.3 1.4 2.1 ~/Work/dataset/fastq2cloud/word1-word2-count-ind1 ~/Work/dataset/fastq2cloud/corpus/corpus__ind-counts-word.tsv | 
tr ' ' '\t' > ~/Work/dataset/fastq2cloud/word1-word2-count-ind1-ind2
echo 'Done'

echo 'Build cooccurrence matrix ...'
python ./module_3/cooccurrence_matrix.py  ~/Work/dataset/fastq2cloud/word1-word2-count-ind1-ind2 $dim
echo 'Done'

