#!/usr/bin/env bash

fastq=$1			# full path to the input fastq.gz file
chunkSize=$2			# how many lines per chunk (keep it < 50K or you will have memory issues while counting pairs)
bindir=/home/garner1/Work/pipelines/fastq2cloud
exp=`echo $fastq | rev | cut -d'/' -f1 | rev | cut -d'.' -f1`
datadir=/home/garner1/Work/dataset/fastq2cloud/"$exp"

echo 'Running ' $exp ' with output files in ' $datadir

# echo "Create MC model ..." 	# approx 4m to run
# Rscript "$bindir"/module_1/MC_model_from_fasta.R

# # GENERATE FASTQ FILES USING NEAT OR VARSIM_RUN                                                                                                                                              
# # !!! check directories and parameters !!!!                                                                                                                                                 
# # parallel "python ~/tools/neat-genreads/genReads.py -r ~/igv/genomes/hg19.fasta -R 50 -o ~/Work/dataset/genome_segmentation/simulated_data -c 1 --job {} 32" ::: `seq 32`
# # python ~/tools/neat-genreads/mergeJobs.py -i  /home/garner1/tool_test/neat -o  /home/garner1/tool_test -s ~/tools/samtools-1.2

# FASTQ.gz TO READS
echo "Preapare reads ..."
mkdir -p "$datadir"/chuncks
rm -f "$datadir"/chuncks/*
chuncksPrefix="$datadir"/chuncks/chunck_
rm -f "$chuncksPrefix"*
time bash $bindir/module_2/parse_fastq.sh $fastq $chunkSize $chuncksPrefix
echo -e "Done\n"
echo "######"

# These are pickle files TO BE USED IN GENSIM WORD2VEC sequentially, do not concatenate with cat
echo "Prepare corpus ..."	
model=/media/DS2415_seq/silvanog/Genome_segmentation/transitionMatrix_3to3.csv
echo $chuncks
time parallel python $bindir/module_2/create_corpus_2.py {} $model ::: "$datadir"/chuncks/chunck_*
echo -e "Done\n"
echo "######"

echo "Move corpus into specific directory and make the vocabulary"
mkdir -p "$datadir"/corpus
rm -f "$datadir"/corpus/*
mv "$chuncksPrefix"*_sentences.txt "$datadir"/corpus
mkdir -p "$datadir"/corpus_summary
rm -f "$datadir"/corpus_summary/*
corpus="$datadir"/corpus/*
cat $corpus | tr -d "'[]," | tr ' ' '\n' | LC_ALL=C sort -u > "$datadir"/corpus_summary/vocabulary.txt
echo -e 'Done\n'
echo "######"

echo 'Build the Document by Term matrix ...'
time parallel python "$bindir"/module_3/termDocumentMatrix.py {} "$datadir"/corpus_summary/vocabulary.txt ::: "$datadir"/corpus/*_sentences.txt 
mkdir -p "$datadir"/pickle
rm -f "$datadir"/pickle/*
mv "$datadir"/corpus/*_sparseDocTermMat.pickle "$datadir"/pickle
time python $bindir/module_3/mergeDocTermMat.py "$datadir"/pickle
echo -e 'Done\n'
echo "######"
################################
# echo 'Parse cooccurrences ...'
# cat "$datadir"/corpus_summary/*_counter.txt | tr -d "(),'" | datamash -t' ' --sort groupby 1,2 sum 3 | tr ' ' '\t' | 
# LC_ALL=C sort -k1,1  > "$datadir"/corpus_summary/word1-word2-count

# LC_ALL=C join -1 1 -2 2 -o 1.1 1.2 1.3 2.1 "$datadir"/corpus_summary/word1-word2-count "$datadir"/corpus_summary/corpus__ind-word.tsv |
# tr ' ' '\t'|LC_ALL=C sort -k2,2 > "$datadir"/corpus_summary/word1-word2-count-ind1

# LC_ALL=C join -1 2 -2 2 -o 1.1 1.2 1.3 1.4 2.1 "$datadir"/corpus_summary/word1-word2-count-ind1 "$datadir"/corpus_summary/corpus__ind-word.tsv | 
# tr ' ' '\t' > "$datadir"/corpus_summary/word1-word2-count-ind1-ind2

# threshold_repeat=`cat "$datadir"/corpus_summary/word1-word2-count-ind1-ind2 | datamash --sort groupby 1 sum 3 | datamash --sort groupby 2 count 1 | datamash -f max 2 | cut -f1`
# time cat "$datadir"/corpus_summary/word-counts.tsv | awk -v t=$threshold_repeat '$2<t' | cut -f1 |
# parallel --pipe -L1000 --round-robin --compress grep -wFf - "$datadir"/corpus_summary/word1-word2-count-ind1-ind2 | LC_ALL=C sort -u > "$datadir"/corpus_summary/sparse_word1-word2-count-ind1-ind2
# time LC_ALL=C comm -2 -3 <(LC_ALL=C sort "$datadir"/corpus_summary/word1-word2-count-ind1-ind2) <(LC_ALL=C sort $datadir/corpus_summary/sparse_word1-word2-count-ind1-ind2) > $datadir/corpus_summary/dense_word1-word2-count-ind1-ind2

# echo 'Done'
###################################
# Before building the co-occurrence matrix it is important to think about the usefulness of filtering the corpus. 
# By plotting ##w vs #w (i.e.: numb of words which are repeated #w times, for different #w) we see a max value in ##w. 
# The correspoding #v value can be chosen as an automated threshold to separate a sparser regime (left of #w) from a less sparse regime (rigth of #w) in the co-occurrence matrix: M = S+D

# The two matrices S and D can be analyzed separately
######################################

# echo 'Build cooccurrence matrix ...'
# dim=`cat "$datadir"/corpus_summary/corpus__ind-word.tsv | wc -l`
# rank=2
# time python "$bindir"/module_3/cooccurrence_AsymMatrix.py  "$datadir"/corpus_summary/word1-word2-count-ind1-ind2 $dim $rank $datadir
# echo 'Done'

####################
# echo 'Run gensim word2vec implementation ...'
# python ./module_3/word2vector.py "$datadir"/corpus
# echo 'Done!
###########################################
# cat "$datadir"/corpus_summary/word-counts.tsv | awk -v t=$threshold_repeat '$2<t' | cut -f1 | split -l 1000 --additional-suffix=.input - "$datadir"/corpus_summary/sparse_vocabulary_
# cat $datadir/corpus_summary/word1-word2-count-ind1-ind2 | parallel 'grep -f {} - > {.}.output' ::: "$datadir"/corpus_summary/sparse_vocabulary_*.input
# cat "$datadir"/corpus_summary/sparse_vocabulary_*.output | LC_ALL=C sort -u > "$datadir"/corpus_summary/sparse_word1-word2-count-ind1-ind2
# rm -f "$datadir"/corpus_summary/sparse_vocabulary_*
