#!/usr/bin/env bash

fastq=$1		# full path to the input fastq.gz file
linesperfile=$2		# numb of lines per chunk: min{50000,#reads/20}
dim=$3		        # MC in-dimension parameter
keep_N=$4		# 0 if you want to discard reads with N; 1 if you want to randomly replace them

bindir=/home/garner1/Work/pipelines/fastq2cloud
exp=`echo $fastq | rev | cut -d'/' -f1 | rev | cut -d'.' -f1`
datadir=/home/garner1/Work/dataset/fastq2cloud/"$exp"/

echo 'Running ' $exp ' with output files in ' $datadir

echo "Create MC model ..." 	
kmer=`echo $(( 3 + $dim ))`
time bash "$bindir"/segmentation/jelly.sh $fastq $kmer $dim $datadir
echo "Done"

echo "Preapare reads ..."
mkdir -p "$datadir"/chuncks && rm -f "$datadir"/chuncks/*
chuncksPrefix="$datadir"/chuncks/chunck_ && rm -f "$chuncksPrefix"*
time bash $bindir/corpus/parse_fastq_1.sh $fastq $linesperfile $chuncksPrefix
echo "Done"

echo "Prepare corpus ..."	
model=$datadir/transitionMatrix_fromFastq_"$dim"to3.csv
time parallel python $bindir/corpus/create_corpus.py {} $model $dim $keep_N ::: "$datadir"/chuncks/chunck_*
echo "Done"

echo "Move corpus into specific directory and make the vocabulary"
mkdir -p "$datadir"/corpus
rm -f "$datadir"/corpus/*
mv "$chuncksPrefix"*_sentences.txt "$datadir"/corpus
mkdir -p "$datadir"/corpus_summary
rm -f "$datadir"/corpus_summary/*
corpus="$datadir"/corpus/*

time cat $corpus | tr -d "'[]," | tr ' ' '\n' | LC_ALL=C sort | LC_ALL=C uniq -c | awk '{print $1"\t"$2}' > "$datadir"/corpus_summary/count_word.txt
bash ./functions/word_frequency.sh "$datadir"/corpus_summary/count_word.txt

awk '{print $2}' "$datadir"/corpus_summary/count_word.txt > "$datadir"/corpus_summary/vocabulary.txt
awk '{print $1}' "$datadir"/corpus_summary/count_word.txt | LC_ALL=C sort -nr|cat -n|awk '{print $1"\t"$2}'|
datamash -s groupby 2 max 1 > "$datadir"/corpus_summary/rank_frequency.dat

awk '{print length($1)}' "$datadir"/corpus_summary/vocabulary.txt | LC_ALL=C sort -n | LC_ALL=C uniq -c | 
gnuplot -p -e 'set terminal pdf; set output "wordlen-freq.pdf";set logscale y; plot "/dev/stdin" using 2:1'
mv wordlen-freq.pdf "$datadir"/corpus_summary
echo 'Done'

# echo 'Build the Document by Term matrix ...'
# time parallel python "$bindir"/structure/termDocumentMatrix.py {} "$datadir"/corpus_summary/vocabulary.txt $word_len_threshold ::: "$datadir"/corpus/*_sentences.txt 
# mkdir -p "$datadir"/pickle
# rm -f "$datadir"/pickle/*
# mv "$datadir"/corpus/*_sparseDocTermMat.pickle "$datadir"/pickle
# time python $bindir/structure/mergeDocTermMat.py "$datadir"/pickle
# rm -f "$datadir"/pickle/chunck*pickle
# echo 'Done'

# echo 'Build the co-occurrence matrix ...'
# # The total number of docs can be < than the initial numb of reads because some reads might contain too few words
# time python $bindir/structure/cooccurrenceMat.py "$datadir"/pickle/DTM.pickle
# echo 'Done'

#######################################################
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
# time python "$bindir"/structure/cooccurrence_AsymMatrix.py  "$datadir"/corpus_summary/word1-word2-count-ind1-ind2 $dim $rank $datadir
# echo 'Done'

####################
# echo 'Run gensim word2vec implementation ...'
# python ./structure/word2vector.py "$datadir"/corpus
# echo 'Done!
###########################################
# cat "$datadir"/corpus_summary/word-counts.tsv | awk -v t=$threshold_repeat '$2<t' | cut -f1 | split -l 1000 --additional-suffix=.input - "$datadir"/corpus_summary/sparse_vocabulary_
# cat $datadir/corpus_summary/word1-word2-count-ind1-ind2 | parallel 'grep -f {} - > {.}.output' ::: "$datadir"/corpus_summary/sparse_vocabulary_*.input
# cat "$datadir"/corpus_summary/sparse_vocabulary_*.output | LC_ALL=C sort -u > "$datadir"/corpus_summary/sparse_word1-word2-count-ind1-ind2
# rm -f "$datadir"/corpus_summary/sparse_vocabulary_*
