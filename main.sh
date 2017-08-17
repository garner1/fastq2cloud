#!/usr/bin/env bash

fastq=$1			# full path to the input fastq.gz file
chunkSize=$2			# how many reads per chunk
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
time bash ./module_2/parse_fastq.sh $fastq $chunkSize $chuncksPrefix
echo "Done"

# These are pickle files TO BE USED IN GENSIM WORD2VEC sequentially, do not concatenate with cat
echo "Prepare corpus ..."	
model=/media/DS2415_seq/silvanog/Genome_segmentation/transitionMatrix_3to3.csv
echo $chuncks
time parallel python $bindir/module_2/create_corpus.py {} $model ::: "$datadir"/chuncks/chunck_*
echo "Done"

echo "Move corpus into specific directory"
mkdir -p "$datadir"/corpus
rm -f "$datadir"/corpus/*
mv "$chuncksPrefix"*_sentences.txt "$datadir"/corpus
echo "Done"

echo 'Parse corpus ...'
corpus="$datadir"/corpus/*
mkdir -p "$datadir"/corpus_summary
parallel "cat {} | LC_ALL=C grep \"'\" | cut -d\"'\" -f2|LC_ALL=C sort |LC_ALL=C uniq -c |awk '{print \$2,\$1}'|tr ' ' '\t' > {.}.counts" ::: "$datadir"/corpus/chunck_*_sentences.txt
cat "$datadir"/corpus/chunck_*_sentences.counts | datamash --sort groupby 1 sum 2 > "$datadir"/corpus_summary/word-counts.tsv
cat $corpus | grep "'" | tr -d '^a' | tr -d '^S' | LC_ALL=C sort -u | 
cat -n | awk '{OFS="\t";print$1,$2}'|LC_ALL=C sort -k2,2 | tr -d "'" > "$datadir"/corpus_summary/corpus__ind-word.tsv
echo 'Done'

echo 'Count pairs ...'
time parallel python "$bindir"/module_3/count_pairs.py {} ::: "$datadir"/corpus/*_sentences.txt 
mv "$datadir"/corpus/*_counter.txt "$datadir"/corpus_summary
echo 'Done'

echo 'Parse cooccurrences ...'
cat "$datadir"/corpus_summary/*_counter.txt | tr -d "(),'" | datamash -t' ' --sort groupby 1,2 sum 3 | tr ' ' '\t' | 
LC_ALL=C sort -k1,1  > "$datadir"/corpus_summary/word1-word2-count

LC_ALL=C join -1 1 -2 2 -o 1.1 1.2 1.3 2.1 "$datadir"/corpus_summary/word1-word2-count "$datadir"/corpus_summary/corpus__ind-word.tsv |
tr ' ' '\t'|LC_ALL=C sort -k2,2 > "$datadir"/corpus_summary/word1-word2-count-ind1

LC_ALL=C join -1 2 -2 2 -o 1.1 1.2 1.3 1.4 2.1 "$datadir"/corpus_summary/word1-word2-count-ind1 "$datadir"/corpus_summary/corpus__ind-word.tsv | 
tr ' ' '\t' > "$datadir"/corpus_summary/word1-word2-count-ind1-ind2
echo 'Done'

# echo 'Build cooccurrence matrix ...'
# dim=`cat "$datadir"/corpus_summary/corpus__ind-word.tsv | wc -l`
# rank=2
# time python "$bindir"/module_3/cooccurrence_AsymMatrix.py  "$datadir"/corpus_summary/word1-word2-count-ind1-ind2 $dim $rank $datadir
# echo 'Done'

####################
# echo 'Run gensim word2vec implementation ...'
# python ./module_3/word2vector.py "$datadir"/corpus
# echo 'Done!
