#!/usr/bin/env bash

# simulation_flag=$1		# to simulate (0) sequencing or not (1); remember that convention is different from other languages

# if [[ $simulation_flag == 0 ]]; then # simulate sequencing
#     rnd_seed=$2			     # random seed for the sequencing simulation
#     genome_dir=$3		     # directory location of the reference genome
#     genome_name=$4		     # name of the genome
#     coverage=$5			     # desired coverage
#     keep_N=$6                        # 0 if you want to discard reads with N; 1 if you want to randomly replace them
#     exp="simseq_"$coverage"X_"$genome_name"_rs"$rnd_seed # name of the experiment
# fi

# if [[ $simulation_flag == 1 ]]; then # do not simulate sequencing
#     # fastq=$2                         # full path to the input fastq.gz file
#     # keep_N=$3                        # 0 if you want to discard reads with N; 1 if you want to randomly replace them
#     # exp=`echo $fastq | rev | cut -d'/' -f1 | rev | cut -d'.' -f1` # name of the experiment
# fi

exp=hg19
bindir=/home/garner1/Work/pipelines/fastq2cloud        # location of the executable files
datadir=/home/garner1/Work/dataset/fastq2cloud/"$exp"/ # location of the output files
# mkdir -p "$datadir"fastq			       

echo 'Running ' $exp ' with output files in ' $datadir

# if [[ $simulation_flag == 1 ]]; then # do not simulate sequencing
#     # gunzip -c $fastq > "$datadir"fastq/unzip.fastq
#     # echo "Quality control on the fastq file ..."
#     # sickle se -f "$datadir"fastq/unzip.fastq -t sanger -o "$datadir"fastq/trimmed.fastq -n
#     # echo "Done"
#     # fastq="$datadir"fastq/trimmed.fastq
# fi

# if [[ $simulation_flag == 0 ]]; then # simulate sequencing
#     # echo "Simulate sequencing ..."
#     # parallel ~/tools/ART/art_bin_MountRainier/art_illumina -rs $rnd_seed -ss HS25 -i {} -o {.}.1X_SE -l 150 -f $coverage ::: "$genome_dir"/*.fa #with -f coverage
#     # mkdir -p "$datadir"fastq && rm -f "$datadir"fastq/*
#     # mv "$genome_dir"/*.{fq,aln} "$datadir"fastq && rm -f "$datadir"fastq/merged.fq
#     # cat "$datadir"fastq/*.fq > "$datadir"fastq/merged.fastq && rm -f "$datadir"fastq/*.fq
#     # parallel gzip {} ::: "$datadir"fastq/*.aln 
#     # fastq="$datadir"/fastq/merged.fastq
#     # echo "Done"
# fi

# echo "Preapare reads ..."
# # readsOnly="$datadir"chuncks/
# # mkdir -p "$datadir"/chuncks && rm -f "$datadir"/chuncks/*
# # bash $bindir/corpus/parse_fastq.sh $fastq $readsOnly
# echo "Done"

# echo "Split the fastq file ..."
# # rm -rf "$datadir"splitFastq && mkdir -p "$datadir"splitFastq
# # cat $fastq | split -l 4000000 - "$datadir"splitFastq/ # split to process 1M reads per time with jellyfish
# echo "Done"

#!!! an option is to use kmc3
#!!!CONSIDER BUILDING THE MODEL FROM A REFERENCE GENOME
# echo "Create MC model ..." 	# it takes 30min on BB23

######################################################
# # USE THIS IF MODEL IS CREATED FROM FASTQ FILE:
# rm -f "$datadir"??.jf
# for dim in `seq 6 2 10`; do
#     echo "count"
#     for file in `ls "$datadir"splitFastq/`; do
# 	input="$datadir"splitFastq/"$file"
# 	bash "$bindir"/segmentation/jelly_f1.sh $dim $input "$datadir""$file".jf
#     done 
#     echo "merge"
#     jellyfish merge -o "$datadir"merged_"$dim".jf "$datadir"??.jf
#     rm -f "$datadir"??.jf
# done
# time parallel bash "$bindir"/segmentation/jelly_f2.sh "$datadir"merged_{}.jf {} 1.0 ::: 6 8 10
######################################################

##########################################################
# USE THIS IF MODEL IS CREATED FROM REFERENCE GENOME:
# jellyfish count -m 9 -s 3G -t 32 --disk -o "$datadir"9mer_hg19.jf ~/igv/genomes/hg19.fasta
# jellyfish count -m 10 -s 3G -t 32 --disk -o "$datadir"10mer_hg19.jf ~/igv/genomes/hg19.fasta
# jellyfish count -m 11 -s 3G -t 32 --disk -o "$datadir"11mer_hg19.jf ~/igv/genomes/hg19.fasta
# jellyfish count -m 12 -s 3G -t 32 --disk -o "$datadir"12mer_hg19.jf ~/igv/genomes/hg19.fasta
# parallel bash "$bindir"/segmentation/jelly_f2.sh "$datadir"{}mer_hg19.jf {} 1.0 ::: 9 10 11 12

# rm -rf "$datadir"MCmodel		# clean directory
# mkdir -p "$datadir"MCmodel 
# mv "$datadir"*.csv "$datadir"MCmodel
# rm -rf "$datadir"*.{fa,tsv,jf}
# echo "Done"
##########################################################

# echo "Prepare corpus ..."
# g++ -std=c++11 $bindir/corpus/tokenizer_withMean.cpp -o $bindir/corpus/mean & pid1=$! # compile the tokenizer
# g++ -std=c++11 $bindir/corpus/tokenizer_withMean_RC.cpp -o $bindir/corpus/meanRC & pid2=$! # compile the tokenizer
# wait $pid1
# wait $pid2
# time parallel "./corpus/mean {} '$datadir'MCmodel/ > {}_sentences" ::: "$datadir"chuncks/??
# time parallel "./corpus/meanRC {} '$datadir'MCmodel/ > {}_RCsentences" ::: "$datadir"chuncks/??
# echo "Done"

# echo "Move corpus into specific directory and make the vocabulary"
# mkdir -p "$datadir"corpus && rm -f "$datadir"corpus/*
# mv "$readsOnly"*_sentences "$datadir"corpus
# mkdir -p "$datadir"corpus_summary && rm -f "$datadir"corpus_summary/*

# corpus="$datadir"corpus/*
# cat $corpus | 
# # cut -d',' -f2-| 	# remove first word
# tr -d "'[]" | tr ',' '\n' | LC_ALL=C sort | LC_ALL=C uniq -c | 
# awk '{print $1"\t"$2}' > "$datadir"corpus_summary/count_word.txt 

# bash ./functions/word_frequency.sh "$datadir"corpus_summary/count_word.txt

# awk '{print $2}' "$datadir"corpus_summary/count_word.txt > "$datadir"corpus_summary/vocabulary.txt

# awk '{print length($1)}' "$datadir"corpus_summary/vocabulary.txt | LC_ALL=C sort -n | LC_ALL=C uniq -c | 
# gnuplot -p -e 'set terminal pdf; set output "wordlen-freq.pdf";set logscale y; plot "/dev/stdin" using 2:1'
# mv wordlen-freq.pdf "$datadir"corpus_summary
# echo 'Done'

# echo 'Build the Document by Term matrix ...'
# word_len_threshold=1000		# max word length
# # python "$bindir"/structure/termDocumentMatrix.py "$datadir"corpus/Y_sentences "$datadir"corpus_summary/vocabulary.txt $word_len_threshold
# time parallel python "$bindir"/structure/termDocumentMatrix.py {} "$datadir"corpus_summary/vocabulary.txt $word_len_threshold ::: "$datadir"corpus/??_sentences

# mkdir -p "$datadir"pickle && rm -f "$datadir"pickle/*
# mv "$datadir"corpus/*_sparseDocTermMat.pickle "$datadir"pickle

# time python $bindir/structure/mergeDocTermMat.py "$datadir"pickle

# rm -f "$datadir"pickle/*_sentences_sparseDocTermMat.pickle
# echo 'Done'

# echo 'Build the co-occurrence matrix ...'
# time python $bindir/structure/cooccurrenceMat.py "$datadir"pickle/DTM.pickle
# echo 'Done'

####################
# echo 'Run gensim word2vec implementation ...'
# python ./structure/word2vector.py "$datadir"corpus
# echo 'Done!
###########################################
