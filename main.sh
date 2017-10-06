#!/usr/bin/env bash

simulation_flag=$1		# to simulate (0) sequencing or not (1); remember that convention is different from other languages

if [[ $simulation_flag == 0 ]]; then
    rnd_seed=$2			# random seed for the sequencing simulation
    genome_dir=$3		# directory location of the reference genome
    genome_name=$4		# name of the genome
    coverage=$5			# desired coverage
    keep_N=$6                   # 0 if you want to discard reads with N; 1 if you want to randomly replace them
    exp="simseq_"$coverage"X_"$genome_name"_rs"$rnd_seed # name of the experiment
fi

if [[ $simulation_flag == 1 ]]; then
    fastq=$2           # full path to the input fastq.gz file
    linesperfile=$3    # numb of lines per chunk: min{50000,#reads/20}
    keep_N=$4          # 0 if you want to discard reads with N; 1 if you want to randomly replace them
    exp=`echo $fastq | rev | cut -d'/' -f1 | rev | cut -d'.' -f1` # name of the experiment
fi

bindir=/home/garner1/Work/pipelines/fastq2cloud        # location of the executable files
datadir=/home/garner1/Work/dataset/fastq2cloud/"$exp"/ # location of the output files

echo 'Running ' $exp ' with output files in ' $datadir

if [[ $simulation_flag == 0 ]]; then
    echo "Simulate sequencing ..."
    parallel ~/tools/ART/art_bin_MountRainier/art_illumina -rs $rnd_seed -ss HS25 -i {} -o {.}.1X_SE -l 150 -f $coverage ::: "$genome_dir"/*.fa #with -f coverage
    mkdir -p $datadir/fastq && rm -f $datadir/fastq/*
    mv "$genome_dir"/*.{fq,aln} $datadir/fastq && rm -f $datadir/fastq/merged.fq
    cat $datadir/fastq/*.fq > $datadir/fastq/merged.fastq && rm -f $datadir/fastq/*.fq
    parallel gzip {} ::: $datadir/fastq/*.aln $datadir/fastq/merged.fastq
    linesperfile=50000
    fastq="$datadir"/fastq/merged.fastq.gz
    echo "Done"
fi

echo "Preapare reads ..."
mkdir -p "$datadir"/chuncks && rm -f "$datadir"/chuncks/*
chuncksPrefix="$datadir"/chuncks/chunck_ 
time bash $bindir/corpus/parse_fastq_1.sh $fastq $linesperfile $chuncksPrefix
echo "Done"

mkdir -p "$datadir"splitFastq
zcat $fastq | split -l 4000000 - "$datadir"splitFastq/
echo "Create MC model ..." 
for dim in `seq 10 2 14`; do
    for file in `ls "$datadir"splitFastq/`; do
	input="$datadir"splitFastq/"$file"
	time bash "$bindir"/segmentation/jelly_f1.sh $dim $input "$datadir""$file".jf
    done 
    jellyfish merge -o "$datadir"merged_"$dim".jf "$datadir"a*.jf
    rm -f "$datadir""$file".jf
done
time parallel bash "$bindir"/segmentation/jelly_f2.sh "$datadir"/merged_{}.jf {} 1.0 ::: 10 12 14
rm -rf $datadir/MCmodel		# rm old directory
mkdir -p $datadir/MCmodel && mv "$datadir"*.csv $datadir/MCmodel
rm -f $datadir/{_*,merged_*}	# clean up
modelDirectory=$datadir/MCmodel # create the directory storing the MC models
echo "Done"
# !!!INTRODUCE THE RIGHT TO LEFT AND LEFT TO RIGTH REPRESENTATION!!!

# echo "Prepare corpus ..."
# g++ -std=c++11 $bindir/corpus/tokenizer_withMean.cpp -o $bindir/corpus/mean & pid1=$! # compile the tokenizer
# g++ -std=c++11 $bindir/corpus/tokenizer_withMedian.cpp -o $bindir/corpus/median & pid2=$! # compile the tokenizer
# wait $pid1
# wait $pid2
# # time parallel "./corpus/mean {} '$datadir'MCmodel/ > {}_MeanSentences" ::: "$datadir"chuncks/chunck_*
# time parallel "./corpus/median {} '$datadir'MCmodel/ > {}_MeanSentences" ::: "$datadir"chuncks/chunck_*
# echo "Done"

################################################################
#!!! NEED TO REMOVE FIRST WORD AND SHORT WORDS IN DOCUMENTS  !!!
################################################################

# echo "Move corpus into specific directory and make the vocabulary"
# mkdir -p "$datadir"/corpus && rm -f "$datadir"/corpus/*
# mv "$chuncksPrefix"*_sentences.txt "$datadir"/corpus
# mkdir -p "$datadir"/corpus_summary && rm -f "$datadir"/corpus_summary/*
# corpus="$datadir"/corpus/*

# time cat $corpus | tr -d "'[]," | tr ' ' '\n' | LC_ALL=C sort | LC_ALL=C uniq -c | awk '{print $1"\t"$2}' > "$datadir"/corpus_summary/count_word.txt
# bash ./functions/word_frequency.sh "$datadir"/corpus_summary/count_word.txt

# awk '{print $2}' "$datadir"/corpus_summary/count_word.txt > "$datadir"/corpus_summary/vocabulary.txt
# awk '{print $1}' "$datadir"/corpus_summary/count_word.txt | LC_ALL=C sort -nr|cat -n|awk '{print $1"\t"$2}'|
# datamash -s groupby 2 max 1 > "$datadir"/corpus_summary/rank_frequency.dat

# awk '{print length($1)}' "$datadir"/corpus_summary/vocabulary.txt | LC_ALL=C sort -n | LC_ALL=C uniq -c | 
# gnuplot -p -e 'set terminal pdf; set output "wordlen-freq.pdf";set logscale y; plot "/dev/stdin" using 2:1'
# mv wordlen-freq.pdf "$datadir"/corpus_summary
# echo 'Done'

# echo 'Build the Document by Term matrix ...'
# word_len_threshold=100		# max word length
# time parallel python "$bindir"/structure/termDocumentMatrix.py {} "$datadir"/corpus_summary/vocabulary.txt $word_len_threshold ::: "$datadir"/corpus/*_sentences.txt 
# mkdir -p "$datadir"/pickle && rm -f "$datadir"/pickle/*
# mv "$datadir"/corpus/*_sparseDocTermMat.pickle "$datadir"/pickle
# time python $bindir/structure/mergeDocTermMat.py "$datadir"/pickle
# rm -f "$datadir"/pickle/chunck*pickle
# echo 'Done'

# echo 'Build the co-occurrence matrix ...'
# # The total number of docs can be < than the initial numb of reads because some reads might contain too few words
# time python $bindir/structure/cooccurrenceMat.py "$datadir"/pickle/DTM.pickle
# echo 'Done'

####################
# echo 'Run gensim word2vec implementation ...'
# python ./structure/word2vector.py "$datadir"/corpus
# echo 'Done!
###########################################
