#!/usr/bin/env bash
#TO BE CALLED FROM MAIN DIRECTORY
file=$1				# text file with count-word structure

name=`echo $file | rev | cut -d'.' -f2- | rev`
path=`echo $file | rev | cut -d'/' -f2- | rev`

python $PWD/functions/word_frequency.py $file > "$name"_power-law_parameters.txt

alpha_dagger=`cat "$name"_power-law_parameters.txt | grep C | cut -d' ' -f3`

tail -n+2 "$path"/count_word.txt | cut -f1 | LC_ALL=C sort -nr | cat -n | datamash -s -g 2 max 1 > "$name"_freq-rank.dat
# cat "$name"_freq-rank.dat | gnuplot -p -e 'set logscale x;set logscale y;set xlabel "frequency"; set ylabel "numb of words occurring at least x";plot "/dev/stdin" using 1:2'
cat "$name"_freq-rank.dat | 
gnuplot -e 'set terminal pdf; set output "cumulative.pdf";set logscale x;set logscale y;set xlabel "frequency"; set ylabel "numb of words occurring at least x"; plot "/dev/stdin" using 1:2' 
mv cumulative.pdf "$name"_cumulative.pdf

# xpdf $path/freq-prob.pdf 
