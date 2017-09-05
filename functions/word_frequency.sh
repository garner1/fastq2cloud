#!/usr/bin/env bash

file=$1				# text file with count-word structure


path=`echo $file | rev | cut -d'/' -f2- | rev`

python word_frequency.py $file > "$path"/power-law_parameters.txt
alpha_dagger=`cat "$path"/power-law_parameters.txt | grep C | cut -d' ' -f3`
tail -n+2 "$path"/count_word.txt | cut -f1 | LC_ALL=C sort -nr | cat -n | datamash -s -g 2 max 1 > "$path"/freq-rank.dat
cat "$path"/freq-rank.dat | 
gnuplot -p -e 'set logscale x;set logscale y;set xlabel "frequency"; set ylabel "numb of words occurring at least x";plot "/dev/stdin" using 1:2'
cat "$path"/freq-rank.dat | 
gnuplot -e 'set terminal pdf; set output "cumulative.pdf";set logscale x;set logscale y;set xlabel "frequency"; set ylabel "numb of words occurring at least x"; plot "/dev/stdin" using 1:2' 
mv cumulative.pdf $path

xpdf $path/freq-prob.pdf 
