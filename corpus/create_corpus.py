#!/usr/bin/env python
#smoothed z-score thresholding algorithm

import sys
import csv
import numpy as np
import pandas as pd
from scipy.signal import argrelextrema
from itertools import combinations
from collections import Counter
import os.path
import pickle
import pylab
import itertools
import random

def thresholding_algo(y, lag, threshold, influence):
    signals = np.zeros(len(y))
    filteredY = np.array(y)
    avgFilter = [0]*len(y)
    stdFilter = [0]*len(y)
    avgFilter[lag - 1] = np.mean(y[0:2*lag]) # set the mean to that of the window of 2lag size
    stdFilter[lag - 1] = np.std(y[0:2*lag])  # set the std to that of the window of 2lag size
    for i in range(lag, len(y)-lag):
        if (y[i] - avgFilter[i - 1]) > threshold * stdFilter [i-1] and y[i-1] < y[i]: # condition to call a stop word
            signals[i] = 1      # call signal
            filteredY[i] = influence * y[i] + (1 - influence) * filteredY[i-1] # update local data points
            avgFilter[i] = np.mean(filteredY[(i-lag):(i+lag)]) # update the local mean
            stdFilter[i] = np.std(filteredY[(i-lag):(i+lag)])  # update the local std
        else:                                                  # if not stop word
            signals[i] = 0                                     # is a word
            filteredY[i] = y[i]                                # update local data
            avgFilter[i] = np.mean(filteredY[(i-lag):(i+lag)]) # update local mean
            stdFilter[i] = np.std(filteredY[(i-lag):(i+lag)])  # update local std

    return dict(signals = np.asarray(signals), # return list of words or stop words
                avgFilter = np.asarray(avgFilter), # list of local averages
                stdFilter = np.asarray(stdFilter)) # list of local stds

#######################################################################################################

chunck = str(sys.argv[1])       # list of reads
model_file = str(sys.argv[2])   # file with MC parameters
dim = int(sys.argv[3])          # MC dimension
keep_N = int(sys.argv[4])      # if True (1) randomly replace letter N; if False (0) discard read

column_names = ['kmer','information']
model = pd.read_csv(model_file, sep='\t', names=column_names) # MC dataframe

column_names = ['read']
reads = pd.read_csv(chunck, header=None, names=column_names) # reads dataframe

sentences = []
alphabet = ['C','G','A','T']

# Settings parameters
lag = dim+3                     # window size where to consider local average and std
threshold = 1                   # above this consider it a stop word
influence = 0.5                 # influence of the stop-word on the mooving average
sentences = []
for index, row in reads.iterrows():
    print row['read']
    if keep_N :                 # BUG IN THIS LOOP IN CASE OF FALSE!!!!
        read = row['read'].replace('N',random.choice(alphabet)) # random substitution of letter N (OR IS IT BETTER TO DISCARD THE READS?)
    if 'N' in row['read'] and not keep_N :
        continue                # discard this read
    else:
        read = row['read']
    kmers = pd.DataFrame([read[i:i+dim+3] for i in range(len(read)-(dim+3-1))]) # create list of kmers from read
    kmers.columns = ['kmer']
    infodata = kmers.merge(model, on='kmer', how='left').fillna(0) # JOIN KMERS DF AND ENTROPY VALUES DF(model)
    y = infodata.values[:,1].astype(np.float)
    y = np.lib.pad(y, (dim,dim-1), 'constant', constant_values=(0)) # set to 0 the borders of the reads

    result = thresholding_algo(y, lag=lag, threshold=threshold, influence=influence)
    stopword = [i for i, e in enumerate(result["signals"]) if e != 0] # find the positions of the stop words
    for pos in stopword:        # place a comma in the read at the position of the stop word
        read = read[:pos] + ',' + read[pos + 1:]
    sentence = [word for word in read.split(',')[1:-1] if len(word) > 3] # build sentence discarding words at the border and short words
    sentences.append(sentence)

    # if index == 1:
    #     # Plot result
    #     pylab.plot(np.arange(1, len(y)+1), y)
    #     pylab.plot(np.arange(1, len(y)+1), result["avgFilter"], color="cyan", lw=2)
    #     pylab.plot(np.arange(1, len(y)+1), result["avgFilter"] + threshold * result["stdFilter"], color="green", lw=2)
    #     pylab.plot(np.arange(1, len(y)+1), result["avgFilter"] - threshold * result["stdFilter"], color="green", lw=2)
    #     pylab.scatter(np.arange(1, len(y)+1), result["signals"]*y, color="red", lw=1)
    #     pylab.show()

thefilename = open(chunck + "_sentences.txt",'wb')
for item in sentences:
    thefilename.write("%s\n" % item)    
thefilename.close()
