#!/usr/bin/env python
#smoothed z-score thresholding algorithm

import sys
import csv
import numpy as np
import pandas as pd
from scipy.signal import argrelextrema
from itertools import combinations
from collections import Counter
from os.path import isfile, join
import pickle
import pylab
import itertools
import random
from os import listdir

def thresholding_algo(y, lag, threshold, influence,dim):
    signals = np.zeros(len(y))  # initialize the stop-words vector
    filteredY = np.array(y)     # initialize the filtered vector to the entire data
    avgFilter = [0]*len(y)      # init the vector of mean values
    stdFilter = [0]*len(y)      # init the vector of std values
    avgFilter[lag - 1] = np.mean(y[0:dim]) # set the first mean to that of the window of dim size
    stdFilter[lag - 1] = np.std(y[0:dim])  # set the first std to that of the window of dim size
    for i in range(lag, len(y) - (dim-1)):
        if (y[i] - avgFilter[i - 1]) > threshold * stdFilter [i-1] and y[i-1] < y[i]: # condition to call a stop word
            signals[i] = 1      # call signal
            filteredY[i] = influence * y[i] + (1 - influence) * filteredY[i-1] # update local data points
            avgFilter[i] = np.mean(filteredY[(i-dim):(i+dim)]) # update the local mean
            stdFilter[i] = np.std(filteredY[(i-dim):(i+dim)])  # update the local std
        else:                                                  # if not stop word
            signals[i] = 0                                     # is a word
            filteredY[i] = y[i]                                # update local data
            avgFilter[i] = np.mean(filteredY[(i-dim):(i+dim)]) # update local mean
            stdFilter[i] = np.std(filteredY[(i-dim):(i+dim)])  # update local std

    return dict(signals = np.asarray(signals), # return list of words or stop words
                avgFilter = np.asarray(avgFilter), # list of local averages
                stdFilter = np.asarray(stdFilter)) # list of local stds

#######################################################################################################
chunck = str(sys.argv[1])       # list of reads
mypath = str(sys.argv[2])
keep_N = int(sys.argv[3])      # if True (1) randomly replace letter N; if False (0) discard read with N

files = [join(mypath, f) for f in listdir(mypath) if isfile(join(mypath, f))] # files with MC parameters
dims = [files[ind].split('/')[-1].split('t')[0].split('_')[-1] for ind in range(len(files))] # MC kmer size
column_names = ['kmer','information']
model = [pd.read_csv(model_file, sep='\t', names=column_names) for model_file in files] # MC dataframes

column_names = ['read']
reads = pd.read_csv(chunck, header=None, names=column_names) # reads dataframe

sentences = []
alphabet = ['C','G','A','T']

# Settings parameters
threshold = 1                   # above this consider it a stop word
influence = 0.5                 # influence of the stop-word on the mooving average
sentences = []
for index, row in reads.iterrows():
    if keep_N :                 # BUG IN THIS LOOP IN CASE OF FALSE!!!!
        read = row['read'].replace('N',random.choice(alphabet)) # random substitution of letter N (OR IS IT BETTER TO DISCARD THE READS?)
    if 'N' in row['read'] and not keep_N :
        continue                # discard this read
    else:
        read = row['read']
    ind = 0
    dic_of_mean = {}
    dic_of_values = {}
    for dim in dims:
        dim = int(dim)
        kmers = pd.DataFrame([read[i:i+dim] for i in range(len(read) - (dim + 1))]) # create list of kmers from read
        kmers.columns = ['kmer']
        infodata = kmers.merge(model[ind], on='kmer', how='left').fillna(0) # JOIN KMERS DF AND ENTROPY VALUES DF(model)
        ind += 1
        dic_of_values[dim] = infodata.values[:,1].astype(np.float)
        aux = [i for i in dic_of_values[dim] if i > 0] # evaluate the mean only of the positive values, since zeros are due to low statistics
        if len(aux)>0:
            dic_of_mean[dim] = np.mean(aux)
        if len(dic_of_mean) > 0:
            dim = min(dic_of_mean, key=dic_of_mean.get) # choose the model with the minimum average entropy, since this is the best model describing the data from an info theoretic POV
        else:
            continue
    y = dic_of_values[dim]                      # set the values associated to the minimum average entropy

    lag = dim+10                         # lag before initializing local average and std
    result = thresholding_algo(y, lag=lag, threshold=threshold, influence=influence, dim=dim)
    stopword = [i for i, e in enumerate(result["signals"]) if e != 0] # find the positions of the stop words
    for pos in stopword:        # place a comma in the read at the position of the stop word
        read = read[:pos] + ',' + read[pos + 1:]
    sentence = [word for word in read.split(',')[1:-1] if len(word) > 3] # build sentence discarding words at the border and short words
    sentences.append(sentence)

    # if index < 3:
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
