#!/usr/bin/env python

import sys
import csv
import numpy as np
import pandas as pd
from scipy.signal import argrelextrema
from itertools import combinations
from collections import Counter
import os.path
import pickle

chunck = str(sys.argv[1])
model_file = str(sys.argv[2])   # '/media/DS2415_seq/silvanog/Genome_segmentation/transitionMatrix_3to3.csv'

column_names = ['6mer','information']
model = pd.read_csv(model_file, sep=',', names=column_names)
column_names = ['read']
reads = pd.read_csv(chunck, header=None, names=column_names)
sentences = []
for index, row in reads.iterrows():
    kmers = pd.DataFrame([row['read'][i:i+6] for i in range(len(row['read'])-5)]) # the parameter 6 is due to in_dim+out_dim=3+3
    kmers.columns = ['6mer']
    infodata = kmers.merge(model, on='6mer', how='left')
    loc = argrelextrema(infodata.values[:,1], np.greater) # find location of the maxima
    max_df = infodata.iloc[loc]
    new_df = max_df.copy(deep=True)
    index_list = list(max_df.index.values)
    for ind in range(len(index_list)-1):
        deltaLoc = index_list[ind+1] - index_list[ind]
        deltaInfo = float(max_df.loc[index_list[ind+1]]['information']) - float(max_df.loc[index_list[ind]]['information'])
        if (deltaLoc <= 3) and (deltaInfo < 0) and (index_list[ind+1] in new_df.index):
            new_df.drop(index_list[ind+1],inplace=True)
        if (deltaLoc <= 3) and (deltaInfo >= 0) and (index_list[ind] in new_df.index):
            new_df.drop(index_list[ind],inplace=True)
    maxima = list(new_df.index.values)
    if len(maxima)>3:
        read = row['read']
        dictionary = []
        for ind in range(len(maxima)):
            if ind == 0: 
                start = 0
                end = maxima[ind]+3
                dictionary.append(read[start:end])
            if ind > 0 and ind < len(maxima):
                start = maxima[ind-1]+3
                end = maxima[ind]+3
                dictionary.append(read[start:end])
            if ind == len(maxima):
                start = maxima[ind-1]+3
                end = maxima[ind]+3
                dictionary.append(read[start:end])
        start = maxima[-1]+3
        end = len(read)
        dictionary.append(read[maxima[-1]+3:])
        if len(dictionary) >= 2:
            sentences.append(dictionary)

thefilename = open(chunck + "_sentences.txt",'wb')
for item in sentences:
    thefilename.write("%s\n" % item)    
thefilename.close()
            

