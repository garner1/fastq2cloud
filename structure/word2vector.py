#!/usr/bin/env python

import sys
import csv
import numpy as np
import pandas as pd
import os.path
import pickle
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import gensim, logging
from os import listdir
from os.path import isfile, join

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

mypath = sys.argv[1]

class MySentences(object):
    def __init__(self, dirname):
        self.dirname = dirname
 
    def __iter__(self):
        for fname in os.listdir(self.dirname):
            for line in open(os.path.join(self.dirname, fname)):
                yield line.split()
 
# sentences = MySentences(mypath) # a memory-friendly iterator
# model = gensim.models.Word2Vec(sentences,min_count=5,workers=32,window=10)

# model.save('/home/garner1/Work/dataset/fastq2cloud/hg19/word2vec_dir/word2vec')

model = gensim.models.Word2Vec.load('/home/garner1/Work/dataset/fastq2cloud/hg19/word2vec_dir/word2vec')
# V = model.wv.vocab.keys()
# V = model.wv.vocab
# print(V)
# print (model.wv['GAGATGGA'])
print (model.wv.most_similar('GAGATGGA'))

# onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
# sentences = []
# for f in listdir(mypath):       # sequentially add all files in directory
#     filename = join(mypath, f)
#     with open(filename,'rb') as f:
#         sentences.extend(pickle.load(f))
#     print filename, len(sentences)
# model = gensim.models.Word2Vec(sentences, min_count=10, workers=32)

# tsne = TSNE(n_components=2)
# X_tsne = tsne.fit_transform(X)

# idx = np.random.randint(X.shape[0], size=10000)
# X_tsne = tsne.fit_transform(X[idx,:])

# plt.scatter(X_tsne[:, 0], X_tsne[:, 1])
# plt.show()

########################################
# model = gensim.models.Word2Vec(iter=1)  # an empty model, no training yet
# model.build_vocab(sentences)  # can be a non-repeatable, 1-pass generator

# with open('/home/garner1/Work/dataset/genome_segmentation/sentence_directory/xab.inputReads_sentences.txt','rb') as f:
#     new_sentences = pickle.load(f)
# model.train(new_sentences,total_examples=model.corpus_count, epochs=model.iter)  # can be a non-repeatable, 1-pass generator

# with open('/home/garner1/Work/dataset/genome_segmentation/sentence_directory/xac.inputReads_sentences.txt','rb') as f:
#     new_sentences = pickle.load(f)
# model.train(new_sentences,total_examples=model.corpus_count, epochs=model.iter)  # can be a non-repeatable, 1-pass generator

