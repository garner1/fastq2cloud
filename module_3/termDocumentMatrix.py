#! /usr/bin/python

from itertools import chain, combinations
from collections import Counter
import os.path
import sys
import re
import cPickle as pickle
import numpy as np
import scipy.sparse as sp
from sklearn.feature_extraction.text import CountVectorizer

file_name = str(sys.argv[1])    # ex: $datadir/segmented_chr??.dic

docs = []
with open(file_name,"rb") as f:
    for sentence in f:
        docs.append(sentence)

voc = str(sys.argv[2])          # vocabulary txt file
vocabulary = []
with open(voc) as f:
    for word in f:
        vocabulary.append(str(word).rstrip())

vectorizer = CountVectorizer(lowercase=False,vocabulary=vocabulary)
X = vectorizer.fit_transform(docs)

with open(file_name + "_sparseDocTermMat.pickle", 'wb') as f:
    pickle.dump(X, f)

# print(X.nnz)
# print(np.shape(X))
# print(np.shape(vocabulary))


