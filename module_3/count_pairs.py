#! /usr/bin/python

from itertools import chain, combinations
from collections import Counter
import os.path
import sys
import re
import cPickle as pickle
import numpy as np
import scipy.sparse as sp

file_name = str(sys.argv[1])    # ex: $datadir/segmented_chr??.dic

corpus = pickle.load( open( file_name, "rb" ) )

###########################
# Note that this will only work if the elements are in the same order in every list, 
# e.g. Counter(chain.from_iterable(combinations(x,2) for x in [[1,2],[2,1]])) produces Counter({(1, 2): 1, (2, 1): 1}). 
# If you want to count each pairing regardless of order, turn each list into a set first: 
# Counter(chain.from_iterable(combinations(x,2) for x in [set([1,2]),set([2,1])])) produces Counter({(1, 2): 2}):
##########################
# counts = Counter(chain.from_iterable(combinations(sentence, 2) for sentence in corpus))
counts = Counter(chain.from_iterable(combinations(set(sentence), 2) for sentence in corpus))

# mat = sp.dok_matrix((9807,9807), dtype=np.int8)
# for pair, value in counts.items():
#     mat[pair[0],pair[1]] = value
# print mat.shape
##################################

with open(file_name + "_counter.txt",'w') as f:
    for k,v in  counts.most_common():
        f.write( "{} {}\n".format(k,v) )

