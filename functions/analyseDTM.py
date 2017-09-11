#! /usr/bin/python

import sys
from os import listdir
from os.path import isfile, join
import cPickle as pickle
from scipy.sparse import vstack
import numpy as np
import matplotlib.pyplot as plt

mypath = str(sys.argv[1])
dtm = pickle.load( open(mypath + "/DTM.pickle", 'rb') )

x = np.squeeze(np.asarray(dtm.sum(0)))

# plt.hist(x, normed=True, bins=500)
# plt.ylabel('Probability')
# plt.xlim(0,500)
# plt.show()

import collections
counter=collections.Counter(x)
# print(counter)
# # Counter({1: 4, 2: 4, 3: 2, 5: 2, 4: 1})
# print(counter.values())
# # [4, 4, 2, 1, 2]
# print(counter.keys())
# # [1, 2, 3, 4, 5]
# print(counter.most_common(3))
# # [(1, 4), (2, 4), (3, 2)]
with open(mypath + "/countDoc-freq.txt",'wb') as f:
    for k,v in  counter.most_common():
        f.write( "{} {}\n".format(k,v) )
