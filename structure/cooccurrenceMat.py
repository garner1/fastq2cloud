#! /usr/bin/python

import sys
from os import listdir
from os.path import isfile, join
import cPickle as pickle
from scipy.sparse import vstack
import numpy as np

file_name = str(sys.argv[1])

dtm = pickle.load( open( file_name, "rb" ) )

cooc = dtm.transpose()*dtm

with open(file_name + "_cooccurrenceMat.pickle", 'wb') as f:
    pickle.dump(cooc, f)
