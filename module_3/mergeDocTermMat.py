#! /usr/bin/python

import sys
from os import listdir
from os.path import isfile, join
import cPickle as pickle
from scipy.sparse import vstack

mypath = str(sys.argv[1])
filepaths = [join(mypath, f) for f in listdir(mypath) if isfile(join(mypath, f))]

blocks = []
for file_name in filepaths:
    dtm = pickle.load( open( file_name, "rb" ) )
    blocks.append(dtm)

dtm = vstack(blocks)

with open(mypath + "/DTM.pickle", 'wb') as f:
    pickle.dump(dtm, f)
