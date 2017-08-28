#! /usr/bin/python

from scipy.sparse.linalg import svds
from scipy.sparse import *
import os.path
import sys
from numpy import *
 
file_name = str(sys.argv[1])    # "$datadir"/corpus_summary/word1-word2-count-ind1-ind2
dim = int(sys.argv[2])
rank = int(sys.argv[3])
datadir = str(sys.argv[4])


coomat = dok_matrix((dim, dim), dtype=int)
with open(file_name, 'rb') as infile:
    for line in infile:
        lista = line.strip().split('\t')
        coomat[int(lista[3])-1,int(lista[4])-1] = int(lista[2])

print('Matrix built')
D = coomat.nnz                                               # number of pairs
rowsum = coomat.sum(axis=0)                                  # sum along rows
colsum = coomat.sum(axis=1)                                  # sum along columns

del(coomat)
print('Matrix deleted')
coomat = dok_matrix((dim, dim))
with open(file_name, 'rb') as infile:
    for line in infile:
        lista = line.strip().split('\t')
        data = log(D*float(lista[2])*1.0/float(rowsum[0,int(lista[4])-1]*colsum[int(lista[3])-1,0]))
        if data > 0:
            coomat[int(lista[3])-1,int(lista[4])-1] = data
print('Matrix rebuilt')
# coomat = coomat + transpose(coomat)
# Build the PMI:
# coomat.data = ma.log(coomat.data) #use ma library to mask invalid values, like log(0)
# Build the PPMI:
# coomat.data = ma.masked_less(coomat.data, 0) # mask values less than 0

[u,s,vt] = svds(coomat, rank, which='LM', return_singular_vectors=True)

u.dump(datadir+"/corpus_summary/U.dat")
# u = load(datadir+"/corpus_summary/u.dat")
s.dump(datadir+"/corpus_summary/S.dat")
# s = load(datadir+"/corpus_summary/S.dat")
vt.dump(datadir+"/corpus_summary/Vt.dat")
# vt = load(datadir+"/corpus_summary/Vt.dat")

# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# elements = random.choice(dim,(1000,1),replace=True) # pick 1000 random points to plot
# x = u[elements,0]*sqrt(s[0])
# y = u[elements,1]*sqrt(s[1])
# z = u[elements,2]*sqrt(s[2])
# # x = u[:,0]*sqrt(s[0])
# # y = u[:,1]*sqrt(s[1])
# # z = u[:,2]*sqrt(s[2])

# ax.set_title(file_name)
# ax.set_xlabel('X Label')
# ax.set_ylabel('Y Label')
# ax.set_zlabel('Z Label')
# ax.set_xlim([-0.1,+0.1])
# ax.set_ylim([-0.1,+0.1])
# ax.set_zlim([0.0,+0.1])
# ax.scatter(x, y, z)

# # fig = plt.figure()
# # ax = fig.add_subplot(111)
# # # x = u[:,0]*sqrt(s[0])
# # # y = u[:,1]*sqrt(s[1])
# # # z = u[:,2]*sqrt(s[2])
# # # elements = random.choice(dim,(10000,1),replace=True) # pick 1000 random points to plot
# # # x = u[elements,0]*sqrt(s[0])
# # # y = u[elements,1]*sqrt(s[1])
# # x = u[:,0]*sqrt(s[0])
# # y = u[:,1]*sqrt(s[1])
# # ax.set_title(file_name)
# # ax.set_xlabel('X Label')
# # ax.set_ylabel('Y Label')
# # ax.scatter(x, y)


# # # import seaborn as sns; sns.set(color_codes=True)
# # # ax = sns.regplot(x,y,fit_reg=False)
# # # x=range(rank)
# # # ax = sns.regplot(array(x),array(s),fit_reg=False)

# plt.show()
# # # plt.savefig(file_name+'.pdf')
