#! /usr/bin/python

from scipy.sparse.linalg import svds
from scipy.sparse import coo_matrix
import os.path
import sys
from numpy import *

 
file_name = str(sys.argv[1])
dim = int(sys.argv[2])
# k = int(sys.argv[3])
# rank = int(sys.argv[4])

r = []
c = []
d = []
with open(file_name, 'rb') as infile:
    for line in infile:
        lista = line.strip().split('\t')
        r.append(int(lista[3])-1) # minus 1 because python counts from 0
        c.append(int(lista[4])-1)
        d.append(int(lista[2]))

coomat = coo_matrix((d, (r, c)), shape=(dim, dim),dtype=int) # sparse mat in coordinate format
coomat = coomat + transpose(coomat) 
# !!!! NEED TO DIVIDE BY 2 THE DIAGONAL, OTHERWISE DOUBLE COUNTING !!!!!
# !!! NEED TO NORMALIZE AND TAKE LOGARITHM
print(coomat)
# coomat.data = ma.log(0.5*coomat.nnz*coomat.data) # SPPMI
# coomat.data = ma.masked_less(coomat.data, 0)

# [u,s,vt] = svds(coomat, rank, which='LM', return_singular_vectors=True)

# savetxt(file_name + '_U.txt', u, fmt='%1.4e')
# savetxt(file_name + '_S.txt', s, fmt='%1.4e')
# savetxt(file_name + '_Vt.txt', transpose(vt), fmt='%1.4e') # useless because symm

# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# # fig = plt.figure()
# # ax = fig.add_subplot(111, projection='3d')
# # elements = random.choice(dim,(1000,1),replace=True) # pick 1000 random points to plot
# # x = u[elements,0]*sqrt(s[0])
# # y = u[elements,1]*sqrt(s[1])
# # z = u[elements,2]*sqrt(s[2])
# # ax.set_title(file_name)
# # ax.set_xlabel('X Label')
# # ax.set_ylabel('Y Label')
# # ax.set_zlabel('Z Label')
# # ax.scatter(x, y)

# fig = plt.figure()
# ax = fig.add_subplot(111)
# # x = u[:,0]*sqrt(s[0])
# # y = u[:,1]*sqrt(s[1])
# # z = u[:,2]*sqrt(s[2])
# elements = random.choice(dim,(10000,1),replace=True) # pick 1000 random points to plot
# # x = u[elements,0]*sqrt(s[0])
# # y = u[elements,1]*sqrt(s[1])
# x = u[:,0]*sqrt(s[0])
# y = u[:,1]*sqrt(s[1])
# ax.set_title(file_name)
# ax.set_xlabel('X Label')
# ax.set_ylabel('Y Label')
# ax.scatter(x, y)


# # import seaborn as sns; sns.set(color_codes=True)
# # ax = sns.regplot(x,y,fit_reg=False)
# # x=range(rank)
# # ax = sns.regplot(array(x),array(s),fit_reg=False)

# plt.show()
# # plt.savefig(file_name+'.pdf')
