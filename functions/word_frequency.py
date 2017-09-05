#!/usr/bin/env python

import os
import sys
import numpy as np

filename = str(sys.argv[1])
path2file = os.path.dirname(filename)

f = open(filename,"r")
lines = f.readlines()
counts = []
for x in lines:
    counts.append(int(x.split('\t')[0]))
f.close()

import collections
import pandas as pd
counter = collections.Counter(counts)
df = pd.DataFrame.from_dict(counter, orient='index').reset_index()
df.columns = ['frequency','multiplicity']
xmin = df.frequency.quantile(0.0)
xmax = df.frequency.quantile(1.0)
print "xmin = " + str(xmin)
print "xmax = " + str(xmax)
filtered_df = df.loc[df['frequency'] <= xmax]
filtered_df = filtered_df.loc[filtered_df['frequency'] >= xmin]
filtered_df['multiplicity'] = filtered_df.multiplicity/filtered_df.multiplicity.sum() # normalize the y axis
filtered_df = filtered_df.sort_values(by='frequency')

df.to_csv(os.path.join(path2file, "freq-multiplicity" + "." + "csv"))
######
# To estimate scaling exponent of the cumulative frequency plot use formula (5) from "Power laws, Pareto distributions and Zipfs laws" by Newman
######
alpha = 1.0+filtered_df.shape[0]*1.0/np.log10(filtered_df.frequency/xmin).sum()
C = alpha - 1.0
print "alpha = " + str(alpha)
print "C = " + str(C)
import matplotlib.pyplot as plt
plt.plot(filtered_df.frequency,C*filtered_df.frequency**(-alpha),'r--')
# plt.plot(filtered_df.frequency,C*filtered_df.frequency**(-alpha+0.5),'b--')
# plt.plot(filtered_df.frequency,C*filtered_df.frequency**(-alpha-0.5),'g--')
plt.plot(filtered_df.frequency,filtered_df.multiplicity,'ro') 
plt.xscale('log')# the log is base 10
plt.yscale('log')
plt.xlabel('Frequency')
plt.ylabel('Probability')
# plt.show()
plt.savefig(os.path.join(path2file, "freq-prob.pdf"))
