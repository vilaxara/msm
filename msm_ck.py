import pyemma
import os
import matplotlib.pyplot as plt
import numpy as np
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from pyemma import config

dtrajs=np.load('./dtrajs.npy',allow_pickle=True)
dtrajs=list(dtrajs)

lag=10
n_states=5

M = msm.estimate_markov_model(dtrajs, lag=lag)
#M=msm.bayesian_markov_model(dtrajs, lag=lag)
M.save('M_10.h5',overwrite=True)

plt.figure(figsize=(10,8))
mplt.plot_cktest(M.cktest(n_states))
plt.savefig('./ck10_5.png',dpi=110)
plt.clf()


print('fraction of states used = ', M.active_state_fraction)
print('fraction of counts used = ', M.active_count_fraction)


np.savetxt('S1.txt',M.metastable_sets[0]) 
# To save the clusters in a perticular metastable state. Use awk.sh to sort the clusters according to the populations and cluster_extractiion.py to extract the perticluar cluster of interest
