import pyemma
import os
import matplotlib.pyplot as plt
import numpy as np
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from pyemma import config

'''

Script to divide the data into clusters. We use k-means clustering. If your data is high dimensional data, use t-ICA to reduce the dimensions, else proceed with clustering.

'''

X =np.load('Y.npy',allow_pickle=True)
X = list(X)

n_clusters = 500 # number of clusters

clustering = coor.cluster_kmeans(X, k=n_clusters, max_iter=1000, tolerance=1e-10, fixed_seed=True)

clustering.save('500C.h5')

dtrajs = clustering.dtrajs
np.save('dtrajs', dtrajs)

# Cluster population

histogram = np.bincount(np.concatenate(dtrajs), minlength=len(cluster.clustercenters))
np.savetxt('cluster_counts.txt',histogram,fmt='%d')

# Plotting clusters on top of free-energy surface

mplt.plot_free_energy(np.vstack(X)[:,0], np.vstack(X)[:,1])
cc_x=clustering.clustercenters[:,0]
cc_y=clustering.clustercenters[:,1]
plt.plot(cc_x,cc_y, linewidth=0, marker='o', markersize=3, color='black')
plt.savefig("fes_cc.png")
plt.clf()

