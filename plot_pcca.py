import pyemma
import os
import matplotlib.pyplot as plt
import numpy as np
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from pyemma import config

X=np.load('Y.npy',allow_pickle=True)
X = list(X)

M=pyemma.load('M.h5')

nstates = 4
M.pcca(nstates)

dtrajs=np.load('dtrajs.npy',allow_pickle=True)

metastable_traj = M.metastable_assignments[np.concatenate(dtrajs)]

plt.figure(figsize=(10,8))
pyemma.plots.plot_state_map(
    *X_concatenated[:, :2].T, metastable_traj)
#plt.xlabel('TIC1')
#plt.ylabel('TIC2')
#plt.show()
plt.savefig('pcca.png',dpi=110)
