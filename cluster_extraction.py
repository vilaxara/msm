import pyemma
import os
import matplotlib.pyplot as plt
import numpy as np
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from pyemma import config

topfile='*.pdb'
feat = coor.featurizer(topfile)

traj_list=[] # Create a list of xtc file locations

inp = coor.source(traj_list, feat)
#X = inp.get_output()

dtrajs = np.load('dtrajs.npy',allow_pickle=True)
dtrajs = list(dtrajs)

# Cluster of interest

fr19=[]
for i in range(0,len(dtrajs)):
    a=np.where(dtrajs[i]==19)
    if np.shape(a)[1]==0:
        continue
    else:
        for j in list(a[0]):
            fr19.append([i,j])


coor.save_trajs(inp, np.array(np.array(fr19)), outfiles=['./C_19.xtc'])
