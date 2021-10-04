import pyemma
import os
import matplotlib.pyplot as plt
import numpy as np
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from pyemma import config

X = np.load('X.npy')
data=list(X)



state_no = 1; state_index = np.int32(state_no); # Stae of interest
x_min = 0.285; x_max = 0.290; y_min = 0.20; y_max = 0.21; # Boundary values 
X_min = x_min; X_max = x_max; Y_min = y_min; Y_max = y_max;
cnt = 0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            cnt = cnt+1

print(cnt)

state = np.zeros((cnt,4))
var=0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            state[var,0] = data[i][j][0]
            state[var,1] = data[i][j][1]
            state[var,2] = np.int32(i)
            state[var,3] = np.int32(j)
            var = var+1

print(np.min(state[:,2]),np.max(state[:,2]),np.min(state[:,3]),np.max(state[:,3]))
np.savetxt('./state_'+str(state_index)+'.txt',state[::10,0:2],fmt='%16.8lf')

samples_pre = state[::1000,2:].astype('int32')
np.savetxt('./samples_1.txt',samples_pre,fmt='%6d')

pyemma.coordinates.save_trajs(inp,samples,outfiles=['./state_{}.xtc'.format(n,n) for n in {state_index,state_index}]) # Generate inp as mentioned in cluster_extraction.py 
