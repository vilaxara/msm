import pyemma
import os
import matplotlib.pyplot as plt
import numpy as np
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from pyemma import config

'''

To compute netwrok path flux

'''
M = pyemma.load('M_250_1ns.h5')
M.pcca(3)
A=M.metastable_sets[1]
B=M.metastable_sets[2]
tpt=msm.tpt(M,A,B)
cg, cgflux = tpt.coarse_grain(M.metastable_sets)

paths, path_fluxes = cgflux.pathways(fraction=0.99)
print('percentage       \tpath')
print('-------------------------------------')
for i in range(len(paths)):
    print(np.round(path_fluxes[i] / np.sum(path_fluxes), 3),' \t', paths[i] + 1)


(paths,pathfluxes) = cgflux.pathways()
cumflux = 0
print ("Path flux\t\t%path\t%of total\tpath")
for i in range(len(paths)):
    cumflux += pathfluxes[i]
    print (pathfluxes[i],'\t','%3.1f'%(100.0*pathfluxes[i]/tpt.total_flux),'%\t','%3.1f'%(100.0*cumflux/tpt.total_flux),'%\t\t',paths[i])
