import pyemma
import os
import matplotlib.pyplot as plt
import numpy as np
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from pyemma import config

M=pyemma.load('M.h5')

nstates = 4
M.pcca(nstates)

import pyemma
import os
import matplotlib.pyplot as plt
import numpy as np
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from pyemma import config

M=pyemma.load('M.h5')

nstates = 4
M.pcca(nstates)

mfpt = np.zeros((nstates, nstates))
for i in range(nstates):
    for j in range(nstates):
        mfpt[i, j] = (M.mfpt(
            M.metastable_sets[i],
            M.metastable_sets[j])*2*(10**(-6)))

inverse_mfpt = np.zeros_like(mfpt)
nz = mfpt.nonzero()
inverse_mfpt[nz] = 1.0 / mfpt[nz]


plt.figure(figsize=(10,8))
pos=np.array([[-1,0], [0, -1], [1, 0], [0, 0]])
pyemma.plots.plot_network(inverse_mfpt,arrow_label_format='%.3f \u03BCs',arrow_labels=mfpt,size=12)
plt.savefig('network_mfpt.png',dpi=210)


for i, s in enumerate(M.metastable_sets):
    print('SP_{} = {:f}'.format(i + 1, M.pi[s].sum()))



for i in range(0,len(mfpt)):
        for j in range(i+1,len(mfpt)):
            if (i==j) :
                continue
            else:

                print(i+1,' --->',j+1,':', round(mfpt[i][j],3))
                print(j+1,' --->',i+1,':', round(mfpt[j][i],3))
