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

p_res=138 # NUmber of protein amino-acid residues

ind_arr = np.zeros((p_res,2))

for i in range(0,138) :

        ind_arr[i][0]=138
        ind_arr[i][1]=i


feat.add_residue_mindist(residue_pairs=ind_arr, scheme='closest-heavy', threshold=0.5)

print(feat.describe())

traj_list=[] # Create a list of xtc file locations

inp = coor.source(traj_list, feat)
X = inp.get_output()


# Dimensional reduction
lag = 100  # tica lagtime
tica_obj = coor.tica(X, lag4=tica_lag)
print ('Retained dimension: ', tica_obj.dimension())
Y = tica_obj.get_output()

np.save('Y',Y)

plt.figure(figsize=(10,8))
mplt.plot_free_energy(np.vstack(Y)[:, 0], np.vstack(Y)[:, 1])
plt.xlabel('independent component 1'); plt.ylabel('independent component 2')
plt.Text(0, 0.5, 'independent component 2')
plt.savefig("t_all_fes100l.png")
plt.clf()
