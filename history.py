156/1: import pyemma
156/2: its = pyemma.load("its_macro_five_states.h5")
156/3: model = its.models[1]
156/4: model.stationary_distribution
157/1: import pyemma
157/2: its = pyemma.load("its_macro_five_states.h5")
157/3: model = its.models[1]
157/4: model
157/5: model.stationary_distribution
157/6: model.nstates
157/7: model.dtrajs_full
157/8: model.dtrajs_full.shape
157/9: len(model.dtrajs_full)
157/10: model.dtrajs_full[0].shape
157/11: model.dtrajs_full[419].shape
157/12: model.dtrajs_full[1].shape
157/13: model.dtrajs_full[2].shape
157/14: import numpy as np
157/15: np.where(model.dtrajs_full[0]==1)
157/16: np.where(model.dtrajs_full[1]==1)
157/17: np.where(model.dtrajs_full[2]==1)
157/18: model.active_state_indexes
158/1: import pyemma
158/2: model = pyemma.load("M_new_5ns.h5")
158/3: model.active_state_indexes
158/4: model.active_set
159/1: import pyemma
159/2: its = pyemma.load("../its2.h5")
159/3: its.models[15]
159/4: its.models[14]
160/1: import pyemma
160/2: its = pyemma.load("its_macro_five_states.h5")
160/3: model = its.models[1]
160/4: model
160/5: model = its.models[2]
160/6: model = its.models[0]
160/7: model
160/8: model = its.models[1]
160/9: model.stationary_distribution
160/10: import msmtools
160/11: msmtools.analysis.mfpt?
160/12: msmtools.analysis.mfpt(model.transition_matrix, 0)
160/13: msmtools.analysis.mfpt(model.transition_matrix, 1)
160/14: msmtools.analysis.mfpt(model.transition_matrix, 2)
160/15: msmtools.analysis.mfpt(model.transition_matrix, 3)
160/16: msmtools.analysis.mfpt(model.transition_matrix, 4)
160/17: model.mfpt?
160/18: model.mfpt(0, 1)
160/19: model.mfpt([0],[1])*10*1e-12*1e9
160/20: model.mfpt([0],[2])*10*1e-12*1e9
150/314: pwd
150/315: ls
150/316:
#fig, pos = mplt.plot_markov_model(M_new, pos=pos, state_sizes=state_sizes
#plt.savefig('network.pdf')
150/317: pos
161/1:
import pyemma                                                                       
import os                                                                           
import matplotlib.pyplot as plt                                                    
import numpy as np                                                                 
import pyemma.coordinates as coor 
import pyemma.msm as msm 
import pyemma.plots as mplt 
from pyemma import config
161/2: its = pyemma.load('../its.h5')
161/3: dtrajs = np.load('../dtrajs.npy')
161/4: dtrajs = np.load('../dtrajs.npy',allow_pickle=True)
161/5: dtrajs = list(dtrajs)
161/6: M = msm.estimate_markov_model(dtrajs,2000)
161/7:
print('fraction of states used = ', M.active_state_fraction)
print('fraction of counts used = ', M.active_count_fraction)
161/8:
dtrajs_concatenated = np.concatenate(dtrajs)
data = np.load('../../X7.npy')
data = list(data)
data_concatenated = np.concatenate(data)

nstates = 5
M.pcca(nstates)
metastable_traj = M.metastable_assignments[dtrajs_concatenated]

pyemma.plots.plot_state_map(
    *data_concatenated[:, :2].T, metastable_traj)
plt.xlabel('RMSD 1')
plt.ylabel('RMSD 2')
plt.show()
for i, s in enumerate(M.metastable_sets):
    print('π_{} = {:f}'.format(i + 1, M.pi[s].sum()));plt.savefig('5state.png')
161/9:
dtrajs_concatenated = np.concatenate(dtrajs)
data = np.load('../../X7.npy',allow_pickle=True)
data = list(data)
data_concatenated = np.concatenate(data)

nstates = 5
M.pcca(nstates)
metastable_traj = M.metastable_assignments[dtrajs_concatenated]

pyemma.plots.plot_state_map(
    *data_concatenated[:, :2].T, metastable_traj)
plt.xlabel('RMSD 1')
plt.ylabel('RMSD 2')
plt.show()
for i, s in enumerate(M.metastable_sets):
    print('π_{} = {:f}'.format(i + 1, M.pi[s].sum()));plt.savefig('5state.png')
161/10:
j1=0;j2=0
for i in range(np.shape(M.eigenvectors_right())[0]):
    if(M.eigenvectors_right()[i,1] > 0.0):
        j1=j1+1

    if(M.eigenvectors_right()[i,1] < 0.0):
        j2=j2+1
print(j1,j2)
Mp = np.zeros((j1,4))
Mn = np.zeros((j2,4))
j1=0;j2=0
for i in range(np.shape(M.eigenvectors_right())[0]):
    if(M.eigenvectors_right()[i,1] > 0.0):
        Mp[j1][0] = i+1
        Mp[j1][1] = i
        Mp[j1][2] = M.eigenvectors_right()[i,1]
        Mp[j1][3] = M.eigenvectors_left()[0,i]
        j1=j1+1
    if(M.eigenvectors_right()[i,1] < 0.0):
        Mn[j2][0] = i+1
        Mn[j2][1] = i
        Mn[j2][2] = M.eigenvectors_right()[i,1]
        Mn[j2][3] = M.eigenvectors_left()[0,i]
        j2=j2+1

print("p:  ",np.sum(Mp[:,3])*100 )
print("n:  ",np.sum(Mn[:,3])*100 )
161/11:
j1=0;j2=0
for i in range(np.shape(M.eigenvectors_right())[0]):
    if(M.eigenvectors_right()[i,1] > 0.0) and (M.eigenvectors_right()[i,2] > 0.0):
        j1=j1+1

    if(M.eigenvectors_right()[i,1] > 0.0) and (M.eigenvectors_right()[i,2] < 0.0):
        j2=j2+1
print(j1,j2)
Mpp = np.zeros((j1,4))
Mpn = np.zeros((j2,4))
j1=0;j2=0
for i in range(np.shape(M.eigenvectors_right())[0]):
    if(M.eigenvectors_right()[i,1] > 0.0) and (M.eigenvectors_right()[i,2] > 0.0):
        Mpp[j1][0] = i+1
        Mpp[j1][1] = i
        Mpp[j1][2] = M.eigenvectors_right()[i,1]
        Mpp[j1][3] = M.eigenvectors_left()[0,i]
        j1=j1+1
    if(M.eigenvectors_right()[i,1] > 0.0) and (M.eigenvectors_right()[i,2] < 0.0):
        Mpn[j2][0] = i+1
        Mpn[j2][1] = i
        Mpn[j2][2] = M.eigenvectors_right()[i,1]
        Mpn[j2][3] = M.eigenvectors_left()[0,i]
        j2=j2+1
print("pp:  ",np.sum(Mpp[:,3])*100 )
print("pn:  ",np.sum(Mpn[:,3])*100 )
161/12:

j1=0;j2=0
for i in range(np.shape(M.eigenvectors_right())[0]):
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] > 0.0):
        j1=j1+1

    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] < 0.0):
        j2=j2+1

print(j1,j2)
Mnp = np.zeros((j1,4))
Mnn = np.zeros((j2,4))

j1=0;j2=0
for i in range(np.shape(M.eigenvectors_right())[0]):
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] > 0.0):
        Mnp[j1][0] = i+1
        Mnp[j1][1] = i
        Mnp[j1][2] = M.eigenvectors_right()[i,1]
        Mnp[j1][3] = M.eigenvectors_left()[0,i]
        j1=j1+1
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] < 0.0):
        Mnn[j2][0] = i+1
        Mnn[j2][1] = i
        Mnn[j2][2] = M.eigenvectors_right()[i,1]
        Mnn[j2][3] = M.eigenvectors_left()[0,i]
        j2=j2+1

print("np:  ",np.sum(Mnp[:,3])*100 )
print("nn:  ",np.sum(Mnn[:,3])*100 )
161/13:
np.savetxt('pp.txt',Mpp[:,1],'%6.0d')
np.savetxt('pn.txt',Mpn[:,1],'%6.0d')
np.savetxt('np.txt',Mnp[:,1],'%6.0d')
np.savetxt('nn.txt',Mnn[:,1],'%6.0d')
161/14:
j1=0;j2=0
for i in range(np.shape(M.eigenvectors_right())[0]):
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] > 0.0) and (M.eigenvectors_right()[i,3] > 0.0):
        j1=j1+1

    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] > 0.0) and (M.eigenvectors_right()[i,3] < 0.0):
        j2=j2+1
print(j1,j2)
Mnpp = np.zeros((j1,4))
Mnpn = np.zeros((j2,4))

j1=0;j2=0
for i in range(np.shape(M.eigenvectors_right())[0]):
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] > 0.0) and (M.eigenvectors_right()[i,3] > 0.0):
        Mnpp[j1][0] = i+1
        Mnpp[j1][1] = i
        Mnpp[j1][2] = M.eigenvectors_right()[i,1]
        Mnpp[j1][3] = M.eigenvectors_left()[0,i]
        j1=j1+1
    if(M.eigenvectors_right()[i,1] < 0.0) and (M.eigenvectors_right()[i,2] > 0.0) and (M.eigenvectors_right()[i,3] < 0.0):
        Mnpn[j2][0] = i+1
        Mnpn[j2][1] = i
        Mnpn[j2][2] = M.eigenvectors_right()[i,1]
        Mnpn[j2][3] = M.eigenvectors_left()[0,i]
        j2=j2+1

print("npp:  ",np.sum(Mnpp[:,3])*100 )
print("npn:  ",np.sum(Mnpn[:,3])*100 )
161/15:
np.savetxt('npp.txt',Mnpp[:,1],'%6.0d')
np.savetxt('npn.txt',Mnpn[:,1],'%6.0d')
161/16: M
161/17: dtrajs2 = np.load('dtrajs_2d.npy',allow_pickle=True)
161/18: M_new = msm.estimate_markov_model(dtrajs2, 2000)
161/19: dtrajs2 = list(dtrajs2)
161/20: M_new = msm.estimate_markov_model(dtrajs2, 2000)
161/21: M_new
161/22: M_new.stationary_distribution
161/23: np.shape(M.metastable_sets)
161/24: M.metastable_sets = []
161/25: M.metastable_sets[0] = []
161/26: pp = np.loadtxt('./pp.txt',dtype='int64')
161/27: np.shape(M.metastable_sets)
161/28: np.shape(M.metastable_sets[0])
161/29: np.shape(M.metastable_sets[1])
161/30: M.metastable_sets[1] = []
161/31: M.metastable_sets[2] = []
161/32: M.metastable_sets[3] = []
161/33: M.metastable_sets[4] = []
161/34: pn = np.loadtxt('./pn.txt',dtype='int64')
161/35: npp = np.loadtxt('./npp.txt',dtype='int64')
161/36: npn = np.loadtxt('./npn.txt',dtype='int64')
161/37: nn = np.loadtxt('./nn.txt',dtype='int64')
161/38: M.metastable_sets[0] = pp
161/39: M.metastable_sets[1] = pn
161/40: M.metastable_sets[2] = npp
161/41: M.metastable_sets[3] = npn
161/42: M.metastable_sets[4] = nn
161/43: np.shape(M.metastable_sets)
161/44: np.shape(M.metastable_sets[0])
161/45: np.shape(M.metastable_sets[1])
161/46: np.shape(M.metastable_sets[2])
161/47: np.shape(M.metastable_sets[3])
161/48: np.shape(M.metastable_sets[4])
161/49: s5 = np.loadtxt('./5s.txt',dtype='int64')
161/50: metastable_traj = []
161/51: metastable_traj = s5
161/52:
pyemma.plots.plot_state_map( 
    *data_concatenated[:, :2].T, metastable_traj) 
plt.xlabel('RMSD 1') 
plt.ylabel('RMSD 2') 
plt.show() 
for i, s in enumerate(M.metastable_sets): 
    print('π_{} = {:f}'.format(i + 1, M.pi[s].sum()));plt.savefig('5state_new.png')
161/53:
states = [c1, c2, c3, c4, c5]
for i in [1,2,3,4,5]:
    for j in [1,2,3,4,5]:
        if (i!=j):
            print('c'+str(i),"-->", 'c'+str(j), int(M_new.mfpt(states[i-1], states[j-1]))*10*0.001)
161/54: c1 = 0
161/55: c2 = 1
161/56: c3 = 2
161/57: c4 = 3
161/58: c5 = 5
161/59:
states = [c1, c2, c3, c4, c5]
for i in [1,2,3,4,5]:
    for j in [1,2,3,4,5]:
        if (i!=j):
            print('c'+str(i),"-->", 'c'+str(j), int(M_new.mfpt(states[i-1], states[j-1]))*10*0.001)
161/60:
for i in [1,2,3,4,5]:
    for j in [1,2,3,4,5]:
        if (i!=j):
            print('c'+str(i),"-->", 'c'+str(j), int(M_new.mfpt(M.metastable_sets[i-1], M.metastable_sets[j-1]))*10*0.001)
161/61:
nstates = 5
mfpt = np.zeros((nstates, nstates)) 
for i in range(nstates): 
    for j in range(nstates): 
        mfpt[i, j] = 10*1e-12*M_new.mfpt( 
            M.metastable_sets[i], 
            M.metastable_sets[j]) 
             
        print(i,j,mfpt[i,j])
161/62:
nstates = 5
mfpt = np.zeros((nstates, nstates)) 
for i in range(nstates): 
    for j in range(nstates): 
        mfpt[i, j] = 10*0.001*M.mfpt( 
            M.metastable_sets[i], 
            M.metastable_sets[j]) 
             
        print(i,j,mfpt[i,j])
161/63:
pyemma.plots.plot_state_map( 
    *data_concatenated[:, :2].T, metastable_traj) 
plt.xlabel('RMSD 1') 
plt.ylabel('RMSD 2') 
plt.show() 
for i, s in enumerate(M.metastable_sets): 
    print('π_{} = {:f}'.format(i + 1, M.pi[s].sum()));plt.savefig('5state_new.pdf')
161/64:
clustering = pyemma.load('../clustering.h5')
highest_membership = M.metastable_distributions.argmax(1)
coarse_state_centers5 = clustering.clustercenters[M.active_set[highest_membership]]
161/65: coarse_state_centers5
161/66: M_new.metastable_distributions
161/67: coarse_state_centers5
161/68: pos = [[0.27337772, 0.2001345 ],[0.11118733, 0.10308006],[0.4887169 , 0.1104988 ],[0.5556613 , 0.36805862],[0.77665424, 0.60384786]]
161/69: M_new.stationary_distribution
161/70: state_sizes = M_new.stationary_distribution
161/71:
fig, pos = mplt.plot_markov_model(M_new, pos=pos, state_sizes=state_sizes)
plt.savefig('network.pdf')
161/72: state_sizes = list(state_sizes)
161/73:
fig, pos = mplt.plot_markov_model(M_new, pos=pos, state_sizes=state_sizes)
plt.savefig('network.pdf')
161/74: pos
161/75: np.shape(pos)
161/76: pos = list(pos)
161/77:
fig, pos = mplt.plot_markov_model(M_new, pos=pos, state_sizes=state_sizes)
plt.savefig('network.pdf')
161/78: pos
161/79: state_sizes
161/80: state_sizes = [[0.1159849201174332],[0.6104111856351302],[0.16435720632421782],[0.06683501303412505],[0.04241167488909371]]
161/81: np.shape(state_sizes)
161/82:
fig, pos = mplt.plot_markov_model(M_new, pos=pos, state_sizes=state_sizes)
plt.savefig('network.pdf')
161/83: state_sizes = M_new.stationary_distribution
161/84:
print('probabilities', M_new.stationary_distribution)
print('plotting sizes', state_sizes)
161/85:
fig, pos = mplt.plot_markov_model(M_new, pos=pos, state_sizes=state_sizes)
plt.savefig('network.pdf')
161/86: pos = list(pos)
161/87:
fig, pos = mplt.plot_markov_model(M_new, pos=pos, state_sizes=state_sizes)
plt.savefig('network.pdf')
161/88: pos
161/89:
fig, pos = mplt.plot_markov_model(M_new, pos=pos)
plt.savefig('network.png')
161/90:
fig, pos = mplt.plot_markov_model(M_new)
plt.savefig('network.png')
161/91: np.shape(state_sizes)
161/92: np.shape(pos)
161/93: pos
161/94: coarse_state_centers5
161/95: np.shape(coarse_state_centers5)
161/96:
fig, pos = mplt.plot_markov_model(M_new)
plt.savefig('network.pdf')
161/97:
fig, pos = mplt.plot_markov_model(M_new,state_sizes=state_sizes)
plt.savefig('network.png')
161/98:
fig, pos = mplt.plot_markov_model(M_new,state_sizes=state_sizes)
plt.savefig('network.pdf')
161/99:
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2) 
plt.ylim(0, )
plt.xlim(1,5000 )
plt.savefig('./its_new.png')
161/100:
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2) 
plt.ylim(0, )
plt.xlim(1,500 )
plt.savefig('./its_new.png')
161/101:
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2) 
plt.ylim(0, )
plt.xlim(0,3000 )
plt.savefig('./its_new2.png')
161/102: plt.clf()
161/103:
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2) 
plt.ylim(0, )
plt.xlim(0,3000 )
plt.savefig('./its_new2.png')
161/104: plt.clf()
161/105:
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2) 
plt.ylim(0, )
plt.xlim(0,5000 )
plt.savefig('./its_new.png')
161/106:
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2) 
plt.ylim(0, )
plt.xlim(0,5000 );plt.figaspect(640*480)
plt.savefig('./its_new.png')
161/107:
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2) 
plt.ylim(0, )
plt.xlim(0,5000 );plt.figure(figsize=640*480)
plt.savefig('./its_new.png')
161/108:
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2) 
plt.ylim(0, )
plt.xlim(0,5000 )
plt.figure(figsize=640*480)
plt.savefig('./its_new.png')
161/109:
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2) 
plt.ylim(0, )
plt.xlim(0,5000 )
plt.show()
161/110:
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2) 
plt.ylim(0, )
plt.xlim(10,5000 )
plt.savefig('its_new.png')
161/111:
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2) 
plt.ylim(0, )
plt.xlim(10,5000 )
plt.figure(figsize=(640,480))
plt.savefig('its_new.png')
162/1:
import pyemma                                                                       
import os                                                                           
import matplotlib.pyplot as plt                                                    
import numpy as np                                                                 
import pyemma.coordinates as coor 
import pyemma.msm as msm 
import pyemma.plots as mplt 
from pyemma import config
162/2: its = pyemma.load('../its.h5')
162/3:
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2) 
plt.ylim(0, )
plt.xlim(1,5000 )
plt.savefig('its_new.png')
162/4:
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2) 
plt.ylim(0, )
plt.xlim(1,3000 )
plt.savefig('its_new2.png')
162/5:
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2) 
plt.ylim(0, )
plt.xlim(1,4000 )
plt.savefig('its_new3.png')
162/6: X = np.load('../../X7.npy')
162/7: X = np.load('../../X7.npy',allow_pickle=True)
162/8: X = list(X)
162/9: dtrajs = np.load('../dtrajs.npy')
162/10: dtrajs = np.load('../dtrajs.npy',allow_pickle=True)
162/11: dtrajs = list(dtrajs)
162/12: dtrajs_concatenated = np.concatenate(dtrajs)
162/13: M_new = msm.estimate_markov_model(dtrajs, 2000)
162/14:
pyemma.plots.plot_free_energy(*data_concatenated.T, weights=np.concatenate(M.trajectory_weights()), ax=axs[0],kT=0.593,cbar_label='kcal')
axs[1].set_xlabel('RMSD 1')
axs[1].set_ylabel('RMSD 2')
axs[1].set_title('Reweighted free energy')
plt.savefig('fes_reweight.png')
162/15: data_concatenated = np.concatenate(X)
162/16:
pyemma.plots.plot_free_energy(*data_concatenated.T, weights=np.concatenate(M.trajectory_weights()), ax=axs[0],kT=0.593,cbar_label='kcal')
axs[1].set_xlabel('RMSD 1')
axs[1].set_ylabel('RMSD 2')
axs[1].set_title('Reweighted free energy')
plt.savefig('fes_reweight.png')
162/17:
pyemma.plots.plot_free_energy(*data_concatenated.T, weights=np.concatenate(M_new.trajectory_weights()), ax=axs[0],kT=0.593,cbar_label='kcal')
axs[1].set_xlabel('RMSD 1')
axs[1].set_ylabel('RMSD 2')
axs[1].set_title('Reweighted free energy')
plt.savefig('fes_reweight.png')
162/18:
pyemma.plots.plot_free_energy(*data_concatenated.T, weights=np.concatenate(M_new.trajectory_weights()),kT=0.593,cbar_label='kcal')
axs[1].set_xlabel('RMSD 1')
axs[1].set_ylabel('RMSD 2')
axs[1].set_title('Reweighted free energy')
plt.savefig('fes_reweight.png')
162/19:
pyemma.plots.plot_free_energy(*data_concatenated.T, weights=np.concatenate(M_new.trajectory_weights()),kT=0.593,cbar_label='kcal')
plt.Text(0, 0.5, 'independent componen 2')
plt.xlabel('RMSD 1')
plt.ylabel('RMSD 2')
plt.savefig('fes_reweight.png')
162/20:
pyemma.plots.plot_free_energy(*data_concatenated.T, weights=np.concatenate(M_new.trajectory_weights()),kT=0.593,cbar_label='kcal',cmap='jet_r')
plt.Text(0, 0.5, 'independent componen 2')
plt.xlabel('RMSD 1')
plt.ylabel('RMSD 2')
plt.savefig('fes_reweight.pdf')
162/21:
pyemma.plots.plot_cktest(M_new.cktest(5), units='ps');
plt.savefig('ck5.png')
162/22:
pyemma.plots.plot_cktest(bayesian_msm.cktest(5), units='ps');
plt.savefig('ck5_bayes.png')
162/23: M_new.save('M_new_2000cls_20nslagtime.h5')
162/24:
traj_list = [] 
indir = '../../../data5/XTC2/'  
topfile =  indir+'../ref_1phc_backbone.pdb'  
for i in range(1,10): 
    file = '../../../data5/XTC2/no_pbc2_0'+str(i)+'.xtc' 
    print (file) 
    traj_list.append(file)

for i in range(10,21): 
    file = '../../../data5/XTC2/no_pbc2_'+str(i)+'.xtc' 
    print (file) 
    traj_list.append(file)

for i in range(0,10): 
    file = '../../../data5/XTC2/no_pbc2_00'+str(i)+'.xtc' 
    print (file) 
    traj_list.append(file)

for i in range(10,100): 
    file = '../../../data5/XTC2/no_pbc_0'+str(i)+'.xtc' 
    print (file) 
    traj_list.append(file)

for i in range(100,400): 
    file = '../../../data5/XTC2/no_pbc2_'+str(i)+'.xtc' 
    print (file) 
    traj_list.append(file)
162/25: np.shape(traj_list)
162/26:
feat = coor.featurizer(topfile)
inp = coor.source(traj_list, feat)
162/27:
traj_list = [] 
indir = '../../../data5/XTC2/'  
topfile =  indir+'../ref_1phc_backbone.pdb'  
for i in range(1,10): 
    file = '../../../data5/XTC2/no_pbc2_0'+str(i)+'.xtc' 
    print (file) 
    traj_list.append(file)

for i in range(10,21): 
    file = '../../../data5/XTC2/no_pbc2_'+str(i)+'.xtc' 
    print (file) 
    traj_list.append(file)

for i in range(0,10): 
    file = '../../../data5/XTC2/no_pbc2_00'+str(i)+'.xtc' 
    print (file) 
    traj_list.append(file)

for i in range(10,100): 
    file = '../../../data5/XTC2/no_pbc2_0'+str(i)+'.xtc' 
    print (file) 
    traj_list.append(file)

for i in range(100,400): 
    file = '../../../data5/XTC2/no_pbc2_'+str(i)+'.xtc' 
    print (file) 
    traj_list.append(file)
162/28: np.shape(traj_list)
162/29:
feat = coor.featurizer(topfile)
inp = coor.source(traj_list, feat)
162/30: data = X
162/31: np.shape(data)
162/32:
state_no = 0; state_index = np.int32(state_no);
x_min = 0.1; x_max = 0.102; y_min = 0.05; y_max = 0.06; #State1 final
X_min = x_min; X_max = x_max; Y_min = y_min; Y_max = y_max;
cnt = 0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            cnt = cnt+1
162/33:
state_no = 0; state_index = np.int32(state_no);
x_min = 0.1; x_max = 0.102; y_min = 0.05; y_max = 0.06; #State1 final
X_min = x_min; X_max = x_max; Y_min = y_min; Y_max = y_max;
cnt = 0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            cnt = cnt+1
162/34:
state_no = 0; state_index = np.int32(state_no);
x_min = 0.11; x_max = 0.102; y_min = 0.05; y_max = 0.06; #State1 final
X_min = x_min; X_max = x_max; Y_min = y_min; Y_max = y_max;
cnt = 0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            cnt = cnt+1
            
print(cnt)
162/35:
state_no = 0; state_index = np.int32(state_no);
x_min = 0.1; x_max = 0.102; y_min = 0.05; y_max = 0.06; #State1 final
X_min = x_min; X_max = x_max; Y_min = y_min; Y_max = y_max;
cnt = 0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            cnt = cnt+1
            
print(cnt)
162/36:
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
162/37:
samples_pre = state[::1000,2:].astype('int32')
np.savetxt('./samples_0.txt',samples_pre,fmt='%6d')
samples = np.array([np.array(samples_pre)])
162/38: pyemma.coordinates.save_trajs(inp,samples,outfiles=['./state_{}.xtc'.format(n,n) for n in {state_index,state_index}])
162/39:
state_no = 0; state_index = np.int32(state_no);
x_min = 0.28; x_max = 0.29; y_min = 0.20; y_max = 0.22; #State1 final
X_min = x_min; X_max = x_max; Y_min = y_min; Y_max = y_max;
cnt = 0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            cnt = cnt+1
            
print(cnt)
162/40:
state_no = 0; state_index = np.int32(state_no);
x_min = 0.285; x_max = 0.290; y_min = 0.20; y_max = 0.21; #State1 final
X_min = x_min; X_max = x_max; Y_min = y_min; Y_max = y_max;
cnt = 0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            cnt = cnt+1
            
print(cnt)
162/41:
state_no = 1; state_index = np.int32(state_no);
x_min = 0.285; x_max = 0.290; y_min = 0.20; y_max = 0.21; #State1 final
X_min = x_min; X_max = x_max; Y_min = y_min; Y_max = y_max;
cnt = 0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            cnt = cnt+1
            
print(cnt)
162/42:
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
162/43:
samples_pre = state[::1000,2:].astype('int32')
np.savetxt('./samples_1.txt',samples_pre,fmt='%6d')
samples = np.array([np.array(samples_pre)])
162/44: pyemma.coordinates.save_trajs(inp,samples,outfiles=['./state_{}.xtc'.format(n,n) for n in {state_index,state_index}])
162/45:
state_no = 3; state_index = np.int32(state_no);
x_min = 0.45; x_max = 0.455; y_min = 0.10; y_max = 0.11; #State1 final
X_min = x_min; X_max = x_max; Y_min = y_min; Y_max = y_max;
cnt = 0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            cnt = cnt+1
            
print(cnt)
162/46:
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
162/47:
samples_pre = state[::1000,2:].astype('int32')
np.savetxt('./samples_1.txt',samples_pre,fmt='%6d')
samples = np.array([np.array(samples_pre)])
162/48: pyemma.coordinates.save_trajs(inp,samples,outfiles=['./state_{}.xtc'.format(n,n) for n in {state_index,state_index}])
162/49:
state_no = 0; state_index = np.int32(state_no);
x_min = 0.1; x_max = 0.102; y_min = 0.105; y_max = 0.106; #State1 final
X_min = x_min; X_max = x_max; Y_min = y_min; Y_max = y_max;
cnt = 0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            cnt = cnt+1
            
print(cnt)
162/50:
state_no = 0; state_index = np.int32(state_no);
x_min = 0.1; x_max = 0.102; y_min = 0.10; y_max = 0.105; #State1 final
X_min = x_min; X_max = x_max; Y_min = y_min; Y_max = y_max;
cnt = 0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            cnt = cnt+1
            
print(cnt)
162/51:
state_no = 0; state_index = np.int32(state_no);
x_min = 0.1; x_max = 0.102; y_min = 0.104; y_max = 0.106; #State1 final
X_min = x_min; X_max = x_max; Y_min = y_min; Y_max = y_max;
cnt = 0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            cnt = cnt+1
            
print(cnt)
162/52:
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
162/53:
samples_pre = state[::1000,2:].astype('int32')
np.savetxt('./samples_1.txt',samples_pre,fmt='%6d')
samples = np.array([np.array(samples_pre)])
162/54: pyemma.coordinates.save_trajs(inp,samples,outfiles=['./state_{}.xtc'.format(n,n) for n in {state_index,state_index}])
162/55:
state_no = 0; state_index = np.int32(state_no);
x_min = 0.05; x_max = 0.06; y_min = 0.10; y_max = 0.102; #State1 final
X_min = x_min; X_max = x_max; Y_min = y_min; Y_max = y_max;
cnt = 0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            cnt = cnt+1
            
print(cnt)
162/56:
state_no = 0; state_index = np.int32(state_no);
x_min = 0.05; x_max = 0.06; y_min = 0.10; y_max = 0.105; #State1 final
X_min = x_min; X_max = x_max; Y_min = y_min; Y_max = y_max;
cnt = 0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            cnt = cnt+1
            
print(cnt)
162/57:
state_no = 0; state_index = np.int32(state_no);
x_min = 0.05; x_max = 0.06; y_min = 0.07; y_max = 0.08; #State1 final
X_min = x_min; X_max = x_max; Y_min = y_min; Y_max = y_max;
cnt = 0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            cnt = cnt+1
            
print(cnt)
162/58:
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
np.savetxt('./state_'+str(state_index)+'_2.txt',state[::10,0:2],fmt='%16.8lf')
162/59:
samples_pre = state[::10,2:].astype('int32')
np.savetxt('./samples_0_2.txt',samples_pre,fmt='%6d')
samples = np.array([np.array(samples_pre)])
162/60: pyemma.coordinates.save_trajs(inp,samples,outfiles=['./state_{}_2.xtc'.format(n,n) for n in {state_index,state_index}])
162/61:
state_no = 1; state_index = np.int32(state_no);
x_min = 0.30; x_max = 0.310; y_min = 0.25; y_max = 0.26; #State1 final
X_min = x_min; X_max = x_max; Y_min = y_min; Y_max = y_max;
cnt = 0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            cnt = cnt+1
            
print(cnt)
162/62:
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
np.savetxt('./state_'+str(state_index)+'_2.txt',state[::10,0:2],fmt='%16.8lf')
162/63:
samples_pre = state[::100,2:].astype('int32')
np.savetxt('./samples_1_2.txt',samples_pre,fmt='%6d')
samples = np.array([np.array(samples_pre)])
162/64: pyemma.coordinates.save_trajs(inp,samples,outfiles=['./state_{}_2.xtc'.format(n,n) for n in {state_index,state_index}])
162/65:
state_no = 2; state_index = np.int32(state_no);
x_min = 0.50; x_max = 0.55; y_min = 0.09; y_max = 0.10; #State1 final
X_min = x_min; X_max = x_max; Y_min = y_min; Y_max = y_max;
cnt = 0
for i in range(np.shape(data)[0]):
    for j in range(np.shape(data[i])[0]):
        if( (data[i][j][0] > X_min) and (data[i][j][0] < X_max) and (data[i][j][1] > Y_min) and (data[i][j][1] < Y_max) ):
            cnt = cnt+1
            
print(cnt)
162/66:
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
np.savetxt('./state_'+str(state_index)+'_2.txt',state[::10,0:2],fmt='%16.8lf')
162/67:
samples_pre = state[::100,2:].astype('int32')
np.savetxt('./samples_2_2.txt',samples_pre,fmt='%6d')
samples = np.array([np.array(samples_pre)])
162/68: pyemma.coordinates.save_trajs(inp,samples,outfiles=['./state_{}_2.xtc'.format(n,n) for n in {state_index,state_index}])
   1: %history -g

In [2]: exit  
