import pyemma
import os
import matplotlib.pyplot as plt
import numpy as np
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from pyemma import config

dtrajs = np.load('dtrajs.npy',allow_pickle=True)
dtrajs = list(dtrajs)

lags=1000
its = pyemma.msm.its(dtrajs, lags=lags , nits=10)

# Plotting its plots

plt.figure(figsize=(10,8))
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2)
plt.ylim(0, )
plt.xlim(0, 10)
plt.grid(axis='both')
plt.savefig('./its_1.png')
plt.clf()

plt.figure(figsize=(10,8))
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2)
plt.ylim(0, )
plt.xlim(0, 25)
plt.grid(axis='both')
plt.savefig('./its_2.png')
plt.clf()


plt.figure(figsize=(10,8))
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2)
plt.ylim(0, )
plt.xlim(0, 50)
plt.grid(axis='both')
plt.savefig('./its_3.png')
plt.clf()

plt.figure(figsize=(10,8))
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2)
plt.ylim(0, )
plt.xlim(0, 100)
plt.grid(axis='both')
plt.savefig('./its_4.png')
plt.clf()

plt.figure(figsize=(10,8))
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2)
plt.ylim(0, )
plt.xlim(0, 200)
plt.grid(axis='both')
plt.savefig('./its_5.png')
plt.clf()

plt.figure(figsize=(10,8))
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2)
plt.ylim(0, )
plt.xlim(0, 500)
plt.grid(axis='both')
plt.savefig('./its_6.png')
plt.clf()

plt.figure(figsize=(10,8))
mplt.plot_implied_timescales(its, ylog=True, show_mean=False, units='10ps', linewidth=2)
plt.ylim(0, )
plt.xlim(0, 1000)
plt.grid(axis='both')
plt.savefig('./its_7.png')
plt.clf()
