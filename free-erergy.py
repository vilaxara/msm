# coding: utf-8
import numpy as np
import scipy

'''

Computing free-energies and rates from stationary populations

'''
c = 10.2*10**(-3)  #10.2mM/L
RT = 0.6163314 # RT in kcal/mol at 310.15K
dGv = -(RT*np.log(473789/1661.0)) # Volume correction of energy in kcal/mol

#####Free-energy from populations#####

S_ub = input('Stationary population of UB : ')
S_ub = float(S_ub)

S_b = input('Stationary population of B : ')
S_b = float(S_b)

dG = -(RT*np.log(S_b/S_ub))

dG0 = dG + dGv

print ('Result] \u0394G with populations : ', np.round(dG0,3),' kcal/mol')

#####Binding Rates#####

M_on = input('MFPT ON (\u03BCs) :')
M_on = float(M_on)*(10**(-6))

M_off = input('MFPT OFF (\u03BCs) :')
M_off = float(M_off)*(10**(-6))

k_on = 1.0/(M_on*c)
k_off = 1.0/(M_off)

kd = k_off/k_on

dG_kd = RT*np.log(kd)

print ('Result] k_on : ', "{:.3e}".format(k_on))
print ('Result] k_off : ', "{:.3e}".format(k_off))

print ('Result] \u0394G with kd_simu :', np.round(dG_kd,3), ' kcal/mol')
