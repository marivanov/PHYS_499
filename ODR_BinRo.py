#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 23:38:09 2018

@author: mivanov
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.odr import ODR, Model, RealData
np.set_printoptions(threshold=np.nan)


### Importing Data ###
data_full = np.genfromtxt('heinke_final_data.csv', delimiter=',', 
                          skip_header=1, usecols=(3,4,6,7,8,9,10,
                                                  11,13,14,15,17))

age = data_full[:,0]
age_err = data_full[:,1]
metal = data_full[:,2]
metal_err = data_full[:,3]
data_type = data_full[:,4]
binfrac = data_full[:,5]
binfracerr = data_full[:,6]
binfractype = data_full[:,7]
ro = data_full[:,8]
mass = data_full[:,9]
lxmax = data_full[:,10]
lxmin = data_full[:,11]

lxmean = (lxmax+lxmin)/2
lx_m = lxmean/mass
lx_m[32] = 9e27
lx_m_err = ((lxmax-lxmean)*lx_m)/lxmean
lx_m_err[32] = 3e27


#%% Cleaning up the data a bit (keeping density under 10000 and no null data points for binary fraction)

def funcLin(p,x):
    return p[0]+p[1]*x

index=[]
for i in range(len(ro)):
    if ro[i]<10000 and binfrac[i]!=0:
        index.append(i)

Bin_I = binfrac[index]
BinErr_I = binfracerr[index]
Ro_I = ro[index]

### And now to calculate error in log-space
Ro_Log = np.log10(Ro_I)
Bin_Log = np.log10(Bin_I)
BinErr_Log = np.log10(Bin_I+BinErr_I) - np.log10(Bin_I)
xn = np.linspace(min(Bin_Log-BinErr_Log),max(Bin_Log+BinErr_Log),100)

#%% Setting up the ODR process, the plot arrays and chi squared ###

data = RealData(Bin_Log[0:14],Ro_Log[0:14],sx=BinErr_Log[0:14])
modelLM = Model(funcLin)
odrLM = ODR(data,modelLM,[0.70,-1.95])
odrLM.set_job(fit_type=0)
outLM = odrLM.run()
yn = funcLin(outLM.beta,xn)

#%% Plotting binary fraction vs density cleaned up

plt.figure(2)
plt.clf()
plt.plot(Bin_Log,Ro_Log,'k.')
plt.errorbar(Bin_Log,Ro_Log,xerr=BinErr_Log, ecolor='k',ls='none',elinewidth=0.5,capsize=2,markeredgewidth=0.5)
plt.plot(xn,yn,'k-',linewidth=0.5,label='Linear ODR Fit')
plt.title(r'Log-Log space of Binary fraction vs. Density')
plt.xlabel('Binary fraction, $log_{10}(\mathrm{binfrac})$')
plt.ylabel(r'Density, $log_{10}(\rho)$')
plt.grid()


print('Following y = mx + b for the Linear Fit of Lx/M vs. Metallicity')
print('     m =', ("%.2f" % outLM.beta[1]))
print('     b =', ("%.2f" % outLM.beta[0]))
