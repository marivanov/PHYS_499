#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 09:58:51 2017

@author: mivanov
"""

from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
import corner
from bayes_linear import bayes_linear
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

#%% ### Initial math ###

#for i in range(len(lxmin)):
#    if lxmin[i]==0:
#        lxmin[i]=0.9*lxmax

lxmean = (lxmax+lxmin)/2
lx_m = lxmean/mass
lx_m[32] = 9e27
lx_m_err = ((lxmax-lxmean)*lx_m)/lxmean
lx_m_err[32] = 3e27

#%% ### Changing powers ###

#n = 1.7
#agetop = (age+age_err)**n
#agebot = (age-age_err)**n
#age = (agetop+agebot)/2
#age_err = agetop - age
#
#m = 0.3
#pos=[]
#neg=[]
#for i in range(len(metal)):
#    if metal[i]>0:
#        pos.append(i)
#    else:
#        neg.append(i)
#
#metaltop = np.zeros(len(metal))
#metalbot = np.zeros(len(metal))
#metaltop[pos] = (metal[pos]+metal_err[pos])**m
#metalbot[pos] = (metal[pos]-metal_err[pos])**m
#metaltop[neg] = -1*(np.abs((metal[neg]+metal_err[neg]))**m)
#metalbot[neg] = -1*(np.abs((metal[neg]-metal_err[neg]))**m)
#metal = (metaltop+metalbot)/2
#metal_err = metaltop-metal
#
#l = 0.3
#binfractop = (binfrac+binfracerr)**l
#binfracbot = (binfrac-binfracerr)**l
#binfrac = (binfractop+binfracbot)/2
#binfracerr = binfractop - binfrac

#%% ### Defining variables used to plot and to work with ###

# In this block, follow these definitions:
#           LxM = Luminosity per unit Mass in log10
#           LxMerr = Error in Luminosity per unit Mass in log10
#           A_M = Age times Metallicity
#           A_Merr = Error in Age times Metallicity

A_Mindex=[]
for i in range(len(ro)):
    if ro[i]<10000 and age[i]!=0 and lx_m_err[i]!=0 and metal[i]!=0 and binfrac[i]!=0:
        A_Mindex.append(i)
        
datatype = data_type[A_Mindex]
GCunder6=[]
GCover6=[]
notGC=[]
for i in range(len(datatype)):
    if datatype[i] == 1:
        GCunder6.append(i)
    if datatype[i] == 2:
        GCover6.append(i)
    if datatype[i] == 3:
        notGC.append(i)

LxM = np.log10(lx_m[A_Mindex])
err = lx_m[A_Mindex]+lx_m_err[A_Mindex]
LxMerr = np.log10(err)-LxM

A_M = metal[A_Mindex]*binfrac[A_Mindex]/age[A_Mindex]
alpha = age[A_Mindex]
alpha_err = age_err[A_Mindex]
beta = metal[A_Mindex]
beta_err = metal_err[A_Mindex]
gamma = binfrac[A_Mindex]
gamma_err = binfracerr[A_Mindex]
A_Merr = A_M*np.sqrt((alpha_err/alpha)**2 + 
                     (beta_err/beta)**2 +
                     (gamma_err/gamma)**2)

#%% ### Starting the emcee part ###

# Calling in the data
xerr, yerr = A_Merr, LxMerr
x_obs, y_obs = A_M, LxM

N = len(A_M)
estimates = np.array([60,22])

params, error_int, samples = bayes_linear(x_obs,y_obs,xerr,yerr,p1=estimates,nWalkers=200,nBurn=200,nSample=4000)

m = params[0]
b = params[1]
m_err = error_int[0,:]
b_err = error_int[1,:]

x = np.linspace(min(x_obs+xerr),max(x_obs+xerr),100)
yn = m*x+b

# Plot results
fig1 = plt.figure(1,figsize=(10, 5))

ax1 = fig1.add_subplot(121)
ax1.hist(samples[0,:], 50, histtype="step", color="k")
ax1.axvline(estimates[0])
ax1.set_yticklabels([])
ax1.set_xlabel("$m$")

ax2 = fig1.add_subplot(122)
ax2.hist(samples[1,:], 50, histtype="step", color="k")
ax2.axvline(estimates[1])
ax2.set_yticklabels([])
ax2.set_xlabel("$b$")

fig2 = plt.figure(2)
plt.clf()
#frame1 = fig2.add_axes((0.1,0.3,0.8,0.6))

mMax = m_err[1]
mMin = m_err[0]
bMax = b_err[1]
bMin = b_err[0]
yMax = mMax*x + bMin
yMin = mMin*x + bMax

plt.errorbar(x_obs,y_obs,yerr=yerr,xerr=xerr,ecolor='k',ls='none',
             elinewidth=0.5,capsize=2,markeredgewidth=0.5)
plt.plot(x,yn,'k-')
plt.plot(x,yMax,ls='none')
plt.plot(x,yMin,ls='none')
plt.fill_between(x,yMin,yMax,color='blue',alpha=0.2)
plt.plot(x_obs[GCunder6],y_obs[GCunder6],'or',ms=5,label='GC < 6kpc')
plt.plot(x_obs[GCover6],y_obs[GCover6],'sg',ms=5,label='GC > 6kpc')
plt.plot(x_obs[notGC],y_obs[notGC],'dy',ms=5,label='Open Clusters')
plt.title('X-ray Luminosity per Unit Mass vs (Binary Fraction x Metallicity)/Age')
#plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
plt.ylabel(r'X-ray Luminosity per Unit Mass, $log_{10}(L_{X}/M_{\bigodot})$')
plt.axis([min(x_obs+xerr),max(x_obs+xerr), 24, 29])
plt.legend(loc='best')
plt.grid()

# Residual Calculation
#Res = ((x_obs*m) + b) - y_obs
#StDev = Res/yerr
#
#frame2 = fig2.add_axes((0.1,0.1,0.8,0.2))
#plt.errorbar(x_obs,StDev,yerr=yerr,xerr=xerr,ecolor='k',ls='none',
#             elinewidth=0.5,capsize=2,markeredgewidth=0.5,label='GC and Non-GC')
#plt.plot(x_obs[GCunder6],StDev[GCunder6],'or',ms=5)
#plt.plot(x_obs[GCover6],StDev[GCover6],'sg',ms=5)
#plt.plot(x_obs[notGC],StDev[notGC],'dy',ms=5)
#xn = np.linspace(min(x_obs+xerr),max(x_obs+xerr),100)
#yn = np.zeros(100)
#plt.plot(xn,yn,'-k')
plt.xlabel(r'(Binary Fraction $\times$ Metallicity)/Age, $(\mathrm{binfrac}\times[Fe/H])/Gyr$')
#plt.axis([min(x_obs+xerr),max(x_obs+xerr),-15, 15])
#plt.grid()

fig3 = plt.figure(3)

plt.subplot(2,1,1)
plt.plot(samples[0,:])
plt.ylabel('m')

plt.subplot(2,1,2)
plt.plot(samples[1,:])
plt.ylabel('b')

fig4 = corner.corner(samples.T, labels=["$m$", "$b$"])

print('')
print(' EMCEE ANALYSIS ')
print('')
print('The Slope is', ("%.4f" % m))
print('The Slope Interval is {0:0.4f} to {1:0.4f}'.format(m_err[0], m_err[1]))
print('The Intercept is', ("%.4f" % b))
print('The Intercept Interval is {0:0.4f} to {1:0.4f}'.format(b_err[0], b_err[1]))

#%% ### Calculating the Reduced Chi Squared value

ChiSq = np.sum((((m*x_obs+b)-y_obs)**2)/
                     (yerr**2+(xerr*m)**2))
RedChiSq = ChiSq / (len(x_obs)-2)

print('')
print('The degrees of freedom is', ("%.4f" % (len(x_obs)-2)))
print('The chi squared value is', ("%.4f" % ChiSq))
print('The reduced chi squared value is', ("%.4f" % RedChiSq))
print('')