#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 10:04:02 2017

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


#%% ### Initial math ###

lxmean = (lxmax+lxmin)/2
lx_m = lxmean/mass
#lx_m[32] = 9e27 # Local LxM
lx_m_err = ((lxmax-lxmean)*lx_m)/lxmean
#lx_m_err[32] = 3e27 # Local LxM error


### Defining the Linear model and selecting the data we want ###
def funcLin(p,x):
    return p[0]+p[1]*x

index=[]
for i in range(len(ro)):
    if ro[i]<10000 and binfrac[i]!=0 and lx_m_err[i]!=0:
        index.append(i)

datatype = data_type[index]
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

### Defining variables used to plot and to work with ###
B = binfrac[index]
Berr = binfracerr[index]
LxM = np.log10(lx_m[index])
err = lx_m[index]+lx_m_err[index]
LxMerr = np.log10(err)-LxM
xn = np.linspace(min(B),max(B),100)


### Setting up the ODR process, the plot arrays and chi squared ###
data = RealData(B,LxM,sx=Berr,sy=LxMerr)
modelLM = Model(funcLin)
odrLM = ODR(data,modelLM,[22,2])
odrLM.set_job(fit_type=0)
outLM = odrLM.run()
ynLinearMetal = funcLin(outLM.beta,xn)
Res = funcLin(outLM.beta,B) - LxM
StDev = Res/LxMerr
ChiSq = np.sum(((funcLin(outLM.beta,B)-LxM)/LxMerr)**2)
RedChiSq = ChiSq / (len(B)-len(outLM.beta))


### PLotting ###
fig7 = plt.figure(1)
plt.clf()
#frame1 = fig7.add_axes((0.1,0.3,0.8,0.6))
plt.errorbar(B,LxM,yerr=LxMerr,xerr=Berr,
             ecolor='k',ls='none',elinewidth=0.5,capsize=2,markeredgewidth=0.5)
plt.plot(xn,ynLinearMetal,'k-')
plt.plot(B[GCunder6],LxM[GCunder6],'or',ms=5,label='GC < 6kpc')
plt.plot(B[GCover6],LxM[GCover6],'sg',ms=5,label='GC > 6kpc')
plt.plot(B[notGC],LxM[notGC],'dy',ms=5,label='Non-GC')
plt.title('X-ray Luminosity per Unit Mass vs Binary Fraction')
plt.xlabel(r'Binary Fraction')
plt.axis([min(B-Berr), 0.7, 25, 29])
#plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
plt.ylabel(r'X-ray Luminosity per Unit Mass, $log_{10}(L_{X}/M_{\bigodot})$')
plt.legend(loc='best')
plt.grid()

#frame2 = fig7.add_axes((0.1,0.1,0.8,0.2))
#plt.plot(M,StDev,'ok',ms=3.5)
#yblackline = np.zeros(100)
#plt.plot(xn,yblackline,'-b')
#plt.xlabel(r'Metallicity, $[Fe/H]$')
#plt.axis([-2.5, 0.5, -3, 3])
#plt.grid()

print(' ')
print('The degrees of freedom for the Linear Metallicity test is', len(B) - len(outLM.beta))
print('The Linear Metallicity fit chi squared value is', ("%.2f" % ChiSq))
print('The Linear Metallicity fit reduced chi squared value is', ("%.2f" % RedChiSq))
print(' ')
print('Following y = mx + b for the Linear Fit of Lx/M vs. Metallicity')
print('     m =', ("%.2f" % outLM.beta[1]))
print('     b =', ("%.2f" % outLM.beta[0]))





#%% LIMITING BINARY FRACTION TO UNDER 0.15 OR ALMOST CONSTANT

BinIndex=[]
for i in range(len(ro)):
    if ro[i]<10000 and binfrac[i]!=0 and binfrac[i]<0.15 and lx_m_err[i]!=0:
        BinIndex.append(i)

B = binfrac[BinIndex]
Berr = binfracerr[BinIndex]
LxM = np.log10(lx_m[BinIndex])
err = lx_m[BinIndex]+lx_m_err[BinIndex]
LxMerr = np.log10(err)-LxM
xn = np.linspace(min(B-Berr),max(B+Berr))

datatype = data_type[BinIndex]
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

#### Setting up the ODR process, the plot arrays and chi squared ###
#data = RealData(B,LxM,sx=Berr,sy=LxMerr)
#modelLM = Model(funcLin)
#odrLM = ODR(data,modelLM,[22,2])
#odrLM.set_job(fit_type=0)
#outLM = odrLM.run()
#ynLinearMetal = funcLin(outLM.beta,xn)
#Res = funcLin(outLM.beta,B) - LxM
#StDev = Res/LxMerr
#ChiSq = np.sum(((funcLin(outLM.beta,B)-LxM)/LxMerr)**2)
#RedChiSq = ChiSq / (len(B)-len(outLM.beta))
#
fig7 = plt.figure(2)
plt.clf()
#frame1 = fig7.add_axes((0.1,0.3,0.8,0.6))
plt.errorbar(B,LxM,yerr=LxMerr,xerr=Berr,
             ecolor='k',ls='none',elinewidth=0.5,capsize=2,markeredgewidth=0.5)
#plt.plot(xn,ynLinearMetal,'k-',label='Linear ODR Fit')
plt.plot(B[GCunder6],LxM[GCunder6],'or',ms=5,label='GC < 6kpc')
plt.plot(B[GCover6],LxM[GCover6],'sg',ms=5,label='GC > 6kpc')
plt.plot(B[notGC],LxM[notGC],'dy',ms=5,label='Non-GC')
plt.title('X-ray Luminosity per Unit Mass vs Binary Fraction')
plt.xlabel(r'Binary Fraction')
plt.axis([min(B-Berr), 0.7, 25, 29])
#plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
plt.ylabel(r'X-ray Luminosity per Unit Mass, $log_{10}(L_{X}/M_{\bigodot})$')
plt.legend(loc='best')
plt.grid()






#%% AGE VS LX WITH RESPECT TO BINARY FRACTION BEING ALMOST CONSTANT

A = age[BinIndex]
Aerr = age_err[BinIndex]

#data = RealData(A,LxM,sx=Aerr,sy=LxMerr)
#xn = np.linspace(2,13.5,100)
#modelLA = Model(funcLin)
#odrLA = ODR(data,modelLA,[28,-2])
#odrLA.set_job(fit_type=0)
#outLA = odrLA.run()
#ynLA = funcLin(outLA.beta,xn)
#Res = funcLin(outLA.beta,A) - LxM
#StDev = Res/LxMerr
#ChiSq = np.sum(((funcLin(outLA.beta,A)-LxM)/LxMerr)**2)
#RedChiSq = ChiSq / (len(A)-len(outLA.beta))

fig6 = plt.figure(3)
plt.clf()
#frame1 = fig6.add_axes((0.1,0.3,0.8,0.6))
plt.errorbar(A,LxM,yerr=LxMerr,xerr=Aerr,
             ecolor='k',ls='none',elinewidth=0.5,capsize=2,markeredgewidth=0.5)
#plt.plot(xn,ynLA,'k-',label='Linear ODR Fit')
plt.plot(A[GCunder6],LxM[GCunder6],'or',ms=5,label='GC < 6kpc')
plt.plot(A[GCover6],LxM[GCover6],'sg',ms=5,label='GC > 6kpc')
plt.plot(A[notGC],LxM[notGC],'dy',ms=5,label='Non-GC')
plt.title('X-ray Luminosity per Unit Mass vs Age')
plt.xlabel(r'Age, Gyr')
plt.axis([2, 13.5, 25, 30.5])
#plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
plt.ylabel(r'X-ray Luminosity per Unit Mass, $log_{10}(L_{X}/M_{\bigodot})$')
plt.legend(loc='best')
plt.grid()





#%% METAL VS LX WITH RESPECT TO BINARY FRACTION BEING ALMOST CONSTANT

M = metal[BinIndex]
Merr = metal_err[BinIndex]
xn = np.linspace(-2.5,0,100)

### Setting up the ODR process, the plot arrays and chi squared ###
#data = RealData(M,LxM,sx=Merr,sy=LxMerr)
#modelLM = Model(funcLin)
#odrLM = ODR(data,modelLM,[28,-2])
#odrLM.set_job(fit_type=0)
#outLM = odrLM.run()
#ynLinearMetal = funcLin(outLM.beta,xn)
#Res = funcLin(outLM.beta,M) - LxM
#StDev = Res/LxMerr
#ChiSq = np.sum(((funcLin(outLM.beta,M)-LxM)/LxMerr)**2)
#RedChiSq = ChiSq / (len(M)-len(outLM.beta))

fig7 = plt.figure(4)
plt.clf()
#frame1 = fig7.add_axes((0.1,0.3,0.8,0.6))
plt.errorbar(M,LxM,yerr=LxMerr,xerr=Merr,
             ecolor='k',ls='none',elinewidth=0.5,capsize=2,markeredgewidth=0.5)
#plt.plot(xn,ynLinearMetal,'k-',label='Linear ODR Fit')
plt.plot(M[GCunder6],LxM[GCunder6],'or',ms=5,label='GC < 6kpc')
plt.plot(M[GCover6],LxM[GCover6],'sg',ms=5,label='GC > 6kpc')
plt.plot(M[notGC],LxM[notGC],'dy',ms=5,label='Open clusters')
plt.title('X-ray Luminosity per Unit Mass vs Metallicity')
plt.xlabel(r'Metallicity, $[Fe/H]$')
plt.axis([-2.5, 0, 25, 29])
#plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
plt.ylabel(r'X-ray Luminosity per Unit Mass, $log_{10}(L_{X}/M_{\bigodot})$')
plt.legend(loc='best')
plt.grid()





#%% METAL/AGE VS LX WITH RESPECT TO BINARY FRACTION BEING ALMOST CONSTANT

A_M = metal[BinIndex]/age[BinIndex]
alpha = age[BinIndex]
alpha_err = age_err[BinIndex]
beta = metal[BinIndex]
beta_err = metal_err[BinIndex]
A_Merr = A_M*np.sqrt((alpha_err/alpha)**2 + 
                     (beta_err/beta)**2)

fig5 = plt.figure(5)
plt.clf()
#frame1 = fig7.add_axes((0.1,0.3,0.8,0.6))
plt.errorbar(A_M,LxM,yerr=LxMerr,xerr=A_Merr,
             ecolor='k',ls='none',elinewidth=0.5,capsize=2,markeredgewidth=0.5)
#plt.plot(xn,ynLinearMetal,'k-',label='Linear ODR Fit')
plt.plot(A_M[GCunder6],LxM[GCunder6],'or',ms=5,label='GC < 6kpc')
plt.plot(A_M[GCover6],LxM[GCover6],'sg',ms=5,label='GC > 6kpc')
plt.plot(A_M[notGC],LxM[notGC],'dy',ms=5,label='Open clusters')
plt.title('X-ray Luminosity per Unit Mass vs Metallicity')
plt.xlabel(r'Metallicity/Age, $[Fe/H]/Gyr$')
#plt.axis([-0.4, 0, 25, 28])
#plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
plt.ylabel(r'X-ray Luminosity per Unit Mass, $log_{10}(L_{X}/M_{\bigodot})$')
plt.legend(loc='best')
plt.grid()
