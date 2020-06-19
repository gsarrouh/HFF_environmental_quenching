#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:53:55 2019

@author: gsarrouh
"""
import numpy as np
import matplotlib.pyplot as plt
#
QM_ml1 = 8.96
QM_ml2 = 11.18
Qphi_ml1 = 174.10
Qphi_ml2 = 18.92
Qalpha_ml1 = -0.65
Qalpha_ml2 = -0.7
#
x = np.linspace(SF_midbins[0],SF_midbins[len(SF_midbins)-1],num=200)
y1 = np.log(10)*Qphi_ml1*(10**((x-QM_ml1)*(1+Qalpha_ml1)))*np.exp(-10**(x-QM_ml1))
y2 = np.log(10)*Qphi_ml2*(10**((x-QM_ml2)*(1+Qalpha_ml2)))*np.exp(-10**(x-QM_ml2))
y = y1 + y2
#visualize
plt.close()
SMF = plt.figure(num=1)
plt.plot(x,y)
plt.xscale('linear')
plt.xlim=(8,12.25)
plt.yscale('log')
plt.ylim=(1,750)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off',labelbottom='on')
##################
##################
##################
##################
#
model = np.array([0.0]*len(SF_midbins))
model_sum = np.array([0.0]*len(SF_midbins))
smf = SF_smf
midbins = SF_midbins
#[M_star, phi, alpha] = [8.8264,8.38065,4.83094]
[M_star, phi, alpha] = [8.768857,10.27717,-0.72368]
for ii in range(len(midbins)):
    model[ii] = np.log(10)*phi*(10**((midbins[ii]-M_star)*(1+alpha)))*np.exp(-10**(midbins[ii]-M_star))
    model_sum[ii] = ((smf[ii] - model[ii])**2 / model[ii])