#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 17:48:27 2017

@author: gsarrouh
"""
###################### field_hist.py  #########################
#
## This program plots a histogram of z_phot and magnitude(F160W filter) for the parallel field samples
#
##  Fig.1: z_phot histogram
##  Fig.2: magnitude (F160W) histogram
#
##  Fig.1: z_phot histogram
plt.close()
zphot = Column(np.empty(shape=(30610,1), dtype='float16'))
z_index = 0
counter = 0
size = len(master_cat_par)
while counter < size:
    if master_cat_par[counter]['member'] == 1 and master_cat_par[counter]['sub'] ==2:
        zphot[z_index][0] = master_cat_par[counter]['z_peak']
        z_index +=1
    counter +=1
# put in bins
z_hist, z_bins = np.histogram(zphot, bins=25,range=(0,2))
z_midbins = np.empty_like(z_hist, dtype='float64')
#
size = len(z_bins)-1
counter = 0
while counter < size:
    z_midbins[counter] = (z_bins[counter] + z_bins[counter+1])/2
    counter +=1
# plot
#
fig = plt.figure(num=1)
ax = fig.add_subplot(111)
plt.scatter(z_midbins, z_hist, c='b',marker='.',linewidths=0)
plt.xlabel("$z_{phot}$")
plt.xscale('linear')
plt.xlim(0,2)
plt.ylabel("# count")
plt.yscale('linear')
#plt.ylim(0,1.25)
plt.title("$z_{phot} histogram$")
#plt.legend((memz,outz),('spec population','outliers'),scatterpoints=1,loc='upper left', frameon=False)
#plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on')
plt.minorticks_on()
plt.show()
#
## ##  Fig.2: magnitude (F160W) histogram
#
plt.close()
zees = plt.figure(num=2)
zvz = zees.add_subplot(111)
##plot spec subsample using master_cat['sub']==1 as identifier for loop
#plotzs = np.array
z_ph = []
z_cl = []
other = 0
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['member'] == 1 and master_cat[counter]['sub'] ==2:
        z_ph.append(master_cat[counter]['z_peak'])
        z_cl.append(master_cat[counter]['z_clusterphot'])
    else:
        other +=1
    counter+=1
plt.scatter(z_ph,z_cl,c='k', marker='.', linewidths = 0)
#plt.plot(z_ph,z_cl,'.k',linewidth=0.0)
plt.xlabel("$z_{phot}$")
plt.xscale('linear')
#plt.xlim(0,1.25)
plt.ylabel("$z_{phot} - z_{cl} / (1 + z_{cl})$")
plt.yscale('linear')
#plt.ylim(0,1.25)
plt.title("$z_{phot} vs. z_{clusterphot}$")
#plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on')
plt.minorticks_on()
plt.show()    