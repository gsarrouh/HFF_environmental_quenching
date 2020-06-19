#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 12:18:25 2017

@author: gsarrouh
"""

###  This program plots mass vs apparent magnitude (J restframe "161" from HFF data file; find out what wavelength this corresponds to) for the HFF catalogue
#
## 
## PLOT 1: plot lmass vs m_j (show spec vs phot data)
#
mcol = plt.figure(num=1)
spph = mcol.add_subplot(121)
spec = []
smag = []
phot = []
pmag = []
other = 0
out = 0
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['type'] == 0 or master_cat[counter]['type'] == 3 or master_cat[counter]['sub'] == 0:   #skip stars, outliers, nodata
        out +=1
    elif master_cat[counter]['sub'] == 1 or master_cat[counter]['sub'] == 3:
        spec.append(master_cat[counter]['lmass'])
        smag.append(master_cat[counter]['L_j'])
    elif master_cat[counter]['sub'] ==2:
        phot.append(master_cat[counter]['lmass'])
        pmag.append(master_cat[counter]['L_j'])
    else:
        other +=1
    counter+=1
plt.plot(pmag,phot,'.m',linewidth=0.5)
plt.plot(smag,spec,'.g',linewidth=0.5)
plt.xlabel("$m_{j}$")
plt.xscale('linear')
#plt.gca().invert_xaxis()
plt.xlim(15,40)
plt.ylabel("log(M/M_{\odot})$)")
plt.yscale('linear')
plt.ylim(5,15)
plt.title('mass vs. apparent magnitude')
plt.legend((smag,pmag),('spec','phot'),scatterpoints=1,loc='upper left', frameon=False)
#plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off')
plt.minorticks_on()
plt.show()    
#
# right panel: show SF & Q populations
#
sfq = mcol.add_subplot(122)
sf = []
sfmag = []
q = []
qmag = []
other = 0
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['type'] == 0 or master_cat[counter]['type'] == 3 or master_cat[counter]['sub'] == 0:   #skip stars, outliers, nodata
        out +=1
    elif master_cat[counter]['type'] == 1:
        sf.append(master_cat[counter]['lmass'])
        sfmag.append(master_cat[counter]['L_j'])
    elif master_cat[counter]['type'] ==2:
        q.append(master_cat[counter]['lmass'])
        qmag.append(master_cat[counter]['L_j'])
    else:
        other +=1
    counter+=1
plt.plot(sfmag,sf,'.b',linewidth=0.5)
plt.plot(qmag,q,'.r',linewidth=0.5)
plt.xlabel("$m_{j}$")
plt.xscale('linear')
plt.xlim(15,40)
plt.ylabel("$log(M/M{sol})$")
plt.yscale('linear')
plt.ylim(5,15)
#plt.title('mass vs. apparent magnitude')
plt.legend((sf,q),('SF','Q'),scatterpoints=1,loc='upper left', frameon=False)
#plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',left = 'on', right='on',labelleft='off',labelright='on')
sfq.yaxis.set_label_position("right")
plt.minorticks_on()
plt.show()