#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 21:32:30 2017

@author: gsarrouh
"""
#####################  data_completeness  #####################
#
###  This program creates plots related to completeness of the data.
#
##  Fig.1: apparent magnitude histogram (F160 filter, i.e. 1600nm), 'F160_data' for data
##  Fig.2: 
#
# 
## Fig 1: plot histogram of F160 magnitude (show spec vs phot data)
#
# create list of flux values at 1600nm
# note: plotting all objects for which there are data (i.e. excluding only 'nodata')
smag = []
serr = []
pmag = []
perr = []
other = 0
out1 = 0
out2 = 0
out3 = 0
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['sub'] == 0:   #skip no data
        out1 +=1
    elif master_cat[counter]['sub'] ==1 or master_cat[counter]['sub'] ==3:      #spec & spec_only
        if master_cat[counter]['f_F160W'] >0:
            smag.append(master_cat[counter]['f_F160W'])
            serr.append(master_cat[counter]['e_F160W'])
        else: out2 +=1                  # # of spec objects w/o F160 flux
    elif master_cat[counter]['sub'] ==2:
        if master_cat[counter]['f_F160W'] >0:
            pmag.append(master_cat[counter]['f_F160W'])
            perr.append(master_cat[counter]['e_F160W'])
        else: out3 +=1                  # # of spec objects w/o F160 flux
    else:
        other +=1
    counter+=1
#
# a table of what i'm plotting
names = Column(['Total cat','spec','phot','no data objects','spec w/o F160','phot w/o F160'])
values = Column([size,len(smag),len(pmag),out1,out2,out3])
F160_data = Table([names,values],names=('Property','Value'))
# convert flux to magnitude
f_s = np.empty_like(smag, dtype='float16')
f_p = np.empty_like(pmag, dtype='float16')
counter = 0
size = len(smag)        #spec
while counter < size:
    f_s[counter] = smag[counter]
    smag[counter] = -2.5*np.log10(smag[counter])+25
    counter+=1
counter = 0
size = len(pmag)        #phot
while counter < size:
    f_p[counter] = pmag[counter]
    pmag[counter] = -2.5*np.log10(pmag[counter])+25
    counter+=1
# put in bins
bins = [0,15,17,19,20,21,22,23,23.5,24,24.5,25,25.5,26,26.5,27,27.5,28,29,30]
bins_fine = [0,18,18.5,19,19.25,19.5,19.75,20,20.25,20.5,20.75,21,21.25,21.5,21.75,22,22.25,22.5,22.75,23,24,25,25.25,25.5,25.75,26,26.25,26.5,26.75,27,27.25,27.5,27.75,28,28.5,29,30]
smag_hist, mag_bins = np.histogram(smag, bins=bins,range=(15,30))
pmag_hist, mag_bins = np.histogram(pmag, bins=mag_bins,range=(15,30))
mag_midbins = np.empty_like(smag_hist, dtype='float64')
#smaller bins around peak
smag_hist_fine, mag_bins_fine = np.histogram(smag, bins=bins_fine,range=(15,30))
pmag_hist_fine, mag_bins_fine = np.histogram(pmag, bins=mag_bins_fine,range=(15,30))
mag_midbins_fine = np.empty_like(smag_hist_fine, dtype='float64')
#midbins for both sets
size = len(mag_bins)-1
counter = 0
while counter < size:
    mag_midbins[counter] = (mag_bins[counter] + mag_bins[counter+1])/2
    counter +=1
# finer midbins around peak
size = len(mag_bins_fine)-1
counter = 0
while counter < size:
    mag_midbins_fine[counter] = (mag_bins_fine[counter] + mag_bins_fine[counter+1])/2
    counter +=1
## Plot
# left figure: full set
plt.close()
F160 = plt.figure(num=1)
flux = F160.add_subplot(121)
spec = flux.plot(mag_midbins,smag_hist,'.g',linewidth=0.5)
phot = flux.plot(mag_midbins,pmag_hist,'.m',linewidth=0.5)
plt.xlabel("$m_{F160W}$")
plt.xscale('linear')
#plt.gca().invert_xaxis()
plt.xlim(15,31)
plt.ylabel("count")
plt.yscale('log')
#plt.ylim(5,15)
plt.title('Apparent Magnitude Histogram')
#plt.legend((spec,phot),('spec','phot'),scatterpoints=1,loc='best', frameon=False)
#plt.legend((spec,phot),('spec','phot'),scatterpoints=1,loc='best', frameon=False)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off')
#
# right figure: finer bins
flux = F160.add_subplot(122)
spec_fine = flux.plot(mag_midbins_fine,smag_hist_fine,'.g',linewidth=0.5)
phot_fine = flux.plot(mag_midbins_fine,pmag_hist_fine,'.m',linewidth=0.5)
plt.xlabel("$m_{F160W}$")
plt.xscale('linear')
#plt.gca().invert_xaxis()
plt.xlim(15,31)
plt.ylabel("count")
flux.yaxis.set_label_position("right")
plt.yscale('log')
#plt.ylim(5,15)
plt.title('Hist. w/ finer bins')
#flux.legend((smag_hist_fine.any(),pmag_hist_fine.any()),('spec','phot'),scatterpoints=1,loc='upper left', frameon=False)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on',labelleft='off')
#plt.tick_params(axis='y', which='both', direction='in',color='k',top='on',right='on',labelright='on',labelleft='off')
plt.show()  
#
plt.close()#
## Fig.2: (flux/error) vs. magnitude, all in F160W filter
#lines of best fit for both plots
err_s = f_s/serr    #calculate flux / delta_fluc
err_p = f_p/perr
#spec bestfit
z_spec = np.polyfit(smag, err_s, 1)   # calculate polynomial
f_spec = np.poly1d(z_spec)
smag_new = np.linspace(min(smag), max(smag), 100)  # calculate new x's and y's
err_s_new = f_spec(smag_new)
#phot bestfit
z_phot = np.polyfit(pmag, err_p, 1)   # calculate polynomial
f_phot = np.poly1d(z_phot)
pmag_new = np.linspace(min(pmag), max(pmag), 100)  # calculate new x's and y's
err_p_new = f_phot(pmag_new)
#create figure & plot
e_F160 = plt.figure(num=2)
plt.suptitle('Flux/Error vs. Apparent Magnitude')
spec_error = e_F160.add_subplot(121)
e_spec = spec_error.plot(smag,err_s,'.g',linewidth=0.5)
plt.plot([0,40],[5,5],':r', linewidth=1)
plt.plot(smag_new,err_s_new,'-k',linewidth=1)
plt.xlabel("$m_{F160W}$")
plt.xscale('linear')
#plt.gca().invert_xaxis()
plt.xlim(15,40)
plt.ylabel("$flux_{F160W}/error$")
plt.yscale('log')
#plt.ylim(5,15)
plt.title('Spec')
#plt.legend((e_spec,e_phot),('spec','phot'),scatterpoints=1,loc='best', frameon=False)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off')
#
phot_error = e_F160.add_subplot(122)
e_phot = phot_error.plot(pmag,err_p,'.m',linewidth=0.5)
plt.plot([0,40],[5,5],':r', linewidth=0.8)
plt.plot(pmag_new,err_p_new,'-k',linewidth=1)
plt.xlabel("$m_{F160W}$")
plt.xscale('linear')
#plt.gca().invert_xaxis()
plt.xlim(15,40)
plt.ylabel("$flux_{F160W}/error$")
phot_error.yaxis.set_label_position("right")
plt.yscale('log')
#plt.ylim(5,15)
plt.title('Phot')
#plt.legend((e_spec,e_phot),('spec','phot'),scatterpoints=1,loc='best', frameon=False)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on', labelleft='off')

plt.show()  