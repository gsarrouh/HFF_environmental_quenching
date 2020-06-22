#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 21:32:30 2017

@author: gsarrouh
"""
#####################  data_completeness_4  #####################
#
###  This program creates plots related to MASS completeness of the data.
#
##  v3: calcuates lines of best fit using scipy.optimize.curve_fit for non-linear regression
##  v4: comments out and removes plots for spec population, includes all spec 
##      w/ phot data in phot plots; calculates limiting mag/mass by cluster; 
##      includes all spec & phot objects, not just selected members to better represent 
##      entire data catalogue
#
##  Fig.1(A/B): apparent magnitude histogram (F160 filter, i.e. 1600nm), 'F160_data' for data; (A) auto binning; (B) fine binning
##  Fig.2: (flux/error) vs. magnitude, all in F160W filter; determine limiting mag: (27-0.7)= 26.3
##  Fig.3: mass-to-luminosity plot; determine limiting mass: 
#
# 
from scipy.optimize import curve_fit
from astropy.stats import median_absolute_deviation
## Fig 1: plot histogram of F160 magnitude (show spec vs phot data)
#
# create list of flux values at 1600nm
# note: plotting all objects for which there are data (i.e. excluding only 'nodata'), by cluster
pmag1 = []   #SF magnitude
perr1 = []   #SF error
pmass1 = []  #SF mass
prad1 = []   #SF radius
pmag2 = []   
perr2 = []
pmass2 = []
prad2 = []
pmag3 = []   
perr3 = []
pmass3 = []
prad3 = []
pmag4 = []   
perr4 = []
pmass4 = []
prad4 = []
pmag5 = []   
perr5 = []
pmass5 = []
prad5 = []
pmag6 = []   
perr6 = []
pmass6 = []
prad6 = []
pmag1_mem = []
pmag2_mem = []
pmag3_mem = []
pmag4_mem = []
pmag5_mem = []
pmag6_mem = []

other = 0   #coding checks
out1 = 0
out3 = 0
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['sub'] == 0 or master_cat[counter]['type'] ==0:   #skip no data, stars
        out1 +=1
#    elif master_cat[counter]['sub'] ==1 or master_cat[counter]['sub'] ==3:      #spec & spec_only
#        if master_cat[counter]['member'] ==0:
#            if master_cat[counter]['f_F160W'] >0:
#                smag.append(master_cat[counter]['f_F160W'])            
#                serr.append(master_cat[counter]['e_F160W'])
#                smass.append(master_cat[counter]['lmass'])
#                srad.append(master_cat[counter]['flux_radius'])
#            else: out2 +=1                  # # of spec objects w/o F160 flux
    elif master_cat[counter]['sub'] ==2 or master_cat[counter]['sub']==1:
        if master_cat[counter]['f_F160W'] >0:
            if master_cat[counter]['cluster'] == 1:
                pmag1.append(master_cat[counter]['f_F160W'])             
                perr1.append(master_cat[counter]['e_F160W'])
                prad1.append(master_cat[counter]['flux_radius'])
                if master_cat[counter]['member'] ==0:
                    pmass1.append(master_cat[counter]['lmass'])
                    pmag1_mem.append(master_cat[counter]['f_F160W'])
            elif master_cat[counter]['cluster'] == 2:
                pmag2.append(master_cat[counter]['f_F160W'])             
                perr2.append(master_cat[counter]['e_F160W'])
                prad2.append(master_cat[counter]['flux_radius'])
                if master_cat[counter]['member'] ==0:
                    pmass2.append(master_cat[counter]['lmass'])
                    pmag2_mem.append(master_cat[counter]['f_F160W'])
            elif master_cat[counter]['cluster'] == 3:
                pmag3.append(master_cat[counter]['f_F160W'])             
                perr3.append(master_cat[counter]['e_F160W'])
                prad3.append(master_cat[counter]['flux_radius'])
                if master_cat[counter]['member'] ==0:
                    pmass3.append(master_cat[counter]['lmass'])
                    pmag3_mem.append(master_cat[counter]['f_F160W'])
            elif master_cat[counter]['cluster'] == 4:
                pmag4.append(master_cat[counter]['f_F160W'])             
                perr4.append(master_cat[counter]['e_F160W'])
                prad4.append(master_cat[counter]['flux_radius'])
                if master_cat[counter]['member'] ==0:
                    pmass4.append(master_cat[counter]['lmass'])
                    pmag4_mem.append(master_cat[counter]['f_F160W'])
            elif master_cat[counter]['cluster'] == 5:
                pmag5.append(master_cat[counter]['f_F160W'])             
                perr5.append(master_cat[counter]['e_F160W'])
                prad5.append(master_cat[counter]['flux_radius'])
                if master_cat[counter]['member'] ==0:
                    pmass5.append(master_cat[counter]['lmass'])
                    pmag5_mem.append(master_cat[counter]['f_F160W'])
            elif master_cat[counter]['cluster'] == 6:
                pmag6.append(master_cat[counter]['f_F160W'])             
                perr6.append(master_cat[counter]['e_F160W'])
                prad6.append(master_cat[counter]['flux_radius'])
                if master_cat[counter]['member'] ==0:
                    pmass6.append(master_cat[counter]['lmass'])
                    pmag6_mem.append(master_cat[counter]['f_F160W'])
        else: out3 +=1                  # # of spec objects w/o F160 flux
    else:
        other +=1
    counter+=1
#
# a table of what i'm plotting
names = Column(['phot','no data objects','phot w/o F160'])#'spec','phot','no data objects','spec w/o F160','phot w/o F160'])
values1 = Column([len(pmag1),out1,out3])#[size,len(smag),len(pmag),out1,out2,out3])
values2 = Column([len(pmag2),out1,out3])
values3 = Column([len(pmag3),out1,out3])
values4 = Column([len(pmag4),out1,out3])
values5 = Column([len(pmag5),out1,out3])
values6 = Column([len(pmag6),out1,out3])
F160_data = Table([names,values1,values2,values3,values4,values5,values6],names=('Property','macs0416','macs1149','macs0717','abell370','abell1063','abell2274'))
# convert flux to magnitude
#f_s = np.empty_like(smag, dtype='float16')      #save flux values for signal-to-noise, f_s/_p = spec/phot
f_p1 = np.empty_like(pmag1, dtype='float16')
f_p2 = np.empty_like(pmag2, dtype='float16')
f_p3 = np.empty_like(pmag3, dtype='float16')
f_p4 = np.empty_like(pmag4, dtype='float16')
f_p5 = np.empty_like(pmag5, dtype='float16')
f_p6 = np.empty_like(pmag6, dtype='float16')

#l_s = np.empty_like(smag, dtype='float16')      #save lum values for signal-to-noise
#l_p = np.empty_like(pmag, dtype='float16')
## convert flux to magnitudes for entire catalogue
counter = 0
size = len(pmag1)        #phot
while counter < size:
    f_p1[counter] = pmag1[counter]        #store flux values separately for Fig.2
    pmag1[counter] = -2.5*np.log10(pmag1[counter])+25
    counter+=1
counter = 0
size = len(pmag2)        
while counter < size:
    f_p2[counter] = pmag2[counter]        #store flux values separately for Fig.2
    pmag2[counter] = -2.5*np.log10(pmag2[counter])+25
    counter+=1
counter = 0
size = len(pmag3)        
while counter < size:
    f_p3[counter] = pmag3[counter]        #store flux values separately for Fig.2
    pmag3[counter] = -2.5*np.log10(pmag3[counter])+25
    counter+=1
counter = 0
size = len(pmag4)        
while counter < size:
    f_p4[counter] = pmag4[counter]        #store flux values separately for Fig.2
    pmag4[counter] = -2.5*np.log10(pmag4[counter])+25
    counter+=1
counter = 0
size = len(pmag5)        
while counter < size:
    f_p5[counter] = pmag5[counter]        #store flux values separately for Fig.2
    pmag5[counter] = -2.5*np.log10(pmag5[counter])+25
    counter+=1
counter = 0
size = len(pmag6)        
while counter < size:
    f_p6[counter] = pmag6[counter]        #store flux values separately for Fig.2
    pmag6[counter] = -2.5*np.log10(pmag6[counter])+25
    counter+=1
#
## convert flux to magnitudes for sample only
counter = 0
size = len(pmag1_mem)        #phot
while counter < size:
    pmag1_mem[counter] = -2.5*np.log10(pmag1_mem[counter])+25
    counter+=1
counter = 0
size = len(pmag2_mem)        
while counter < size:
    pmag2_mem[counter] = -2.5*np.log10(pmag2_mem[counter])+25
    counter+=1
counter = 0
size = len(pmag3_mem)        
while counter < size:
    pmag3_mem[counter] = -2.5*np.log10(pmag3_mem[counter])+25
    counter+=1
counter = 0
size = len(pmag4_mem)        
while counter < size:
    pmag4_mem[counter] = -2.5*np.log10(pmag4_mem[counter])+25
    counter+=1
counter = 0
size = len(pmag5_mem)        
while counter < size:
    pmag5_mem[counter] = -2.5*np.log10(pmag5_mem[counter])+25
    counter+=1
counter = 0
size = len(pmag6_mem)        
while counter < size:
    pmag6_mem[counter] = -2.5*np.log10(pmag6_mem[counter])+25
    counter+=1
#
##  compute scatter: Mean Absolute Deviation ('mad') & Standard Deviation ('std_dev')
#
mad = [[0]*6]
mad[0][0] = median_absolute_deviation(pmag1_mem)
mad[0][1] = median_absolute_deviation(pmag2_mem)
mad[0][2] = median_absolute_deviation(pmag3_mem)
mad[0][3] = median_absolute_deviation(pmag4_mem)
mad[0][4] = median_absolute_deviation(pmag5_mem)
mad[0][5] = median_absolute_deviation(pmag6_mem)
std = [[0]*6]
std[0][0] = np.std(pmag1_mem)
std[0][1] = np.std(pmag2_mem)
std[0][2] = np.std(pmag3_mem)
std[0][3] = np.std(pmag4_mem)
std[0][4] = np.std(pmag5_mem)
std[0][5] = np.std(pmag6_mem)
# put in bins
bins = [0,15,17,19,20,21,22,23,23.5,24,24.5,25,25.5,26,26.5,27,27.5,28,29,30]
bins_fine = [0,18,18.5,19,19.25,19.5,19.75,20,20.25,20.5,20.75,21,21.25,21.5,21.75,22,22.25,22.5,22.75,23,23.25,23.5,23.75,24,24.25,24.5,24.75,25,25.25,25.5,25.75,26,26.25,26.5,26.75,27,27.25,27.5,27.75,28,28.5,29,30]
#pmag_hist, mag_bins = np.histogram(pmag, bins=bins,range=(15,30))
#smag_hist, mag_bins = np.histogram(smag, bins=mag_bins,range=(15,30))
pmag_hist1, mag_bins = np.histogram(pmag1, bins=bins,range=(15,30))
pmag_hist2, mag_bins = np.histogram(pmag2, bins=bins,range=(15,30))
pmag_hist3, mag_bins = np.histogram(pmag3, bins=bins,range=(15,30))
pmag_hist4, mag_bins = np.histogram(pmag4, bins=bins,range=(15,30))
pmag_hist5, mag_bins = np.histogram(pmag5, bins=bins,range=(15,30))
pmag_hist6, mag_bins = np.histogram(pmag6, bins=bins,range=(15,30))
mag_midbins = np.empty_like(pmag_hist1, dtype='float64')
#smaller bins around peak
#smag_hist_fine, mag_bins_fine = np.histogram(smag, bins=mag_bins_fine,range=(15,30))
pmag_hist_fine1, mag_bins_fine = np.histogram(pmag1, bins=bins_fine,range=(15,30))
pmag_hist_fine2, mag_bins_fine = np.histogram(pmag2, bins=bins_fine,range=(15,30))
pmag_hist_fine3, mag_bins_fine = np.histogram(pmag3, bins=bins_fine,range=(15,30))
pmag_hist_fine4, mag_bins_fine = np.histogram(pmag4, bins=bins_fine,range=(15,30))
pmag_hist_fine5, mag_bins_fine = np.histogram(pmag5, bins=bins_fine,range=(15,30))
pmag_hist_fine6, mag_bins_fine = np.histogram(pmag6, bins=bins_fine,range=(15,30))
mag_midbins_fine = np.empty_like(pmag_hist_fine1, dtype='float64')
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
#
#
### Plots
#
## Fig. 1A: magnitude histogram, auto binning
plt.close()
F160 = plt.figure(num=1)
plt.suptitle('Apparent Magnitude Histogram (auto bins)')
#cluster 1. macs0416
flux1 = F160.add_subplot(231)
phot1 = flux1.plot(mag_midbins,pmag_hist1,'.m',linewidth=0.5)
plt.plot([26,26],[0,600],':r',label='26',linewidth=1)
plt.plot([27,27],[0,600],':b',label='27',linewidth=1)
plt.plot([28,28],[0,600],':g',label='28',linewidth=1)
plt.xscale('linear')
#plt.gca().invert_xaxis()
plt.xlim=(17,29)
plt.ylabel("count")
plt.yscale('log')
plt.ylim=(1,600)
##plt.title('1. macs0416')
plt.text(17.5,200,'macs0416',fontsize=9)
plt.text(17.5,150,'z = 0.396',fontsize=9)
#plt.legend((m26,m27,m28),('26','27','28'),scatterpoints=1,loc='lower left', frameon=False)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off')
#cluster 2: macs1149
flux2 = F160.add_subplot(232)
phot2 = flux2.plot(mag_midbins,pmag_hist2,'.m',linewidth=0.5)
plt.plot([26,26],[0,600],':r',label='26',linewidth=1)
plt.plot([27,27],[0,600],':b',label='27',linewidth=1)
plt.plot([28,28],[0,600],':g',label='28',linewidth=1)
plt.xscale('linear')
plt.xlim=(17,29)
#plt.ylabel("count")
plt.yscale('log')
plt.ylim=(1,600)
#plt.title('2. macs1149')
plt.text(17.5,200,'macs1149',fontsize=9)
plt.text(17.5,150,'z = 0.543',fontsize=9)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off',labelleft='off')
#cluster 3: macs0717
flux3 = F160.add_subplot(233)
phot3 = flux3.plot(mag_midbins,pmag_hist3,'.m',linewidth=0.5)
plt.plot([26,26],[0,600],':r',label='26',linewidth=1)
plt.plot([27,27],[0,600],':b',label='27',linewidth=1)
plt.plot([28,28],[0,600],':g',label='28',linewidth=1)
plt.xscale('linear')
#plt.gca().invert_xaxis()
plt.xlim=(17,29)
plt.ylabel("count")
flux3.yaxis.set_label_position("right")
plt.yscale('log')
plt.ylim=(1,600)
#plt.title('3. macs0717')
plt.text(17.5,200,'macs0717',fontsize=9)
plt.text(17.5,150,'z = 0.545',fontsize=9)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on',labelleft='off')
#cluster 4: abell370
flux4 = F160.add_subplot(234)
phot4 = flux4.plot(mag_midbins,pmag_hist4,'.m',linewidth=0.5)
plt.plot([26,26],[0,600],':r',label='26',linewidth=1)
plt.plot([27,27],[0,600],':b',label='27',linewidth=1)
plt.plot([28,28],[0,600],':g',label='28',linewidth=1)
plt.xlabel("$m_{F160W}$")
plt.xscale('linear')
#plt.gca().invert_xaxis()
plt.xlim=(17,29)
plt.ylabel("count")
plt.yscale('log')
plt.ylim=(1,600)
#plt.title('4. abell370')
plt.text(17.5,200,'abell370',fontsize=9)
plt.text(17.5,150,'z = 0.375',fontsize=9)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off')
#cluster 5: abell1063
flux5 = F160.add_subplot(235)
phot5 = flux5.plot(mag_midbins,pmag_hist5,'.m',linewidth=0.5)
plt.plot([26,26],[0,600],':r',label='26',linewidth=1)
plt.plot([27,27],[0,600],':b',label='27',linewidth=1)
plt.plot([28,28],[0,600],':g',label='28',linewidth=1)
plt.xlabel("$m_{F160W}$")
plt.xscale('linear')
#plt.gca().invert_xaxis()
plt.xlim=(17,29)
#plt.ylabel("count")
plt.yscale('log')
plt.ylim=(1,600)
#plt.title('5. abell1063')
plt.text(17.5,200,'abell1063',fontsize=9)
plt.text(17.5,150,'z = 0.348',fontsize=9)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off',labelleft='off')
#cluster 6: abell2744
flux6 = F160.add_subplot(236)
phot6 = flux6.plot(mag_midbins,pmag_hist6,'.m',linewidth=0.5)
plt.plot([26,26],[0,600],':r',label='26',linewidth=1)
plt.plot([27,27],[0,600],':b',label='27',linewidth=1)
plt.plot([28,28],[0,600],':g',label='28',linewidth=1)
plt.xlabel("$m_{F160W}$")
plt.xscale('linear')
#plt.gca().invert_xaxis()
plt.xlim=(17,29)
plt.ylabel("count")
flux6.yaxis.set_label_position("right")
plt.yscale('log')
plt.ylim=(1,600)
#plt.title('6. abell2744')
plt.text(17.5,200,'abell2744',fontsize=9)
plt.text(17.5,150,'z = 0.308',fontsize=9)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on',labelleft='off')
#
#
plt.close()
#
## Fig. 1B: magnitude histogram, fine binning
F160f = plt.figure(num=2)
#plt.suptitle('Apparent Magnitude Hist. (fine bins)')
#cluster 1. macs0416
flux1 = F160f.add_subplot(231)
phot1 = flux1.plot(mag_midbins_fine,pmag_hist_fine1,'.m',linewidth=0.5)
#plt.plot([26,26],[0,600],':r',label='26',linewidth=1)
plt.plot([27,27],[0,600],':b',label='27',linewidth=1)
#plt.plot([28,28],[0,600],':g',label='28',linewidth=1)
plt.xscale('linear')
plt.xlim=(17,29)
plt.ylabel("# count")
plt.yscale('log')
plt.ylim=(1,600)
plt.text(17.5,200,'macs0416',fontsize=9)
plt.text(17.5,135,'z = 0.396',fontsize=9)
plt.text(17.5,1.25,'$m_{50\%,comp} = 27$',fontsize=8)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off')
#cluster 2: macs1149
flux2 = F160f.add_subplot(232)
phot2 = flux2.plot(mag_midbins_fine,pmag_hist_fine2,'.m',linewidth=0.5)
#plt.plot([26,26],[0,600],':r',label='26',linewidth=1)
plt.plot([27,27],[0,600],':b',label='27',linewidth=1)
#plt.plot([28,28],[0,600],':g',label='28',linewidth=1)
plt.xscale('linear')
plt.xlim=(17,29)
plt.yscale('log')
plt.ylim=(1,600)
plt.text(17.5,200,'macs1149',fontsize=9)
plt.text(17.5,135,'z = 0.543',fontsize=9)
plt.text(17.5,1.25,'$m_{50\%,comp} = 27$',fontsize=8)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off',labelleft='off')
#cluster 3: macs0717
flux3 = F160f.add_subplot(233)
phot3 = flux3.plot(mag_midbins_fine,pmag_hist_fine3,'.m',linewidth=0.5)
#plt.plot([26,26],[0,600],':r',label='26',linewidth=1)
plt.plot([27,27],[0,600],':b',label='27',linewidth=1)
#plt.plot([28,28],[0,600],':g',label='28',linewidth=1)
plt.xscale('linear')
plt.xlim=(17,29)
plt.ylabel("# count")
flux3.yaxis.set_label_position("right")
plt.yscale('log')
plt.ylim=(1,600)
plt.text(17.5,200,'macs0717',fontsize=9)
plt.text(17.5,135,'z = 0.545',fontsize=9)
plt.text(17.5,1.25,'$m_{50\%,comp} = 27$',fontsize=8)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on',labelleft='off')
#cluster 4: abell370
flux4 = F160f.add_subplot(234)
phot4 = flux4.plot(mag_midbins_fine,pmag_hist_fine4,'.m',linewidth=0.5)
#plt.plot([26,26],[0,600],':r',label='26',linewidth=1)
plt.plot([27,27],[0,600],':b',label='27',linewidth=1)
#plt.plot([28,28],[0,600],':g',label='28',linewidth=1)
plt.xlabel("$m_{F160W}$")
plt.xscale('linear')
plt.xlim=(17,29)
plt.ylabel("# count")
plt.yscale('log')
plt.ylim=(1,600)
plt.text(17.5,200,'abell370',fontsize=9)
plt.text(17.5,135,'z = 0.375',fontsize=9)
plt.text(17.5,1.25,'$m_{50\%,comp} = 27$',fontsize=8)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off')
#cluster 5: abell1063
flux5 = F160f.add_subplot(235)
phot5 = flux5.plot(mag_midbins_fine,pmag_hist_fine5,'.m',linewidth=0.5)
#plt.plot([26,26],[0,600],':r',label='26',linewidth=1)
plt.plot([27,27],[0,600],':b',label='27',linewidth=1)
#plt.plot([28,28],[0,600],':g',label='28',linewidth=1)
plt.xlabel("$m_{F160W}$")
plt.xscale('linear')
plt.xlim=(17,29)
plt.yscale('log')
plt.ylim=(1,600)
plt.text(17.5,200,'abell1063',fontsize=9)
plt.text(17.5,135,'z = 0.348',fontsize=9)
plt.text(17.5,1.25,'$m_{50\%,comp} = 27$',fontsize=8)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off',labelleft='off')
#cluster 6: abell2744
flux6 = F160f.add_subplot(236)
phot6 = flux6.plot(mag_midbins_fine,pmag_hist_fine6,'.m',linewidth=0.5)
#plt.plot([26,26],[0,600],':r',label='26',linewidth=1)
plt.plot([27,27],[0,600],':b',label='27',linewidth=1)
#plt.plot([28,28],[0,600],':g',label='28',linewidth=1)
plt.xlabel("$m_{F160W}$")
plt.xscale('linear')
plt.xlim=(17,29)
plt.ylabel("# count")
flux6.yaxis.set_label_position("right")
plt.yscale('log')
plt.ylim=(1,600)
plt.text(17.5,200,'abell2744',fontsize=9)
plt.text(17.5,135,'z = 0.308',fontsize=9)
plt.text(17.5,1.25,'$m_{50\%,comp} = 27$',fontsize=8)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on',labelleft='off')
plt.show()  
#
#
## Fig.2: (flux/error) vs. magnitude, all in F160W filter
#lines of best fit for both plots   
err_p1 = f_p1/perr1   #calculate flux / delta_fluc = y = 10^(-Ax + B)
err_p2 = f_p2/perr2
err_p3 = f_p3/perr3
err_p4 = f_p4/perr4
err_p5 = f_p5/perr5
err_p6 = f_p6/perr6
#spec/phot bestfit
#def func(x,a,b,c):
#    return a * np.power(10,b*x)+c

#sfit, s_cov = curve_fit(func, smag, err_s,  p0=(10, 0.1,5))
#pfit, p_cov = curve_fit(func, pmag, err_p,  p0=(50, 0.,1000))

#sfit, s_cov = curve_fit(lambda x,a,b,c: np.log10(a)+c,smag,err_s, p0=(1, 0.,10000))
#pfit, p_cov = curve_fit(lambda x,a,b,c: a*np.log10(-b*x)+c,pmag,err_p,  p0=(1, 0.,10000))
#s_yfit = [func(x,sfit[0],sfit[1],sfit[2]) for x in smag]
#p_yfit = [func(x,pfit[0],pfit[1],pfit[2]) for x in pmag]#


#plt.close()
e_F160 = plt.figure(num=3)
plt.suptitle('Flux/Error vs. Apparent Magnitude')
#cluster 1. macs0416
phot_error1 = e_F160.add_subplot(231)
e_phot = phot_error1.scatter(pmag1,err_p1,c='m',marker='.',linewidth=0.5)
#fit_phot = phot_error.scatter(pmag, p_yfit,c='k',marker='.',linewidth=0.05)#,label='Fitted Curve')
plt.plot([0,40],[5,5],':r', linewidth=1)
plt.plot([27,27],[0,6000],':b',label='27',linewidth=1)
#plt.plot(pmag_new,err_p_new,'-k',linewidth=1)
plt.xscale('linear')
plt.xlim=(17,31)
plt.ylabel("$flux_{F160W}/error$")
plt.yscale('log')
plt.ylim=(1,5000)
#plt.legend((e_spec,e_phot),('spec','phot'),scatterpoints=1,loc='best', frameon=False)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off', labelleft='on')
plt.text(24,1000,'macs0416',fontsize=8)
#cluster 2. macs1149
phot_error2 = e_F160.add_subplot(232)
e_phot = phot_error2.scatter(pmag2,err_p2,c='m',marker='.',linewidth=0.5)
#fit_phot = phot_error.scatter(pmag, p_yfit,c='k',marker='.',linewidth=0.05)#,label='Fitted Curve')
plt.plot([0,40],[5,5],':r', linewidth=1)
plt.plot([27,27],[0,6000],':b',label='27',linewidth=1)
#plt.plot(pmag_new,err_p_new,'-k',linewidth=1)
plt.xscale('linear')
plt.xlim=(17,31)
#plt.ylabel("$flux_{F160W}/error$")
plt.yscale('log')
plt.ylim=(1,5000)
#plt.legend((e_spec,e_phot),('spec','phot'),scatterpoints=1,loc='best', frameon=False)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off', labelleft='off')
plt.text(24,1000,'macs1149',fontsize=8)
#cluster 3. macs0717
phot_error3 = e_F160.add_subplot(233)
e_phot = phot_error3.scatter(pmag3,err_p3,c='m',marker='.',linewidth=0.5)
#fit_phot = phot_error.scatter(pmag, p_yfit,c='k',marker='.',linewidth=0.05)#,label='Fitted Curve')
plt.plot([0,40],[5,5],':r', linewidth=1)
plt.plot([27,27],[0,6000],':b',label='27',linewidth=1)
#plt.plot(pmag_new,err_p_new,'-k',linewidth=1)
plt.xscale('linear')
plt.xlim=(17,31)
plt.ylabel("$flux_{F160W}/error$")
phot_error3.yaxis.set_label_position("right")
plt.yscale('log')
plt.ylim=(1,5000)
#plt.legend((e_spec,e_phot),('spec','phot'),scatterpoints=1,loc='best', frameon=False)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on', labelleft='off')
plt.text(24,1000,'macs0717',fontsize=8)
#cluster 4. abell370
phot_error4 = e_F160.add_subplot(234)
e_phot = phot_error4.scatter(pmag4,err_p4,c='m',marker='.',linewidth=0.5)
#fit_phot = phot_error.scatter(pmag, p_yfit,c='k',marker='.',linewidth=0.05)#,label='Fitted Curve')
plt.plot([0,40],[5,5],':r', linewidth=0.8)
plt.plot([27,27],[0,6000],':b',label='27',linewidth=1)
#plt.plot(pmag_new,err_p_new,'-k',linewidth=1)
plt.xlabel("$m_{F160W}$")
plt.xscale('linear')
plt.xlim=(17,31)
plt.ylabel("$flux_{F160W}/error$")
plt.yscale('log')
plt.ylim=(1,5000)
#plt.legend((e_spec,e_phot),('spec','phot'),scatterpoints=1,loc='best', frameon=False)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off', labelleft='on')
plt.text(24,1000,'abell370',fontsize=8)
#cluster 5. abell1063
phot_error5 = e_F160.add_subplot(235)
e_phot = phot_error5.scatter(pmag5,err_p5,c='m',marker='.',linewidth=0.5)
#fit_phot = phot_error.scatter(pmag, p_yfit,c='k',marker='.',linewidth=0.05)#,label='Fitted Curve')
plt.plot([0,40],[5,5],':r', linewidth=0.8)
plt.plot([27,27],[0,6000],':b',label='27',linewidth=1)
#plt.plot(pmag_new,err_p_new,'-k',linewidth=1)
plt.xlabel("$m_{F160W}$")
plt.xscale('linear')
plt.xlim=(17,31)
#plt.ylabel("$flux_{F160W}/error$")
plt.yscale('log')
plt.ylim=(1,5000)
#plt.legend((e_spec,e_phot),('spec','phot'),scatterpoints=1,loc='best', frameon=False)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off', labelleft='off')
plt.text(24,1000,'abell1063',fontsize=8)
#cluster 6. abell2744
phot_error6 = e_F160.add_subplot(236)
e_phot = phot_error6.scatter(pmag6,err_p6,c='m',marker='.',linewidth=0.5)
#fit_phot = phot_error.scatter(pmag, p_yfit,c='k',marker='.',linewidth=0.05)#,label='Fitted Curve')
plt.plot([0,40],[5,5],':r', linewidth=0.8)
plt.plot([27,27],[0,6000],':b',label='27',linewidth=1)
#plt.plot(pmag_new,err_p_new,'-k',linewidth=1)
plt.xlabel("$m_{F160W}$")
plt.xscale('linear')
plt.xlim=(17,31)
plt.ylabel("$flux_{F160W}/error$")
phot_error6.yaxis.set_label_position("right")
plt.yscale('log')
plt.ylim=(1,5000)
#plt.legend((e_spec,e_phot),('spec','phot'),scatterpoints=1,loc='best', frameon=False)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on', labelleft='off')
plt.text(24,1000,'abell2744',fontsize=8)
#
plt.show()  
#
#
##  Fig.3: mass-to-luminosity plot(s); plotting magnitude in place of luminosity
##  (I) full sample
##      (i-vi) plots for each individual cluster
#
#  (I) convert mag to lum (mult by surface area)
#counter = 0
#size = len(l_s)
#while counter < size: 
#    l_s[counter] = f_s[counter]*4*np.pi*srad[counter]**2
#    counter +=1
#counter = 0
#size = len(l_p)
#while counter < size: 
#    l_p[counter] = f_p[counter]*4*np.pi*prad[counter]**2
#    counter +=1
#determine limiting luminosity from limiting magnitude of 26.3

#
plt.close()
mlum = plt.figure(num=4)
#plt.suptitle('Magnitude vs. Mass')
#cluster 3: macs 0717
plum1 = mlum.add_subplot(231)
#plt.plot(l_p,pmass,'.m',linewidth=0.5)
plt.plot(pmag1_mem,pmass1,'.c',linewidth=0.5)
plt.plot([26.3,26.3],[0,15],'--r', linewidth=0.8)
plt.plot([15,35],[7.5,7.5],'-.b', linewidth=0.8)
plt.plot([15,35],[7.1,7.1],'-.b', linewidth=0.8)
plt.xscale('linear')
plt.xlim=(17,30)
plt.ylabel('$log(M/M_{\odot})$')
plt.yscale('linear')
plt.ylim=(5,12.5)
#plt.legend((smass,pmass),('spec','phot'),scatterpoints=1,loc='upper left', frameon=False)
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off', labelleft='on')
plt.minorticks_on()
plt.text(22,11,'macs0416',fontsize=8)
plt.text(22,10.5,'$M/L_{max}: 7.5$',fontsize=8)
plt.text(17.5,6.1,'z = 0.396',fontsize=9)
plt.text(17.5,5.4,'$\sigma = %s$' %('%.2f'%std[0][0]),fontsize=9)
#cluster 2: macs 1149
plum2 = mlum.add_subplot(232)
#plt.title('$m_{F160W, lim}=26.3$',fontsize=10)
plt.plot(pmag2_mem,pmass2,'.c',linewidth=0.5)
plt.plot([26.3,26.3],[0,15],'--r', linewidth=0.8)
plt.plot([15,35],[7.8,7.8],'-.b', linewidth=0.8)
plt.plot([15,35],[6.8,6.8],'-.b', linewidth=0.8)
plt.xscale('linear')
plt.xlim=(17,30)
plt.yscale('linear')
plt.ylim=(5,12.5)
#plt.legend((smass,pmass),('spec','phot'),scatterpoints=1,loc='upper left', frameon=False)
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off', labelleft='off')
plt.minorticks_on()
plt.text(22,11,'macs1149',fontsize=8)
plt.text(22,10.5,'$M/L_{max}: 7.8$',fontsize=8)
plt.text(17.5,6.1,'z = 0.543',fontsize=9)
plt.text(17.5,5.4,'$\sigma = %s$' %('%.2f'%std[0][1]),fontsize=9)
#cluster 3: macs 0717
plum3 = mlum.add_subplot(233)
#plt.plot(l_p,pmass,'.m',linewidth=0.5)
plt.plot(pmag3_mem,pmass3,'.c',linewidth=0.5)
plt.plot([26.3,26.3],[0,15],'--r', linewidth=0.8)
plt.plot([15,35],[8,8],'-.b', linewidth=0.8)
plt.plot([15,35],[6.9,6.9],'-.b', linewidth=0.8)
plt.xscale('linear')
plt.xlim=(17,30)
plt.ylabel('$log(M/M_{\odot})$')
plt.yscale('linear')
plum3.yaxis.set_label_position("right")
plt.ylim=(5,12.5)
#plt.legend((smass,pmass),('spec','phot'),scatterpoints=1,loc='upper left', frameon=False)
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on', labelleft='off')
plt.minorticks_on()
plt.text(22,11,'macs0717',fontsize=8)
plt.text(22,10.5,'$M/L_{max}: 8.0$',fontsize=8)
plt.text(17.5,6.1,'z = 0.545',fontsize=9)
plt.text(17.5,5.4,'$\sigma = %s$' %('%.2f'%std[0][2]),fontsize=9)
#cluster 4: abell370
plum4 = mlum.add_subplot(234)
#plt.plot(l_p,pmass,'.m',linewidth=0.5)
plt.plot(pmag4_mem,pmass4,'.c',linewidth=0.5)
plt.plot([26.3,26.3],[0,15],'--r', linewidth=0.8)
plt.plot([15,35],[7.5,7.5],'-.b', linewidth=0.8)
plt.plot([15,35],[6.4,6.4],'-.b', linewidth=0.8)
plt.xlabel('$m_{F160W}$')
plt.xscale('linear')
plt.xlim=(17,30)
plt.ylabel('$log(M/M_{\odot})$')
plt.yscale('linear')
plum4.yaxis.set_label_position("left")
plt.ylim=(5,12.5)
#plt.legend((smass,pmass),('spec','phot'),scatterpoints=1,loc='upper left', frameon=False)
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off', labelleft='on')
plt.minorticks_on()
plt.text(22,11,'abell370',fontsize=8)
plt.text(22,10.5,'$M/L_{max}: 7.5$',fontsize=8)
plt.text(17.5,6.5,'z = 0.375',fontsize=9)
plt.text(17.5,5.5,'$\sigma = %s$' %('%.2f'%std[0][3]),fontsize=9)
#cluster 5: abell 1063
plum5 = mlum.add_subplot(235)
#plt.plot(l_p,pmass,'.m',linewidth=0.5)
plt.plot(pmag5_mem,pmass5,'.c',linewidth=0.5)
plt.plot([26.3,26.3],[0,15],'--r', linewidth=0.8)
plt.plot([15,35],[7.4,7.4],'-.b', linewidth=0.8)
plt.plot([15,35],[6.3,6.3],'-.b', linewidth=0.8)
plt.xlabel('$m_{F160W}$')
plt.xscale('linear')
plt.xlim=(17,30)
plt.yscale('linear')
plt.ylim=(5,12.5)
#plt.legend((smass,pmass),('spec','phot'),scatterpoints=1,loc='upper left', frameon=False)
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off', labelleft='off')
plt.minorticks_on()
plt.text(22,11,'abell1063',fontsize=8)
plt.text(22,10.5,'$M/L_{max}: 7.4$',fontsize=8)
plt.text(17.5,6.5,'z = 0.348',fontsize=9)
plt.text(17.5,5.5,'$\sigma = %s$' %('%.2f'%std[0][4]),fontsize=9)
#cluster 6: abell 2744
plum6 = mlum.add_subplot(236)
#plt.plot(l_p,pmass,'.m',linewidth=0.5)
plt.plot(pmag6_mem,pmass6,'.c',linewidth=0.5)
plt.plot([26.3,26.3],[0,15],'--r', linewidth=0.8)
plt.plot([15,35],[7.3,7.3],'-.b', linewidth=0.8)
plt.plot([15,35],[6.1,6.1],'-.b', linewidth=0.8)
plt.xlabel('$m_{F160W}$')
plt.xscale('linear')
plt.xlim=(17,30)
plt.ylabel('$log(M/M_{\odot})$')
plt.yscale('linear')
plum6.yaxis.set_label_position("right")
plt.ylim=(5,12.5)
#plt.legend((smass,pmass),('spec','phot'),scatterpoints=1,loc='upper left', frameon=False)
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on', labelleft='off')
plt.minorticks_on()
plt.text(22,11,'abell2744',fontsize=8)
plt.text(22,10.5,'$M/L_{max}: 7.3$',fontsize=8)
plt.text(17.5,6.5,'z = 0.308',fontsize=9)
plt.text(17.5,5.5,'$\sigma = %s$' %('%.2f'%std[0][5]),fontsize=9)
plt.show() 

