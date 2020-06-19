#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 21:32:30 2017

@author: gsarrouh
"""
#####################  data_completeness_3  #####################
#
###  This program creates plots related to completeness of the data.
#
##  v3: calcuates lines of best fit using scipy.optimize.curve_fit for non-linear regression
#
##  Fig.1: apparent magnitude histogram (F160 filter, i.e. 1600nm), 'F160_data' for data
##  Fig.2: (flux/error) vs. magnitude, all in F160W filter; determine limiting mag: (27-0.7)= 26.3
##  Fig.3: mass-to-luminosity plot; determine limiting mass: 
#
# 
from scipy.optimize import curve_fit
## Fig 1: plot histogram of F160 magnitude (show spec vs phot data)
#
# create list of flux values at 1600nm
# note: plotting all objects for which there are data (i.e. excluding only 'nodata')
smag = []   #SF magnitude
serr = []   #SF error
smass = []  #SF mass
srad = []   #SF radius
pmag = []   #same quantities for phot
perr = []
pmass = []
prad = []
other = 0   #coding checks
out1 = 0
out2 = 0
out3 = 0
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['sub'] == 0:   #skip no data
        out1 +=1
    elif master_cat[counter]['sub'] ==1 or master_cat[counter]['sub'] ==3:      #spec & spec_only
        if master_cat[counter]['member'] ==0:
            if master_cat[counter]['f_F160W'] >0:
                smag.append(master_cat[counter]['f_F160W'])            
                serr.append(master_cat[counter]['e_F160W'])
                smass.append(master_cat[counter]['lmass'])
                srad.append(master_cat[counter]['flux_radius'])
            else: out2 +=1                  # # of spec objects w/o F160 flux
    elif master_cat[counter]['sub'] ==2:
        if master_cat[counter]['member'] ==0:
            if master_cat[counter]['f_F160W'] >0:
                pmag.append(master_cat[counter]['f_F160W'])             
                perr.append(master_cat[counter]['e_F160W'])
                pmass.append(master_cat[counter]['lmass'])
                prad.append(master_cat[counter]['flux_radius'])
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
f_s = np.empty_like(smag, dtype='float16')      #save flux values for signal-to-noise, f_s/_p = spec/phot
f_p = np.empty_like(pmag, dtype='float16')
l_s = np.empty_like(smag, dtype='float16')      #save lum values for signal-to-noise
l_p = np.empty_like(pmag, dtype='float16')
counter = 0
size = len(smag)        #spec
while counter < size:
    f_s[counter] = smag[counter]        #store flux values separately for Fig.2
    smag[counter] = -2.5*np.log10(smag[counter])+25
    counter+=1
counter = 0
size = len(pmag)        #phot
while counter < size:
    f_p[counter] = pmag[counter]        #store flux values separately for Fig.2
    pmag[counter] = -2.5*np.log10(pmag[counter])+25
    counter+=1
# convert 'flux_radius' (i.e. flux half-radius) to radius of star
counter = 0
size = len(srad)
while counter < size:
    srad[counter] = srad[counter] / np.sqrt(2)
    counter+=1
counter = 0
size = len(prad)
while counter < size:
    prad[counter] = prad[counter] / np.sqrt(2)
    counter+=1
# put in bins
bins = [0,15,17,19,20,21,22,23,23.5,24,24.5,25,25.5,26,26.5,27,27.5,28,29,30]
bins_fine = [0,18,18.5,19,19.25,19.5,19.75,20,20.25,20.5,20.75,21,21.25,21.5,21.75,22,22.25,22.5,22.75,23,24,25,25.25,25.5,25.75,26,26.25,26.5,26.75,27,27.25,27.5,27.75,28,28.5,29,30]
smag_hist, mag_bins = np.histogram(smag, bins='auto',range=(15,30))
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
## Plots
#
## Fig. 1: magnitude histogram
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
## Fig.2: (flux/error) vs. magnitude, all in F160W filter
#lines of best fit for both plots
err_s = f_s/serr    #calculate flux / delta_fluc = y = 10^(-Ax + B)
err_p = f_p/perr

#spec/phot bestfit
def func(x,a,b,c):
    return a * np.power(10,b*x)+c

sfit, s_cov = curve_fit(func, smag, err_s,  p0=(10, 0.1,5))
pfit, p_cov = curve_fit(func, pmag, err_p,  p0=(50, 0.,1000))

#sfit, s_cov = curve_fit(lambda x,a,b,c: np.log10(a)+c,smag,err_s, p0=(1, 0.,10000))
#pfit, p_cov = curve_fit(lambda x,a,b,c: a*np.log10(-b*x)+c,pmag,err_p,  p0=(1, 0.,10000))
s_yfit = [func(x,sfit[0],sfit[1],sfit[2]) for x in smag]
p_yfit = [func(x,pfit[0],pfit[1],pfit[2]) for x in pmag]#


#plt.close()
e_F160 = plt.figure(num=2)
plt.suptitle('Flux/Error vs. Apparent Magnitude')
spec_error = e_F160.add_subplot(121)
e_spec = spec_error.scatter(smag,err_s,c='g',marker='.',linewidth=0.5)
fit_spec = plt.scatter(smag, s_yfit,c='k',marker='.',linewidth=0.05)#,label='Fitted Curve')
plt.plot([0,40],[5,5],':r', linewidth=1)
#plt.plot(smag_new,err_s_new,'-k',linewidth=1)
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
e_phot = phot_error.scatter(pmag,err_p,c='m',marker='.',linewidth=0.5)
fit_phot = phot_error.scatter(pmag, p_yfit,c='k',marker='.',linewidth=0.05)#,label='Fitted Curve')
plt.plot([0,40],[5,5],':r', linewidth=0.8)
#plt.plot(pmag_new,err_p_new,'-k',linewidth=1)
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
#
plt.show()  
#
#
##  Fig.3: mass-to-luminosity plot(s)
##  (I) full sample
##      (i-vi) plots for each individual cluster
#
#  (I) convert mag to lum (mult by surface area)
counter = 0
size = len(l_s)
while counter < size: 
    l_s[counter] = f_s[counter]*4*np.pi*srad[counter]**2
    counter +=1
counter = 0
size = len(l_p)
while counter < size: 
    l_p[counter] = f_p[counter]*4*np.pi*prad[counter]**2
    counter +=1
#determine limiting luminosity from limiting magnitude of 26.3

#
mlum = plt.figure(num=3)
plt.suptitle('Magnitude vs. Luminosity')
#spec
slum = mlum.add_subplot(121)
#plt.plot(l_s,smass,'.g',linewidth=0.5)
plt.plot(smass,smag,'.g',linewidth=0.5)
plt.plot([0,15],[26.3,26.3],':r', linewidth=0.8)
#plt.xlabel("$luminosity_{F160W}$")
plt.xlabel('$log(M/M_{\odot})$)')
plt.xscale('linear')
#plt.gca().invert_xaxis()
#plt.xlim(15,40)
plt.ylabel('magnitude')
plt.yscale('linear')
#plt.ylim(5,15)
plt.title('spec')
#plt.legend((smass,pmass),('spec','phot'),scatterpoints=1,loc='upper left', frameon=False)
#plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off')
plt.minorticks_on()
#phot
plum = mlum.add_subplot(122)
#plt.plot(l_p,pmass,'.m',linewidth=0.5)
plt.plot(pmass,pmag,'.m',linewidth=0.5)
plt.plot([0,15],[26.3,26.3],':r', linewidth=0.8)
plt.plot([6,6],[15,35],':k', linewidth=0.8)
plt.plot([8,8],[15,35],':k', linewidth=0.8)
plt.xlabel('$log(M/M_{\odot})$)')
plt.xscale('linear')
#plt.gca().invert_xaxis()
#plt.xlim(15,40)
plt.ylabel('$luminosity_{F160W}$')
plt.yscale('linear')
plum.yaxis.set_label_position("right")
#plt.ylim(5,15)
plt.title('phot')
#plt.legend((smass,pmass),('spec','phot'),scatterpoints=1,loc='upper left', frameon=False)
#plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on', labelleft='off')
plt.minorticks_on()
plt.show() 
#  (ii) 

## Fig. 4: mass histrograms by cluster for false positive negative to determine mass bins for completeness corrections
#


