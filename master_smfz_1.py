#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 14:28:35 2017

@author: gsarrouh
"""
## This program will plot the Stellar Mass Function (SMF) and Luminosity 
## ("Schechter") function for the master_dadta file of all six clusters
## (most current version: 'master_data_7); two plots (SF & Q) segregated 
## between spectroscopic and photometric subsamples
#
## v3 includes code to allow for more direct customization of binning, without 
## need to change size of new lists
## v4 removes many of the lists/arrays which track bins values by galaxy type 
## for plotting purposes (e.g. SF_field bin value & corresponding SF_bins mass 
## range stored in new list) and instead simply calls the relevant row of the 
## array idirectly
#
## v4a plots SMF for master_data_7a, with the modified z_phot cut (which 
## reduces the number of phot members selected); 
## lower limit range of bins modified slightly from 5.33 to 5.80 to reflect 
## min value in cluster list (as opposed to 5.33 for field);  # bins set to 20
#
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy
from astropy.table import Table
from astropy.table import Column

# plot stellar mass function:
## (i) collect masses into sorted arrays for SF & Q;
## (ii) bin them as in a histogram
## (iii) find midpoints of bins s.t. the # count from hist & # of bins are equal
## (iv) add error bars to scatter plot
## (v) fit a SCHECHTER function to data
## (vi) plot that shit
#
#
## SECTION (i): collect relevant objects into a single array in order to plot 
## historgram main arrays: (SF/Q)_hist, (SF/Q)_field_hist, creates list of 
## samples to be binned & plotted as histogram/scatterplot
#
BCG_cluster = [[0]*2]
SF_hist = Column(np.zeros(shape=(np.sum(phot_mems)+np.sum(smem),1), dtype='float64'))
SF_index = 0
Q_hist = Column(np.empty(shape=(np.sum(phot_memq)+np.sum(qmem),1), dtype='float64'))
Q_index = 0
other = 0
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==1:      #secure cluster, SF
        SF_hist[SF_index][0] = master_cat[counter]['lmass']
        SF_index +=1
        if master_cat[counter]['BCG'] ==2:
            BCG_cluster[0][0] +=1
    elif master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==2: #secure cluster, Q
        Q_hist[Q_index][0] = master_cat[counter]['lmass']
        Q_index +=1
        if master_cat[counter]['BCG'] ==2:
            BCG_cluster[0][1] +=1
    else:
        other +=1
    counter +=1        
#
## Field sample
BCG_field = [[0]*2]
SF_field_hist = Column(np.zeros(shape=(np.sum(phot_fields)+np.sum(sfield),1), dtype='float64'))
SF_index = 0
Q_field_hist = np.empty(shape=(np.sum(phot_fieldq)+np.sum(qfield),1), dtype='float64')
Q_index = 0
other = 0
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['member'] == 1 and master_cat[counter]['type'] ==1:      #secure field, SF
        SF_field_hist[SF_index][0] = master_cat[counter]['lmass']
        SF_index +=1
        if master_cat[counter]['BCG'] ==2:
            BCG_field[0][0] +=1
    elif master_cat[counter]['member'] == 1 and master_cat[counter]['type'] ==2: #secure field, Q
        Q_field_hist[Q_index][0] = master_cat[counter]['lmass']
        Q_index +=1
        if master_cat[counter]['BCG'] ==2:
            BCG_field[0][1] +=1
    else:
        other +=1
    counter +=1     
#
## to remove entries which return "NAN" (not a number). no idea why this would 
## happen, but this is just a patch so the program runs. have set all such 
## values equal to 0, and will be removed from histogram as they are out of range
counter=0           
size=len(SF_hist)
while counter<size:
    if np.isnan(SF_hist[counter]):
        SF_hist[counter] = 0
    counter+=1
# same for Q
counter=0           
size=len(Q_hist)
while counter<size:
    if np.isnan(Q_hist[counter]):
        Q_hist[counter] = 0
    counter+=1
#
# now same for field sample
#
counter=0           
size=len(SF_field_hist)
while counter<size:
    if np.isnan(SF_field_hist[counter]):
        SF_field_hist[counter] = 0
    counter+=1
# same for Q
counter=0           
size=len(Q_field_hist)
while counter<size:
    if np.isnan(Q_field_hist[counter]):
        Q_field_hist[counter] = 0     
    counter+=1
#
#
## SECTION (ii): sort objects into histogram bins for both SF & Q populations, then sum 
## for 'total' population. use for total pop plot, and for relative fractions
#
## cluster populations arrays: (SF/Q/total)_smf & (SF/Q/total)_field_smf
#
aa = []
aa.append(np.linspace(5.8,8,num=5))
aa.append(np.linspace(8,10,num=5))
SF_smf, SF_bins = np.histogram(SF_hist, bins=aa)#,range=(5.80,12.14))
Q_smf, Q_bins = np.histogram(Q_hist, bins=SF_bins)      #use same bins for both populations to facilitate comparison
total_smf = np.empty_like(SF_smf, dtype='float64')
frac_smf = np.empty(shape=(2,len(SF_bins)-1),dtype='float64')       #1st row: SF fraction; 2nd row Q fraction
total_smf = SF_smf + Q_smf
counter = 0
size = len(SF_smf)#size = frac_smf.shape[1]
while counter < size:             #loop to skip bins with no entires. will probably need to reduce range for completeness
    if total_smf[counter] == 0:
        frac_smf[0][counter] = 0
        frac_smf[1][counter] = 0
    else: 
        frac_smf[0][counter] =  SF_smf[counter] / total_smf[counter]
        frac_smf[1][counter] =  Q_smf[counter] / total_smf[counter]
    counter +=1 
#
## field populations
#
SF_field_smf, SF_field_bins = np.histogram(SF_field_hist, bins=SF_bins,range=(5.80,12.14))
Q_field_smf, Q_field_bins = np.histogram(Q_field_hist, bins=SF_bins)      #use same bins for both populations to facilitate comparison
total_field_smf = np.empty_like(SF_smf, dtype='int8')
frac_field_smf = np.empty(shape=(2,len(SF_bins)-1),dtype='float64')       #1st row: SF fraction; 2nd row Q fraction
total_field_smf = SF_field_smf + Q_field_smf
#
#
## calculate bin fractions
counter = 0
size = len(SF_field_smf)
while counter < size:             #loop to skip bins with no entires. will probably need to reduce range for completeness
    if total_field_smf[counter] == 0:
        frac_field_smf[0][counter] = 0
        frac_field_smf[1][counter] = 0
    else: 
        frac_field_smf[0][counter] =  SF_field_smf[counter] / total_field_smf[counter]
        frac_field_smf[1][counter] =  Q_field_smf[counter] / total_field_smf[counter]
    counter +=1 
#
#
## (iii) find midpoints of spec bins, and split hist. into "x" & "y" for scatter plotting
#
## find midpoint of hist. bins. all populations have been binned identically, 
## so the one 'midbin' will serve for all data arrays to be plotted
#
SF_midbins = np.empty_like(SF_smf, dtype='float64')
#
size = len(SF_bins)-1
counter = 0
while counter < size:
    SF_midbins[counter] = (SF_bins[counter] + SF_bins[counter+1])/2
    counter +=1
#
## (iv) error bars & relative fractions
## Poissonian error bars will be added to the spec. sample only; phot sample 
## requires MCMC error to account for cluster membership correction
#
## Method: treat each bin as its own Poisson distribution, w/ expectation equal
## to count in each respective bin. the error is then the sqrt of that count
#
## *******this is a test code to plot error bars, which will be applied to the 
## entire field sample. adjustment to be made later to isolate spec sample alone
#
## I. error values by sample type, cluster & field
SF_error = np.sqrt(SF_smf)        #assign errors by population as sqrt of count in each mass bin
Q_error = np.sqrt(Q_smf)
SF_field_error = np.sqrt(SF_field_smf)        #assign errors by population as sqrt of count in each mass bin
Q_field_error = np.sqrt(Q_field_smf)
total_field_error = np.sqrt(total_field_smf)
#frac_error = np.sqrt(frac_smf)
#frac_field_error = np.sqrt(frac_field_smf)
#frac_error [0] = (SF_error/SF_smf)
#frac_error [1] = (Q_error/Q_smf)
frac_field_error [0] = (SF_field_error/SF_field_smf)
frac_field_error [1] = (Q_field_error/Q_field_smf)
#
#
###     NOTE: relative fraction calculated above at end of Seciont (ii), and are
###     stored in arrays "frac_smf" & "frac_field_smf"
#
#
## (v) fit a SCHECHTER function to scatter plot
#
#
#
#
## SECTION (vi): create plot
## upper: SMF for cluster, field
## lower: fractions of SF/Q in cluster, field
#
## SMF
## Cluster
SMF = plt.figure(num=1)
gs = gridspec.GridSpec(2,2, width_ratios=[1,1], height_ratios=[2,1])   #make a tiled-plot like vdB2013 w/ fractions below, this line sets the proporitons of plots in the figure
#
cluster = plt.subplot(gs[0])      
plt.plot(SF_midbins,Q_smf,'.r',linewidth=0.5)#Qx_new, Qy_new,'-r',linewidth=0.5)
plt.plot(SF_midbins,total_smf,'.k',linewidth=0.5)#totalx_new, totaly_new,'-k',linewidth=0.5)
plt.errorbar(SF_midbins,SF_smf,fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)#yerr=SF_error,
plt.xscale('linear')
#plt.xlabel('log(M/M_sol)')
plt.xlim(5.75,12.5)
plt.yscale('log')
cluster.minorticks_on()
cluster.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off',labelbottom='off')
cluster.yaxis.set_label_position("left")
plt.ylabel('logarthmic # count')
plt.ylim(0,750)
plt.title('Cluster SMF')
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')
#
## Field
field = plt.subplot(gs[1])      #make a tiled-plot like vdB2013 w/ fractions below
plt.errorbar(SF_midbins,SF_field_smf,yerr=SF_field_error, fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
plt.errorbar(SF_midbins,Q_field_smf,yerr=Q_field_error, fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
plt.errorbar(SF_midbins,total_field_smf,yerr=total_field_error, fmt='.k',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
plt.xscale('linear')
#plt.xlabel('log(M/M_sol)')
plt.xlim(5.75,12.5)
plt.yscale('log')
plt.ylabel('logarthmic # count')
#plt.ylim(15,1200)
field.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on',labelleft='off',labelbottom='off')
field.minorticks_on()
field.yaxis.set_label_position("right")
plt.title('Field SMF')
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')
#plt.text()      #this could replace the legend to avoid visual clutter
#
## cluster fraction
fr_cl = plt.subplot(gs[2])    
plt.errorbar(SF_midbins,frac_smf[0],yerr=frac_error[0], fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
plt.errorbar(SF_midbins,frac_smf[1],yerr=frac_error[1], fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
plt.xscale('linear')
plt.xlabel('log(M/M_sol)')
plt.xlim(5.75,12.5)
plt.yscale('linear')
plt.ylim(-0.1,1.1)
fr_cl.minorticks_on()
fr_cl.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off')
fr_cl.yaxis.set_label_position("left")
plt.ylabel('Galaxy type fraction')
#plt.ylim(1,190)
#plt.title('relative fraction SF & Q')
#
## field fraction
fr_fld = plt.subplot(gs[3])     
plt.errorbar(SF_midbins,frac_field_smf[0],yerr=frac_field_error[0], fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
plt.errorbar(SF_midbins,frac_field_smf[1],yerr=frac_field_error[1], fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
plt.xscale('linear')
plt.xlabel('log(M/M_sol)')
plt.xlim(5.75,12.5)
plt.yscale('linear')
plt.ylim(-0.1,1.1)
fr_fld.minorticks_on()
fr_fld.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on',labelleft='off')
fr_fld.yaxis.set_label_position("right")
plt.ylabel('Galaxy type fraction')


plt.show()
#plt.text(0.25,2,'Passive',fontsize=10)
#plt.text(1.5,0.75,'Star-Forming',fontsize=10)
#
#
## create summary table of binning info - "bin_stats"
bin_stats = Table([SF_midbins,total_smf,SF_smf,Q_smf,total_field_smf,SF_field_smf,Q_field_smf],names=('Bin mass','Total cluster','SF cluster','Q cluster','Total field','SF field','Q field'))
