#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 14:28:35 2017

@author: gsarrouh
"""
################## master_smfz_6 ##################
## This program will plot the Stellar Mass Function (SMF) and Luminosity 
## ("Schechter") function for the master_dadta file of all six clusters
## (most current version: 'master_data_7); two plots (SF & Q) segregated 
## between spectroscopic and photometric subsamples
#
## v2 includes the parallel fields data; 
## v3 commented out to produce conference 
##    plot. to re-insert, remove all #s and change plot panels to include || 
##    field in cenre column (positions 1 & 4)
## v4 creates SMFs by cluster, and calculates completeness correction by cluster for false pos/neg
##    and attempts compute the completeness correction on a cluster-by-cluster basis.
##    this approach was abandoned in favour of a single set of correciton factors 
##    for the entire sample, justified by the fact that mass completeness limits are nearly the same for all clusters. 
## v5 removes individual cluster code for mass correction in section (v), plots
##    correction factors for entire sample treated as single population
## v6 moves all totals and SF/Q fractions until AFTER the correction factors 
##    have been made, as this changes the counts in each hist & associated poissonian error
#
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy
from astropy.table import Table
from astropy.table import Column

# plot stellar mass function:
## (i)      collect masses into sorted arrays for SF & Q;
## (ii)     bin them as in a histogram
## (iii)    find midpoints of bins s.t. the # count from hist & # of bins are equal
## (iv)     add error bars to scatter plot
## (v)      false pos/neg completess correction
## (vi)     fit a SCHECHTER function to data
## (vii)    plot that shit
#
#
## SECTION (i): collect relevant objects into a single array in order to plot 
## historgram main arrays: (SF/Q)_hist, (SF/Q)_field_hist, creates list of 
## samples to be binned & plotted as histogram/scatterplot
#
## Cluster sample
a = 0
b = 0
BCG_cluster = [[0]*2]
SF_hist = []
Q_hist = []
other = 0
counter = 0
SF_1 = []       #lists of SF/Q by cluster, for flase pos/neg spectroscopic completeness corection
SF_2 = []
SF_3 = []
SF_4 = []
SF_5 = []
SF_6 = []
Q_1 = []
Q_2 = []
Q_3 = []
Q_4 = []
Q_5 = []
Q_6 = []
size = len(master_cat)
while counter < size:
    if master_cat[counter]['lmass'] > 7.3:
        if master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==1:      #secure cluster, SF
            SF_hist.append(master_cat[counter]['lmass'])
            if master_cat[counter]['BCG'] ==2:
                BCG_cluster[0][0] +=1
            if master_cat[counter]['cluster'] ==1:          #by cluster for mass completeness correction
                SF_1.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==2:
                SF_2.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==3:
                SF_3.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==4:
                SF_4.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==5:
                SF_5.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==6:
                SF_6.append(master_cat[counter]['lmass'])
        elif master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==2: #secure cluster, Q
            Q_hist.append(master_cat[counter]['lmass'])
            if master_cat[counter]['BCG'] ==2:              #by cluster for mass completeness correction
                BCG_cluster[0][1] +=1
            if master_cat[counter]['cluster'] ==1:
                Q_1.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==2:
                Q_2.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==3:
                Q_3.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==4:
                Q_4.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==5:
                Q_5.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==6:
                Q_6.append(master_cat[counter]['lmass'])
    elif master_cat[counter]['lmass'] < 7.3:
        if master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==1:      #secure cluster, SF
            a+=1
        elif master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==2: #secure cluster, Q
            b+=1
    counter +=1  
#
## Parallel field sample
#BCG_cluster_par = [[0]*2]
#SF_hist_par = Column(np.zeros(shape=(np.sum(phot_mems_par)+np.sum(smem_par),1), dtype='float64'))
#SF_index = 0
#Q_hist_par = Column(np.empty(shape=(np.sum(phot_memq_par)+np.sum(qmem_par),1), dtype='float64'))
#Q_index = 0
#other = 0
#counter = 0
#size = len(master_cat_par)
#while counter < size:
#    if master_cat_par[counter]['member'] == 0 and master_cat_par[counter]['type'] ==1:      #secure cluster, SF
#        SF_hist_par[SF_index][0] = master_cat_par[counter]['lmass']
#        SF_index +=1
#        if master_cat_par[counter]['BCG'] ==2:
#            BCG_cluster_par[0][0] +=1
#    elif master_cat[counter]['member'] == 0 and master_cat_par[counter]['type'] ==2: #secure cluster, Q
#        Q_hist_par[Q_index][0] = master_cat_par[counter]['lmass']
#        Q_index +=1
#        if master_cat_par[counter]['BCG'] ==2:
#            BCG_cluster_par[0][1] +=1
#    else:
#        other +=1
#    counter +=1      
#
## Field sample
BCG_field = [[0]*2]
SF_field_hist = Column(np.zeros(shape=(np.sum(phot_fields)+np.sum(sfield)+np.sum(phot_fields_par)+np.sum(sfield_par),1), dtype='float64'))
SF_index = 0
Q_field_hist = np.empty(shape=(np.sum(phot_fieldq)+np.sum(qfield)+np.sum(phot_fieldq_par)+np.sum(qfield_par),1), dtype='float64')
Q_index = 0
other = 0
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['lmass'] > 7.3:
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
    counter +=1   
#now go through parallel field catalogue to augment the field sample
other = 0
counter = 0
size = len(master_cat_par)
while counter < size:
    if master_cat[counter]['lmass'] > 7.3:
        if master_cat_par[counter]['member'] == 1 and master_cat_par[counter]['type'] ==1:      #secure field, SF
            SF_field_hist[SF_index][0] = master_cat_par[counter]['lmass']
            SF_index +=1
            if master_cat_par[counter]['BCG'] ==2:
                BCG_field[0][0] +=1
        elif master_cat_par[counter]['member'] == 1 and master_cat_par[counter]['type'] ==2: #secure field, Q
            Q_field_hist[Q_index][0] = master_cat_par[counter]['lmass']
            Q_index +=1
            if master_cat_par[counter]['BCG'] ==2:
                BCG_field[0][1] +=1
    counter +=1     
#
## to remove entries which return "NAN" (not a number). no idea why this would 
## happen, but this is just a patch so the program runs. have set all such 
## values equal to 0, and will be removed from histogram as they are out of range
a = 0   #checks to see which lists have NANs in them
b = 0
counter=0           
size=len(SF_hist)
while counter<size:
    if np.isnan(SF_hist[counter]):
        SF_hist[counter] = 0
        a+=1
    counter+=1
# same for Q
counter=0           
size=len(Q_hist)
while counter<size:
    if np.isnan(Q_hist[counter]):
        Q_hist[counter] = 0
        b +=1
    counter+=1
# parallel field
#counter=0           
#size=len(SF_hist_par)
#while counter<size:
#    if np.isnan(SF_hist_par[counter]):
#        SF_hist_par[counter] = 0
#    counter+=1
# same for Q
#counter=0           
#size=len(Q_hist_par)
#while counter<size:
#    if np.isnan(Q_hist_par[counter]):
#        Q_hist_par[counter] = 0
#    counter+=1
# now same for field sample
#
c = 0
d = 0
counter=0           
size=len(SF_field_hist)
while counter<size:
    if np.isnan(SF_field_hist[counter]):
        SF_field_hist[counter] = 0
        c+=1
    counter+=1
# same for Q
counter=0           
size=len(Q_field_hist)
while counter<size:
    if np.isnan(Q_field_hist[counter]):
        Q_field_hist[counter] = 0   
        d+=1
    counter+=1
#
#
## SECTION (ii): sort objects into histogram bins for both SF & Q populations, then sum 
## for 'total' population. use for total pop plot, and for relative fractions
#
## cluster populations arrays: (SF/Q/total)_smf & (SF/Q/total)_field_smf
#
range = [7.3,12.15]     #sets range for all histrograms to be computer: cluster,field,false pos/neg
num_bins = 14
SF_smf, SF_bins = np.histogram(SF_hist, bins=num_bins,range=range)
Q_smf, Q_bins = np.histogram(Q_hist, bins=num_bins)      #use same bins for both populations to facilitate comparison

## smf histograms for individual clusters
#
SF_smf1, SF_bins = np.histogram(SF_1, bins=num_bins,range=range)
SF_smf2, SF_bins = np.histogram(SF_2, bins=num_bins,range=range)
SF_smf3, SF_bins = np.histogram(SF_3, bins=num_bins,range=range)
SF_smf4, SF_bins = np.histogram(SF_4, bins=num_bins,range=range)
SF_smf5, SF_bins = np.histogram(SF_5, bins=num_bins,range=range)
SF_smf6, SF_bins = np.histogram(SF_6, bins=num_bins,range=range)
Q_smf1, SF_bins = np.histogram(Q_1, bins=num_bins,range=range)
Q_smf2, SF_bins = np.histogram(Q_2, bins=num_bins,range=range)
Q_smf3, SF_bins = np.histogram(Q_3, bins=num_bins,range=range)
Q_smf4, SF_bins = np.histogram(Q_4, bins=num_bins,range=range)
Q_smf5, SF_bins = np.histogram(Q_5, bins=num_bins,range=range)
Q_smf6, SF_bins = np.histogram(Q_6, bins=num_bins,range=range)

## parallel field
#SF_smf_par, SF_bins_par = np.histogram(SF_hist_par, bins=SF_bins,range=range)
#Q_smf_par, Q_bins_par = np.histogram(Q_hist_par, bins=SF_bins,range=range)      #use same bins for both populations to facilitate comparison
#total_smf_par = np.empty_like(SF_smf_par, dtype='float64')
#frac_smf_par = np.empty(shape=(2,len(SF_bins)-1),dtype='float64')       #1st row: SF fraction; 2nd row Q fraction
#total_smf_par = SF_smf_par + Q_smf_par
## calculate bin fractions
#counter = 0
#size = len(SF_smf_par)#size = frac_smf.shape[1]
#while counter < size:             #loop to skip bins with no entires. will probably need to reduce range for completeness
#    if total_smf_par[counter] == 0:
#        frac_smf_par[0][counter] = 0
#        frac_smf_par[1][counter] = 0
#    else: 
#        frac_smf_par[0][counter] =  SF_smf_par[counter] / total_smf_par[counter]
#        frac_smf_par[1][counter] =  Q_smf_par[counter] / total_smf_par[counter]
#    counter +=1 
#
#
## Section (iii): find midpoints of spec bins, and split hist. into "x" & "y" for scatter plotting
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
#
#
### Mass corrections
### compute total number of SF/Q galaxies in all 6 clusters, then for individual 
### cluster down to its limiting mass. take the ratio and divide the number count 
### for that galaxy in that mass bin by the ratio to correct for mass comleteness 
### by cluster
#
#
limmass = [7.5,7.8,8.0,7.5,7.4,7.3]  # limiting mass by cluster as determined by magnitude-to-mass plot
## cluster 1: macs1149; limiting mass 7.5
gal_count1 = np.array([0,0],dtype = float)  # [SF, Q]
total_count1 = np.array([0,0],dtype = float)
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==1: #secure member, SF
        if master_cat[counter]['lmass'] > limmass[0]:
            total_count1[0] +=1
            if master_cat[counter]['cluster'] == 2:
                gal_count1[0] +=1
    elif master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==2: #secure member, SF
        if master_cat[counter]['lmass'] > limmass[0]:
            total_count1[1] +=1
            if master_cat[counter]['cluster'] == 2:
                gal_count1[1] +=1
    counter +=1
mass_correction1 = gal_count1/total_count1
SF_smf1[0] = SF_smf1[0]/mass_correction1[0]
Q_smf1[0] = Q_smf1[0]/mass_correction1[1]
#
## cluster 2: macs1149; limiting mass 7.8
gal_count2 = np.array([0,0],dtype = float)  # [SF, Q]
total_count2 = np.array([0,0],dtype = float)
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==1: #secure member, SF
        if master_cat[counter]['lmass'] > limmass[1]:
            total_count2[0] +=1
            if master_cat[counter]['cluster'] == 2:
                gal_count2[0] +=1
    elif master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==2: #secure member, SF
        if master_cat[counter]['lmass'] > limmass[1]:
            total_count2[1] +=1
            if master_cat[counter]['cluster'] == 2:
                gal_count2[1] +=1
    counter +=1
mass_correction2 = gal_count2/total_count2
SF_smf2[0] = SF_smf2[0]/mass_correction2[0]
SF_smf2[1] = SF_smf2[0]/mass_correction2[0]
Q_smf2[0] = Q_smf2[1]/mass_correction2[1]
Q_smf2[1] = Q_smf2[1]/mass_correction2[1]
#
## cluster 3: macs 0717; limiting mass 8.0
gal_count3 = np.array([0,0],dtype = float)  # [SF, Q]
total_count3 = np.array([0,0],dtype = float)
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==1: #secure member, SF
        if master_cat[counter]['lmass'] > limmass[2]:
            total_count3[0] +=1
            if master_cat[counter]['cluster'] == 3:
                gal_count3[0] +=1
    elif master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==2: #secure member, Q
        if master_cat[counter]['lmass'] > limmass[2]:
            total_count3[1] +=1
            if master_cat[counter]['cluster'] == 3:
                gal_count3[1] +=1
    counter +=1
mass_correction3 = gal_count3/total_count3
SF_smf3[0] = SF_smf3[0]/mass_correction3[0]
SF_smf3[1] = SF_smf3[1]/mass_correction3[1]
Q_smf3[0] = Q_smf3[1]/mass_correction3[1]
Q_smf3[1] = Q_smf3[1]/mass_correction3[1]
#
## cluster 4: abell370; limiting mass 7.5
gal_count4 = np.array([0,0],dtype = float)  # [SF, Q]
total_count4 = np.array([0,0],dtype = float)
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==1: #secure member, SF
        if master_cat[counter]['lmass'] > limmass[3]:
            total_count4[0] +=1
            if master_cat[counter]['cluster'] == 4:
                gal_count4[0] +=1
    elif master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==2: #secure member, SF
        if master_cat[counter]['lmass'] > limmass[3]:
            total_count4[1] +=1
            if master_cat[counter]['cluster'] == 4:
                gal_count4[1] +=1
    counter +=1
mass_correction4 = gal_count4/total_count4
SF_smf4[0] = SF_smf4[0]/mass_correction4[0]
Q_smf4[0] = Q_smf4[1]/mass_correction4[1]
#
## cluster 5: abell1063; limiting mass 7.4
gal_count5 = np.array([0,0],dtype = float)  # [SF, Q]
total_count5 = np.array([0,0],dtype = float)
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==1: #secure member, SF
        if master_cat[counter]['lmass'] > limmass[4]:
            total_count5[0] +=1
            if master_cat[counter]['cluster'] == 5:
                gal_count5[0] +=1
    elif master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==2: #secure member, SF
        if master_cat[counter]['lmass'] > limmass[4]:
            total_count5[1] +=1
            if master_cat[counter]['cluster'] == 5:
                gal_count5[1] +=1
    counter +=1
mass_correction5 = gal_count5/total_count5
SF_smf5[0] = SF_smf5[0]/mass_correction5[0]
Q_smf5[0] = Q_smf5[1]/mass_correction5[1]
#
## cluster 6: abell2744; limiting mass 7.3
gal_count6 = np.array([0,0],dtype = float)  # [SF, Q]
total_count6 = np.array([0,0],dtype = float)
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==1: #secure member, SF
        if master_cat[counter]['lmass'] > limmass[5]:
            total_count6[0] +=1
            if master_cat[counter]['cluster'] == 6:
                gal_count6[0] +=1
    elif master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==2: #secure member, SF
        if master_cat[counter]['lmass'] > limmass[5]:
            total_count6[1] +=1
            if master_cat[counter]['cluster'] == 6:
                gal_count6[1] +=1
    counter +=1
mass_correction6 = gal_count6/total_count6
SF_smf6[0] = SF_smf6[0]/mass_correction6[0]
Q_smf6[0] = Q_smf6[1]/mass_correction6[1]
#
##  Section (v): caclulate flase pos/neg fractions by cluster for completeness 
##          corrections in photometric sample
## Fig. 2: mass histrograms by cluster for false positive negative to determine 
##          mass bins for completeness corrections. fraction > 1 means more false
##          positives than negative (i.e. more galaxies erroneously classified as
##          cluster members than members classified as field), and would require a
##          downward correction of the number count in that bin
#
SF_pos1 = []
Q_pos1 = []
SF_neg1 = []
Q_neg1 = []
SF_pos2 = []
Q_pos2 = []
SF_neg2 = []
Q_neg2 = []
SF_pos3 = []
Q_pos3 = []
SF_neg3 = []
Q_neg3 = []
SF_pos4 = []
Q_pos4 = []
SF_neg4 = []
Q_neg4 = []
SF_pos5 = []
Q_pos5 = []
SF_neg5 = []
Q_neg5 = []
SF_pos6 = []
Q_pos6 = []
SF_neg6 = []
Q_neg6 = []
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['member'] == 2:      #false pos
        if master_cat[counter]['type'] == 1:    #SF
            if master_cat[counter]['cluster'] ==1:
                SF_pos1.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==2:
                SF_pos2.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==3:
                SF_pos3.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==4:
                SF_pos4.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==5:
                SF_pos5.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==6:
                SF_pos6.append(master_cat[counter]['lmass'])
        elif master_cat[counter]['type'] == 2:    #Quiescent
            if master_cat[counter]['cluster'] ==1:
                Q_pos1.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==2:
                Q_pos2.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==3:
                Q_pos3.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==4:
                Q_pos4.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==5:
                Q_pos5.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==6:
                Q_pos6.append(master_cat[counter]['lmass'])
    elif master_cat[counter]['member'] == 3:      #false neg
        if master_cat[counter]['type'] == 1:    #SF
            if master_cat[counter]['cluster'] ==1:
                SF_neg1.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==2:
                SF_neg2.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==3:
                SF_neg3.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==4:
                SF_neg4.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==5:
                SF_neg5.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==6:
                SF_neg6.append(master_cat[counter]['lmass'])
        elif master_cat[counter]['type'] == 2:    #Quiescent
            if master_cat[counter]['cluster'] ==1:
                Q_neg1.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==2:
                Q_neg2.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==3:
                Q_neg3.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==4:
                Q_neg4.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==5:
                Q_neg5.append(master_cat[counter]['lmass'])
            elif master_cat[counter]['cluster'] ==6:
                Q_neg6.append(master_cat[counter]['lmass'])
    counter +=1
#
## bin SF & Q for all 6 clusters, then compute false pos/neg fractions by bin for correction factors
num_binsSF = [4,4,4,4,4,4]#   [5,5,5,5,5,5]#
num_binsQ = [6,6,6,6,6,6]#[7,7,7,7,7,7]#
SF_pos_hist1, binsSF1 = np.histogram(SF_pos1, bins=num_binsSF[0],range=range)
SF_neg_hist1, binsSF1 = np.histogram(SF_neg1, bins=num_binsSF[0],range=range)
Q_pos_hist1, binsQ1 = np.histogram(Q_pos1, bins=num_binsQ[0],range=range)
Q_neg_hist1, binsQ1 = np.histogram(Q_neg1, bins=num_binsQ[0],range=range)
SF_pos_hist2, binsSF2 = np.histogram(SF_pos2, bins=num_binsSF[1],range=range)
SF_neg_hist2, binsSF2 = np.histogram(SF_neg2, bins=num_binsSF[1],range=range)
Q_pos_hist2, binsQ2 = np.histogram(Q_pos2, bins=num_binsQ[1],range=range)
Q_neg_hist2, binsQ2 = np.histogram(Q_neg2, bins=num_binsQ[1],range=range)
SF_pos_hist3, binsSF3 = np.histogram(SF_pos3, bins=num_binsSF[2],range=range)
SF_neg_hist3, binsSF3 = np.histogram(SF_neg3, bins=num_binsSF[2],range=range)
Q_pos_hist3, binsQ3 = np.histogram(Q_pos3, bins=num_binsQ[2],range=range)
Q_neg_hist3, binsQ3 = np.histogram(Q_neg3, bins=num_binsQ[2],range=range)
SF_pos_hist4, binsSF4 = np.histogram(SF_pos4, bins=num_binsSF[3],range=range)
SF_neg_hist4, binsSF4 = np.histogram(SF_neg4, bins=num_binsSF[3],range=range)
Q_pos_hist4, binsQ4 = np.histogram(Q_pos4, bins=num_binsQ[3],range=range)
Q_neg_hist4, binsQ4 = np.histogram(Q_neg4, bins=num_binsQ[3],range=range)
SF_pos_hist5, binsSF5 = np.histogram(SF_pos5, bins=num_binsSF[4],range=range)
SF_neg_hist5, binsSF5 = np.histogram(SF_neg5, bins=num_binsSF[4],range=range)
Q_pos_hist5, binsQ5 = np.histogram(Q_pos5, bins=num_binsQ[4],range=range)
Q_neg_hist5, binsQ5 = np.histogram(Q_neg5, bins=num_binsQ[4],range=range)
SF_pos_hist6, binsSF6 = np.histogram(SF_pos6, bins=num_binsSF[5],range=range)
SF_neg_hist6, binsSF6 = np.histogram(SF_neg6, bins=num_binsSF[5],range=range)
Q_pos_hist6, binsQ6 = np.histogram(Q_pos6, bins=num_binsQ[5],range=range)
Q_neg_hist6, binsQ6 = np.histogram(Q_neg6, bins=num_binsQ[5],range=range)
#
## compute false pos/ false neg ratio by cluster
SF_pos_hist = SF_pos_hist1+SF_pos_hist2+SF_pos_hist3+SF_pos_hist4+SF_pos_hist5+SF_pos_hist6
SF_neg_hist = SF_neg_hist1+SF_neg_hist2+SF_neg_hist3+SF_neg_hist4+SF_neg_hist5+SF_neg_hist6
Q_pos_hist = Q_pos_hist1+Q_pos_hist2+Q_pos_hist3+Q_pos_hist4+Q_pos_hist5+Q_pos_hist6
Q_neg_hist = Q_neg_hist1+Q_neg_hist2+Q_neg_hist3+Q_neg_hist4+Q_neg_hist5+Q_neg_hist6
SF_frac = np.empty_like(SF_pos_hist, dtype='float64')
Q_frac = np.empty_like(Q_pos_hist, dtype='float64')
#
#SF_frac1 = np.empty_like(SF_pos_hist1, dtype='float64')
#Q_frac1 = np.empty_like(Q_pos_hist1, dtype='float64')
#SF_frac2 = np.empty_like(SF_pos_hist2, dtype='float64')
#Q_frac2 = np.empty_like(Q_pos_hist2, dtype='float64')
#SF_frac3 = np.empty_like(SF_pos_hist3, dtype='float64')
#Q_frac3 = np.empty_like(Q_pos_hist3, dtype='float64')
#SF_frac4 = np.empty_like(SF_pos_hist4, dtype='float64')
#Q_frac4 = np.empty_like(Q_pos_hist4, dtype='float64')
#SF_frac5 = np.empty_like(SF_pos_hist5, dtype='float64')
#Q_frac5 = np.empty_like(Q_pos_hist5, dtype='float64')
#SF_frac6 = np.empty_like(SF_pos_hist6, dtype='float64')
#Q_frac6 = np.empty_like(Q_pos_hist6, dtype='float64')
#compute fractions
# total sample
zero_falseSF = 0
counter = 0         #SF pop
size = len(SF_pos_hist)
while counter < size:
    if SF_pos_hist[counter] ==0 and SF_neg_hist[counter] ==0:
        SF_frac[counter] = 1
        zero_falseSF +=1
    elif SF_pos_hist[counter] ==0:
        SF_frac[counter] = 10*SF_neg_hist[counter]
    elif SF_neg_hist[counter] ==0:
        SF_frac[counter] = 100*SF_pos_hist[counter]
    else: SF_frac[counter] = SF_pos_hist[counter]/SF_neg_hist[counter]
    counter+=1
zero_falseQ = 0
counter = 0         #Q pop
size = len(Q_pos_hist)
while counter < size:
    if Q_pos_hist[counter] ==0 and Q_neg_hist[counter] ==0:
        Q_frac[counter] = 1
        zero_falseQ +=1
    elif Q_pos_hist[counter] ==0:
        Q_frac[counter] = 10*Q_neg_hist[counter]
    elif Q_neg_hist[counter] ==0:
        Q_frac[counter] = 100*Q_pos_hist[counter]
    else: Q_frac[counter] = Q_pos_hist[counter]/Q_neg_hist[counter]
    counter+=1
#
#
#
# by cluster
# 1.
#counter = 0
#size = len(SF_pos_hist1)
#while counter < size:
#    if SF_pos_hist1[counter] ==0:
#        SF_frac1[counter] = 10*SF_neg_hist1[counter]
#    elif SF_neg_hist1[counter] ==0:
#        SF_frac1[counter] = 100*SF_pos_hist1[counter]
#    else: SF_frac1[counter] = SF_pos_hist1[counter]/SF_neg_hist1[counter]
#    counter+=1
#counter = 0         #Q pop
#size = len(Q_pos_hist1)
#while counter < size:
#    if Q_pos_hist1[counter] ==0:
#        Q_frac1[counter] = 10*Q_neg_hist1[counter]
#    elif Q_neg_hist1[counter] ==0:
#        Q_frac1[counter] = 100*Q_pos_hist1[counter]
#    else: Q_frac1[counter] = Q_pos_hist1[counter]/Q_neg_hist1[counter]
#    counter+=1
# 2.
#counter = 0
#size = len(SF_pos_hist2)
#while counter < size:
#    if SF_pos_hist2[counter] ==0:
#        SF_frac2[counter] = 10*SF_neg_hist2[counter]
#    elif SF_neg_hist2[counter] ==0:
#        SF_frac2[counter] = 100*SF_pos_hist2[counter]
#    else: SF_frac2[counter] = SF_pos_hist2[counter]/SF_neg_hist2[counter]
#    counter+=1
#counter = 0         #Q pop
#size = len(Q_pos_hist2)
#while counter < size:
#    if Q_pos_hist2[counter] ==0:
#        Q_frac2[counter] = 10*Q_neg_hist2[counter]
#    elif Q_neg_hist2[counter] ==0:
#        Q_frac2[counter] = 100*Q_pos_hist2[counter]
#    else: Q_frac2[counter] = Q_pos_hist2[counter]/Q_neg_hist2[counter]
#    counter+=1
#
#
#
# compute midbins
SF_frac_midbins = np.empty_like(SF_frac, dtype='float64')
Q_frac_midbins = np.empty_like(Q_frac, dtype='float64')
#
size = len(binsSF1)-1
counter = 0
while counter < size:
    SF_frac_midbins[counter] = (binsSF1[counter] + binsSF1[counter+1])/2
    counter +=1
size = len(binsQ1)-1
counter = 0
while counter < size:
    Q_frac_midbins[counter] = (binsQ1[counter] + binsQ1[counter+1])/2
    counter +=1
#
## table to summarize binning by cluster
#cluster_names = Column(['1. macs0416','2. macs1149','3. macs0717','4. abell370','5. abell1063std','6. abell2744'])
#SFpos = Column([SF_pos1,SF_pos2,SF_pos3,SF_pos4,SF_pos5,SF_pos6])
#SFneg = Column([SF_neg1,SF_neg2,SF_neg3,SF_neg4,SF_neg5,SF_neg6])
#Qpos = Column([Q_pos1,Q_pos2,Q_pos3,Q_pos4,Q_pos5,Q_pos6])
#Qneg = Column([Q_neg1,Q_neg2,Q_neg3,Q_neg4,Q_neg5,Q_neg6])
#binningSF = Column([num_binsSF[0],num_binsSF[1],num_binsSF[2],num_binsSF[3],num_binsSF[4],num_binsSF[5]])
#binningQ = Column([num_binsQ[0],num_binsQ[1],num_binsQ[2],num_binsQ[3],num_binsQ[4],num_binsQ[5]])
#cluster_bins1 = Table([cluster_names,SFpos,SFneg,binningSF],names=('Cluster','SF pos','SF neg','# SF bins'))
#cluster_bins2 = Table([cluster_names,Qpos,Qneg,binningQ],names=('Cluster','Q pos','Q neg','# Q bins'))
##
#errors
SF_perr = (np.sqrt(SF_pos_hist))/SF_pos_hist    #compute rel. error for false pos/neg
SF_nerr = (np.sqrt(SF_neg_hist))/SF_neg_hist
Q_perr = (np.sqrt(Q_pos_hist))/Q_pos_hist
Q_nerr = (np.sqrt(Q_neg_hist))/Q_neg_hist
SF_frac_error = np.sqrt((SF_perr**2)+(SF_nerr**2))      #computer error of correction factor (sum in quadrature of rel. errors in false pos/neg)
Q_frac_error = np.sqrt((Q_perr**2)+(Q_nerr**2))
#
# sift out 'NAN' from fractions for 0/0
counter=0                   #SF
size=len(SF_frac_error)
while counter<size:
    if np.isnan(SF_frac_error[counter]):
        SF_frac_error[counter] = 0
    counter+=1
counter=0                   #Q
size=len(Q_frac_error)
while counter<size:
    if np.isnan(Q_frac_error[counter]):
        Q_frac_error[counter] = 0
    counter+=1
#
# plot Spectroscopic completion correction factors 
plt.close()
MC = plt.figure(num=2)
MC.suptitle('Mass Completeness Correction Factors')
#plt.plot(SF_midbins,)
plt.errorbar(SF_frac_midbins,SF_frac,yerr=SF_frac_error, fmt='ob',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.8, mfc='none')
plt.errorbar(Q_frac_midbins,Q_frac,yerr=Q_frac_error, fmt='or',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.8, mfc='none')
plt.plot(SF_frac_midbins,SF_frac,'--b', linewidth=0.5, label='Star-forming')
plt.plot(Q_frac_midbins,Q_frac,'--r', linewidth=0.5, label='Quiescent')
plt.plot([0,13],[1,1],'--k',linewidth = 0.5)
plt.legend(loc='upper right', frameon=False)
plt.xlim=(7.5,12)
plt.xlabel('$log(M/M_{\odot})$')
plt.ylim=(-0.5,4.1)
plt.ylabel('Correction factor\n(false pos / false neg)')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on', labelleft='on')
plt.minorticks_on()
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
#
### NOTE: the following is just to interpolate/extrapolate where the bins from 
###       SMF lie within the bin correction range
#add veritcal lines to determine correciton factors for actual bins used in SMF
plt.plot([SF_midbins[0],SF_midbins[0]],[-0.5,6],'--g',linewidth = 0.9)
plt.plot([SF_midbins[1],SF_midbins[1]],[-0.5,6],'--g',linewidth = 0.9)
plt.plot([SF_midbins[2],SF_midbins[2]],[-0.5,6],'--g',linewidth = 0.9)
plt.plot([SF_midbins[3],SF_midbins[3]],[-0.5,6],'--g',linewidth = 0.9)
plt.plot([SF_midbins[4],SF_midbins[4]],[-0.5,6],'--g',linewidth = 0.9)
plt.plot([SF_midbins[5],SF_midbins[5]],[-0.5,6],'--g',linewidth = 0.9)
plt.plot([SF_midbins[6],SF_midbins[6]],[-0.5,6],'--g',linewidth = 0.9)
plt.plot([SF_midbins[7],SF_midbins[7]],[-0.5,6],'--g',linewidth = 0.9)
plt.plot([SF_midbins[8],SF_midbins[8]],[-0.5,6],'--g',linewidth = 0.9)
plt.plot([SF_midbins[9],SF_midbins[9]],[-0.5,6],'--g',linewidth = 0.9)
plt.plot([SF_midbins[10],SF_midbins[10]],[-0.5,6],'--g',linewidth = 0.9)
plt.plot([SF_midbins[11],SF_midbins[11]],[-0.5,6],'--g',linewidth = 0.9)
plt.plot([SF_midbins[12],SF_midbins[12]],[-0.5,6],'--g',linewidth = 0.9)
plt.plot([SF_midbins[13],SF_midbins[13]],[-0.5,6],'--g',linewidth = 0.9)
#plt.plot([SF_midbins[14],SF_midbins[14]],[-0.5,6],'--g',linewidth = 0.9)
#plt.plot([SF_midbins[15],SF_midbins[15]],[-0.5,6],'--g',linewidth = 0.9)
## extrapolating SF at lo-mass
m=(SF_frac[1] - SF_frac[0])/(SF_frac_midbins[1] - SF_frac_midbins[0])
b = SF_frac[1] - (SF_frac_midbins[1]*m)
X = np.linspace(7.4,8.1,50)
y=[]
for x in X:
    y.append(m*x + b)
plt.plot(X,y)
## extrapolating Q at lo-mass
m=(Q_frac[1] - Q_frac[0])/(Q_frac_midbins[1] - Q_frac_midbins[0])
b = Q_frac[1] - (Q_frac_midbins[1]*m)
X = np.linspace(7.4,8.1,50)
y=[]
for x in X:
    y.append(m*x + b)
plt.plot(X,y)
## extrapolating SF at hi-mass
m=(SF_frac[3] - SF_frac[2])/(SF_frac_midbins[3] - SF_frac_midbins[2])
b = SF_frac[3] - (SF_frac_midbins[3]*m)
X = np.linspace(11.6,12.1,50)
y=[]
for x in X:
    y.append(m*x + b)
plt.plot(X,y)
## extrapolating Q at hi-mass
m=(Q_frac[5] - Q_frac[4])/(Q_frac_midbins[5] - Q_frac_midbins[4])
b = Q_frac[5] - (Q_frac_midbins[5]*m)
X = np.linspace(11.6,12.1,50)
y=[]
for x in X:
    y.append(m*x + b)
plt.plot(X,y)
#
# list of correction factors, determined by inspection interpolating SMF mass 
# bins between the mass correction bins
SF_correction = [2.54063,2.30573,2.20829,2.04172,1.87286,1.69589,1.47983,1.26787,1.04994,1,1,1,1,1]      #0.03 & 0.06 cutoff
Q_correction = [0.54764,.47627,.40492,.33385,.48863,.64394,.80382,.96587,1.02776,.88903,.75,.71445,.67861,.64285]
#SF_correction = [2.5859,2.51,2.4295,2.353,2.2734,1.9766,1.585937,1.1953125,1.0625,1.1875,1.3125,1.3625,1.2375,1.1125,0.985,.8625]      for cutoff .03 & .07 (SF/Q)
#Q_correction = [1.20723684,.885,.5625,.2865,.6096,.9309,1.25776,1.5859375,1.625,1.375,1.16183375,1.083706,1.00558,.796875,.578125,.359375]
#
#SF_correction = [2.75189,2.69367,2.6152,2.548,2.47931,2.34368,2.00757,1.67038,1.33367,1.19593,1.12456,1.05547,.981921,.910591]
#Q_correction = [.54758,.476075,.404639,.330849,.487969,.643279,.802955,.9656,1.02748,.889034,.750098,.714258,.677714,.643279]
###     NOTE: relative fraction calculated above at end of Seciont (ii), and are
###     stored in arrays "frac_smf" & "frac_field_smf"
#
total_smf = np.empty_like(SF_smf, dtype='float64')
##  apply correction factors to SMFs:
SF_smf = SF_smf/SF_correction
Q_smf = Q_smf/Q_correction
total_smf = SF_smf + Q_smf
#
#
## field populations
SF_field_smf, SF_field_bins = np.histogram(SF_field_hist, bins=SF_bins,range=range)
Q_field_smf, Q_field_bins = np.histogram(Q_field_hist, bins=SF_bins,range=range)      #use same bins for both populations to facilitate comparison
frac_field_smf = np.empty(shape=(2,len(SF_bins)-1),dtype='float64')       #1st row: SF fraction; 2nd row Q fraction
#apply correction factors to field
SF_field_smf = SF_field_smf/SF_correction
Q_field_smf = Q_field_smf/Q_correction
total_field_smf = SF_field_smf + Q_field_smf
#
frac_smf = np.empty(shape=(2,len(SF_bins)-1),dtype='float64')       #1st row: SF fraction; 2nd row Q fraction
## calculate bin fractions
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
## calculate bin fractions for field
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
## I. error values by sample type, cluster & field
SF_error = np.sqrt(SF_smf)        #assign errors by population as sqrt of count in each mass bin
Q_error = np.sqrt(Q_smf)
total_error = np.sqrt(total_smf)
SF_field_error = np.sqrt(SF_field_smf)        #assign errors by population as sqrt of count in each mass bin
Q_field_error = np.sqrt(Q_field_smf)
total_field_error = np.sqrt(total_field_smf)
frac_error = np.empty_like(frac_smf, dtype='float32')
frac_field_error = np.empty_like(frac_field_smf, dtype='float32')
# cluster SF
counter = 0
size = len(SF_smf)
while counter < size:
    if SF_smf[counter] == 0 or SF_error[counter] ==0:
        frac_error [0][counter] = 0
    elif SF_smf[counter] >0 and SF_smf[counter]<1:
        frac_error [0][counter] = 0
    elif SF_smf[counter] != 0 and SF_error[counter] !=0:
        frac_error [0][counter] = (SF_error[counter]/SF_smf[counter])
    counter +=1
#cluster Q
counter = 0
size = len(Q_smf)
while counter < size:
    if Q_smf[counter] == 0 or Q_error[counter] ==0:
        frac_error [1][counter] = 0
    elif Q_smf[counter] >0 and Q_smf[counter]<1:
        frac_error [1][counter] = 0
    elif Q_smf[counter] != 0 and Q_error[counter] !=0:
        frac_error [1][counter] = (Q_error[counter]/Q_smf[counter])
    counter +=1
# field SF
counter = 0
size = len(SF_field_smf)
while counter < size:
    if SF_field_smf[counter] == 0 or SF_field_error[counter] ==0:
        frac_field_error [0][counter] = 0
    elif SF_field_smf[counter] >0 and SF_field_smf[counter]<1:
        frac_field_error [0][counter] = 0
    elif SF_field_smf[counter] != 0 and SF_field_error[counter] !=0:
        frac_field_error [0][counter] = (SF_field_error[counter]/SF_field_smf[counter])
    counter +=1
#field Q
counter = 0
size = len(Q_field_smf)
while counter < size:
    if Q_field_smf[counter] == 0 or Q_field_error[counter] ==0:
        frac_field_error [1][counter] = 0
    elif Q_field_smf[counter] >0 and Q_field_smf[counter]<1:
        frac_error [1][counter] = 0
    elif Q_field_smf[counter] != 0 and Q_field_error[counter] !=0:
        frac_field_error [1][counter] = (Q_field_error[counter]/Q_field_smf[counter])
    counter +=1

#offset Q data points by 0.02 for clarity
Q_midbins = np.empty_like(SF_midbins, dtype='float64')
counter = 0
size = len(SF_midbins)
while counter < size:
    Q_midbins[counter] = SF_midbins[counter] + 0.02
    counter +=1
#

#
## (vi) fit a SCHECHTER function to scatter plot
#
#
#
#
## SECTION (vii): create plot
## upper: SMF for cluster, field
## lower: fractions of SF/Q in cluster, field
#
## SMF
## TOTAL (all clusters)
plt.close()
SMF = plt.figure(num=1)
gs = gridspec.GridSpec(2,2, width_ratios=[1,1], height_ratios=[2,1])   #make a tiled-plot like vdB2013 w/ fractions below, this line sets the proporitons of plots in the figure
#gs = gridspec.GridSpec(2,3, width_ratios=[1,1,1], height_ratios=[2,1])   #make a tiled-plot like vdB2013 w/ fractions below, this line sets the proporitons of plots in the figure
# Cluster
cluster = plt.subplot(gs[0])      
xa = plt.errorbar(SF_midbins,SF_smf,yerr=SF_error, fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)#yerr=SF_error,
xb = plt.errorbar(Q_midbins,Q_smf,yerr=Q_error,fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
xc = plt.errorbar(SF_midbins,total_smf,yerr=total_error,fmt='.k',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
#plt.plot(SF_midbins,Q_smf,'.r',linewidth=0.5)#Qx_new, Qy_new,'-r',linewidth=0.5)
#plt.plot(SF_midbins,total_smf,'.k',linewidth=0.5)#totalx_new, totaly_new,'-k',linewidth=0.5)
#plt.errorbar(SF_midbins,SF_smf,fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)#yerr=SF_error,
plt.xscale('linear')
#plt.xlabel('log(M/M_sol)')
plt.xlim=(8,12.25)
plt.yscale('log')
cluster.minorticks_on()
cluster.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off',labelbottom='off')
cluster.yaxis.set_label_position("left")
plt.ylabel('# count')
plt.ylim=(1,750)
plt.title('Cluster SMF')
plt.legend((xa,xb,xc),('Star-forming','Quiescent','Total'),scatterpoints=1,loc='lower left', frameon=False, fontsize = 'x-small')
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')
#
# Parallel field
#cluster = plt.subplot(gs[1])      
#plt.plot(SF_midbins,Q_smf_par,'.r',linewidth=0.5)#Qx_new, Qy_new,'-r',linewidth=0.5)
#plt.plot(SF_midbins,total_smf_par,'.k',linewidth=0.5)#totalx_new, totaly_new,'-k',linewidth=0.5)
#plt.errorbar(SF_midbins,SF_smf_par,fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)#yerr=SF_error,
#plt.xscale('linear')
#plt.xlabel('log(M/M_sol)')
#plt.xlim(7,12.5)
#plt.yscale('log')
#cluster.minorticks_on()
#cluster.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off',labelleft='off',labelbottom='off')
#cluster.yaxis.set_label_position("left")
#plt.ylabel('logarthmic # count')
#plt.ylim(0,750)
#plt.title('Parallel SMF')
#plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')
#
## Field
field = plt.subplot(gs[1])      #make a tiled-plot like vdB2013 w/ fractions below
plt.errorbar(SF_midbins,SF_field_smf,yerr=SF_field_error, fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
plt.errorbar(Q_midbins,Q_field_smf,yerr=Q_field_error, fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
plt.errorbar(SF_midbins,total_field_smf,yerr=total_field_error, fmt='.k',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
plt.xscale('linear')
#plt.xlabel('log(M/M_sol)')
plt.xlim=(8,12.25)
plt.yscale('log')
plt.ylabel('# count')
plt.ylim=(1,750)
field.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on',labelleft='off',labelbottom='off')
field.minorticks_on()
field.yaxis.set_label_position("right")
plt.title('Field SMF')
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')
#plt.text()      #this could replace the legend to avoid visual clutter
#
## cluster fraction
fr_cl = plt.subplot(gs[2])    
#plt.plot(SF_midbins,frac_smf[0],'.b',linewidth=0.5)
#plt.plot(SF_midbins,frac_smf[1],'.r',linewidth=0.5)
plt.errorbar(SF_midbins,frac_smf[0],yerr=frac_error[0], fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
plt.errorbar(SF_midbins,frac_smf[1],yerr=frac_error[1], fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
plt.xscale('linear')
plt.xlabel('$log(M/M_{\odot})$')
plt.xlim=(8,12.25)
plt.yscale('linear')
plt.ylim=(-0.3,1.3)
fr_cl.minorticks_on()
fr_cl.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off')
fr_cl.yaxis.set_label_position("left")
plt.ylabel('Galaxy type fraction')
#plt.ylim(1,190)
#plt.title('relative fraction SF & Q')
#
## parallel fraction
#fr_cl = plt.subplot(gs[4])    
#plt.plot(SF_midbins,frac_smf_par[0],'.b',linewidth=0.5)
#plt.plot(SF_midbins,frac_smf_par[1],'.r',linewidth=0.5)
#plt.errorbar(SF_midbins,frac_smf[0],yerr=frac_error[0], fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
#plt.errorbar(SF_midbins,frac_smf[1],yerr=frac_error[1], fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
#plt.xscale('linear')
#plt.xlabel('$log(M/M_{\odot})$')
#plt.xlim(7,12.5)
#plt.yscale('linear')
#plt.ylim(-0.1,1.1)
#fr_cl.minorticks_on()
#fr_cl.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off',labelleft='off')
#fr_cl.yaxis.set_label_position("left")
#plt.ylabel('Galaxy type fraction')
#plt.ylim(1,190)
#
## field fraction
fr_fld = plt.subplot(gs[3])     
plt.errorbar(SF_midbins,frac_field_smf[0],yerr=frac_field_error[0], fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)  #,yerr=frac_field_error[0] insert after x,y at beginning of arguement
plt.errorbar(SF_midbins,frac_field_smf[1],yerr=frac_field_error[1], fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5) #,yerr=frac_field_error[1]
plt.xscale('linear')
plt.xlabel('$log(M/M_{\odot})$')
plt.xlim=(8,12.25)
plt.yscale('linear')
plt.ylim=(-0.3,1.3)
fr_fld.minorticks_on()
fr_fld.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on',labelleft='off')
fr_fld.yaxis.set_label_position("right")
plt.ylabel('Galaxy type fraction')
#
plt.show()
#
#
## Individually by cluster
#
plt.close()
cl_SMF = plt.figure(num=2)
plt.suptitle('SMF by cluster')
# 1. macs0416
smf1 = cl_SMF.add_subplot(231)
plt.scatter(SF_midbins, SF_smf1,c='b',marker='*',linewidths=0.5)
plt.scatter(SF_midbins, Q_smf1,c='r',marker='.',linewidths=0.5)
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off',labelbottom='off')
plt.minorticks_on()
plt.text(9,70,'macs0416',fontsize=9)
plt.text(8,1.5,'macs0416',fontsize=9)
plt.xscale('linear')
plt.xlim(7.5,12)
plt.ylabel("count")
plt.yscale('log')
plt.ylim(1,150)
# 2. macs1149
smf2 = cl_SMF.add_subplot(232)
plt.scatter(SF_midbins, SF_smf2,c='b',marker='*',linewidths=0.5)
plt.scatter(SF_midbins, Q_smf2,c='r',marker='.',linewidths=0.5)
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off',labelleft='off',labelbottom='off')
plt.minorticks_on()
plt.text(8,1.5,'macs1149',fontsize=9)
plt.xscale('linear')
plt.xlim(7.5,12)
plt.yscale('log')
plt.ylim(1,150)
# 3. macs0717
smf3 = cl_SMF.add_subplot(233)
plt.scatter(SF_midbins, SF_smf3,c='b',marker='*',linewidths=0.5)
plt.scatter(SF_midbins, Q_smf3,c='r',marker='.',linewidths=0.5)
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on',labelleft='off',labelbottom='off')
plt.minorticks_on()
plt.text(8,1.5,'macs0717',fontsize=9)
plt.xscale('linear')
plt.xlim(7.5,12)
plt.ylabel("count")
smf3.yaxis.set_label_position("right")
plt.yscale('log')
plt.ylim(1,150)
# 4. abell370
smf4 = cl_SMF.add_subplot(234)
plt.scatter(SF_midbins, SF_smf4,c='b',marker='*',linewidths=0.5)
plt.scatter(SF_midbins, Q_smf4,c='r',marker='.',linewidths=0.5)
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off',labelbottom='on')
plt.minorticks_on()
plt.text(8,1.25,'abell370',fontsize=9)
plt.xscale('linear')
plt.xlabel('$log(M/M_{\odot})$')
plt.xlim(7.5,12)
plt.ylabel("count")
plt.yscale('log')
plt.ylim(1,150)
# 5. abell1063
smf5 = cl_SMF.add_subplot(235)
plt.scatter(SF_midbins, SF_smf5,c='b',marker='*',linewidths=0.5)
plt.scatter(SF_midbins, Q_smf5,c='r',marker='.',linewidths=0.5)
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off',labelbottom='on')
plt.minorticks_on()
plt.text(8,1.5,'abell1063',fontsize=9)
plt.xscale('linear')
plt.xlabel('$log(M/M_{\odot})$')
plt.xlim(7.5,12)
plt.yscale('log')
plt.ylim(1,150)
# 6. abell2744
smf6 = cl_SMF.add_subplot(236)
plt.scatter(SF_midbins, SF_smf6,c='b',marker='*',linewidths=0.5)
plt.scatter(SF_midbins, Q_smf6,c='r',marker='.',linewidths=0.5)
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on',labelleft='off',labelbottom='on')
plt.minorticks_on()
plt.text(8,1.75,'abell2274',fontsize=9)
plt.xscale('linear')
plt.xlabel('$log(M/M_{\odot})$')
plt.xlim(7.5,12)
plt.ylabel("count")
smf6.yaxis.set_label_position("right")
plt.yscale('log')
plt.ylim(1,150)
#
## create summary table of binning info - "bin_stats"
bin_stats = Table([SF_midbins,total_smf,SF_smf,Q_smf,total_field_smf,SF_field_smf,Q_field_smf],names=('Bin mass','Total cluster','SF cluster','Q cluster','Total field','SF field','Q field'))
bin_stats_par = Table([SF_midbins,total_smf_par,SF_smf_par,Q_smf_par],names=('Bin mass','Total || field','SF || field','Q || field'))
#
#
