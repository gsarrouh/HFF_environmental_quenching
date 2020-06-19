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
SF_1 = []       #lists of SF/Q by cluster, for flase pos/neg mass completeness corection
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
    if master_cat[counter]['lmass'] > 7.4:
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
    elif master_cat[counter]['lmass'] < 7.4:
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
    if master_cat[counter]['lmass'] > 7.4:
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
    if master_cat[counter]['lmass'] > 7.4:
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
range = [7.5,12.15]     #sets range for all histrograms to be computer: cluster,field,false pos/neg
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
#create histograms & summarize results by mass bin for false pos/neg
    range = [7.5,12.15]
    num_binsSF = [5,2,2,5,5,5]#[4,4,4,4,4,4]#   
    num_binsQ = [3,4,4,3,3,3]#[6,6,6,6,6,6] #[7,7,7,7,7,7]#
    #
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
    #group clusters: group1: 1&4&5&6, group 2: 2&3 by redshift. see if mass correction bins work out better
    SF_pos_group1 = SF_pos_hist1+SF_pos_hist4+SF_pos_hist5+SF_pos_hist6
    SF_neg_group1 = SF_neg_hist1+SF_neg_hist4+SF_neg_hist5+SF_neg_hist6
    Q_pos_group1 = Q_pos_hist1+Q_pos_hist4+Q_pos_hist5+Q_pos_hist6
    Q_neg_group1 = Q_neg_hist1+Q_neg_hist4+Q_neg_hist5+Q_neg_hist6
    SF_pos_group2 = SF_pos_hist2+SF_pos_hist3
    SF_neg_group2 = SF_neg_hist2+SF_neg_hist3
    Q_pos_group2 = Q_pos_hist2+Q_pos_hist3
    Q_neg_group2 = Q_neg_hist2+Q_neg_hist3
    
    ##create summary tables, one for each cluster
    #
    #midbins for each cluster, SF & Q separately as diff # of bins
    ##SF
    #1
    SF_midbins1 = np.empty_like(SF_pos_hist1, dtype='float64')
    size = len(binsSF1)-1
    counter = 0
    while counter < size:
        SF_midbins1[counter] = (binsSF1[counter] + binsSF1[counter+1])/2
        counter +=1
    #2
    SF_midbins2 = np.empty_like(SF_pos_hist2, dtype='float64')
    size = len(binsSF2)-1
    counter = 0
    while counter < size:
        SF_midbins2[counter] = (binsSF2[counter] + binsSF2[counter+1])/2
        counter +=1
    #3
    SF_midbins3 = np.empty_like(SF_pos_hist3, dtype='float64')
    size = len(binsSF3)-1
    counter = 0
    while counter < size:
        SF_midbins3[counter] = (binsSF3[counter] + binsSF3[counter+1])/2
        counter +=1
    #4
    SF_midbins4 = np.empty_like(SF_pos_hist4, dtype='float64')
    size = len(binsSF4)-1
    counter = 0
    while counter < size:
        SF_midbins4[counter] = (binsSF4[counter] + binsSF4[counter+1])/2
        counter +=1
    #5
    SF_midbins5 = np.empty_like(SF_pos_hist5, dtype='float64')
    size = len(binsSF5)-1
    counter = 0
    while counter < size:
        SF_midbins5[counter] = (binsSF5[counter] + binsSF5[counter+1])/2
        counter +=1
    #6
    SF_midbins6 = np.empty_like(SF_pos_hist6, dtype='float64')
    size = len(binsSF6)-1
    counter = 0
    while counter < size:
        SF_midbins6[counter] = (binsSF6[counter] + binsSF6[counter+1])/2
        counter +=1
    ##Q
    #1
    Q_midbins1 = np.empty_like(Q_pos_hist1, dtype='float64')
    size = len(binsQ1)-1
    counter = 0
    while counter < size:
        Q_midbins1[counter] = (binsQ1[counter] + binsQ1[counter+1])/2
        counter +=1
    #2
    Q_midbins2 = np.empty_like(Q_pos_hist2, dtype='float64')
    size = len(binsQ2)-1
    counter = 0
    while counter < size:
        Q_midbins2[counter] = (binsQ2[counter] + binsQ2[counter+1])/2
        counter +=1
    #3
    Q_midbins3 = np.empty_like(Q_pos_hist3, dtype='float64')
    size = len(binsQ3)-1
    counter = 0
    while counter < size:
        Q_midbins3[counter] = (binsQ3[counter] + binsQ3[counter+1])/2
        counter +=1
    #4
    Q_midbins4 = np.empty_like(Q_pos_hist4, dtype='float64')
    size = len(binsQ4)-1
    counter = 0
    while counter < size:
        Q_midbins4[counter] = (binsQ4[counter] + binsQ4[counter+1])/2
        counter +=1
    #5
    Q_midbins5 = np.empty_like(Q_pos_hist5, dtype='float64')
    size = len(binsQ5)-1
    counter = 0
    while counter < size:
        Q_midbins5[counter] = (binsQ5[counter] + binsQ5[counter+1])/2
        counter +=1
    #6
    Q_midbins6 = np.empty_like(Q_pos_hist6, dtype='float64')
    size = len(binsQ6)-1
    counter = 0
    while counter < size:
        Q_midbins6[counter] = (binsQ6[counter] + binsQ6[counter+1])/2
        counter +=1
    #
    #summary tables for each
    #SF1 = Table([SF_midbins1,SF_pos_hist1,SF_neg_hist1],names=('SF Bins','SF pos','SF neg'))
    #SF2 = Table([SF_midbins2,SF_pos_hist2,SF_neg_hist2],names=('SF Bins','SF pos','SF neg'))
    #SF3 = Table([SF_midbins3,SF_pos_hist3,SF_neg_hist3],names=('SF Bins','SF pos','SF neg'))
    #SF4 = Table([SF_midbins4,SF_pos_hist4,SF_neg_hist4],names=('SF Bins','SF pos','SF neg'))
    #SF5 = Table([SF_midbins5,SF_pos_hist5,SF_neg_hist5],names=('SF Bins','SF pos','SF neg'))
    #SF6 = Table([SF_midbins6,SF_pos_hist6,SF_neg_hist6],names=('SF Bins','SF pos','SF neg'))
    #Q1 = Table([Q_midbins1,Q_pos_hist1,Q_neg_hist1],names=('Q Bins','Q pos','Q neg'))
    #Q2 = Table([Q_midbins2,Q_pos_hist2,Q_neg_hist2],names=('Q Bins','Q pos','Q neg'))
    #Q3 = Table([Q_midbins3,Q_pos_hist3,Q_neg_hist3],names=('Q Bins','Q pos','Q neg'))
    #Q4 = Table([Q_midbins4,Q_pos_hist4,Q_neg_hist4],names=('Q Bins','Q pos','Q neg'))
    #Q5 = Table([Q_midbins5,Q_pos_hist5,Q_neg_hist5],names=('Q Bins','Q pos','Q neg'))
    #Q6 = Table([Q_midbins6,Q_pos_hist6,Q_neg_hist6],names=('Q Bins','Q pos','Q neg'))
    #group summary tables
    SFg1 = Table([SF_midbins1,SF_pos_group1,SF_neg_group1],names=('SF Bins','SF pos','SF neg'))
    SFg2 = Table([SF_midbins2,SF_pos_group2,SF_neg_group2],names=('SF Bins','SF pos','SF neg'))
    Qg1 = Table([Q_midbins1,Q_pos_group1,Q_neg_group1],names=('Q Bins','Q pos','Q neg'))
    Qg2 = Table([Q_midbins2,Q_pos_group2,Q_neg_group2],names=('Q Bins','Q pos','Q neg'))
    

##false pos/neg fractions
SF_frac1 = np.empty_like(SF_pos_hist1, dtype='float64')
SF_frac2 = np.empty_like(SF_pos_hist2, dtype='float64')
Q_frac1 = np.empty_like(Q_pos_hist1, dtype='float64')
Q_frac2 = np.empty_like(Q_pos_hist2, dtype='float64')

zero_falseSF = 0
counter = 0         #SF pop
size = len(SF_pos_group1)
while counter < size:
    if SF_pos_group1[counter] ==0 and SF_neg_group1[counter] ==0:
        SF_frac1[counter] = 1
        zero_falseSF +=1
    else: SF_frac1[counter] = SF_pos_group1[counter]/SF_neg_group1[counter]
    counter+=1
zero_falseQ = 0
counter = 0         #Q pop
size = len(Q_pos_group1)
while counter < size:
    if Q_pos_group1[counter] ==0 and Q_neg_group1[counter] ==0:
        Q_frac1[counter] = 1
        zero_falseQ +=1
    else: Q_frac1[counter] = Q_pos_group1[counter]/Q_neg_group1[counter]
    counter+=1
counter = 0         #SF pop
size = len(SF_pos_group2)
while counter < size:
    if SF_pos_group2[counter] ==0 and SF_neg_group2[counter] ==0:
        SF_frac2[counter] = 1
        zero_falseSF +=1
    else: SF_frac2[counter] = SF_pos_group2[counter]/SF_neg_group2[counter]
    counter+=1
zero_falseQ = 0
counter = 0         #Q pop
size = len(Q_pos_group2)
while counter < size:
    if Q_pos_group2[counter] ==0 and Q_neg_group2[counter] ==0:
        Q_frac2[counter] = 1
        zero_falseQ +=1
    else: Q_frac2[counter] = Q_pos_group2[counter]/Q_neg_group2[counter]
    counter+=1
#errors
SF_perr1 = (np.sqrt(SF_pos_group1))/SF_pos_group1    #compute rel. error for false pos/neg
SF_nerr1 = (np.sqrt(SF_neg_group1))/SF_neg_group1
Q_perr1 = (np.sqrt(Q_pos_group1))/Q_pos_group1
Q_nerr1 = (np.sqrt(Q_neg_group1))/Q_neg_group1
SF_perr2 = (np.sqrt(SF_pos_group2))/SF_pos_group2    #compute rel. error for false pos/neg
SF_nerr2 = (np.sqrt(SF_neg_group2))/SF_neg_group2
Q_perr2 = (np.sqrt(Q_pos_group2))/Q_pos_group2
Q_nerr2 = (np.sqrt(Q_neg_group2))/Q_neg_group2
#
SF_frac_error1 = np.sqrt((SF_perr1**2)+(SF_nerr1**2))      #computer error of correction factor (sum in quadrature of rel. errors in false pos/neg)
Q_frac_error1 = np.sqrt((Q_perr1**2)+(Q_nerr1**2))
SF_frac_error2 = np.sqrt((SF_perr2**2)+(SF_nerr2**2))      #computer error of correction factor (sum in quadrature of rel. errors in false pos/neg)
Q_frac_error2 = np.sqrt((Q_perr2**2)+(Q_nerr2**2))
#
# sift out 'NAN' from fractions for 0/0
counter=0                   #SF
size=len(SF_frac_error1)
while counter<size:
    if np.isnan(SF_frac_error1[counter]):
        SF_frac_error1[counter] = 0
    counter+=1
counter=0                   #Q
size=len(Q_frac_error1)
while counter<size:
    if np.isnan(Q_frac_error1[counter]):
        Q_frac_error1[counter] = 0
    counter+=1
counter=0                   #SF
size=len(SF_frac_error2)
while counter<size:
    if np.isnan(SF_frac_error2[counter]):
        SF_frac_error2[counter] = 0
    counter+=1
counter=0                   #Q
size=len(Q_frac_error2)
while counter<size:
    if np.isnan(Q_frac_error2[counter]):
        Q_frac_error2[counter] = 0
    counter+=1
#plot
plt.close()
MC = plt.figure(num=2)
MC.suptitle('Mass Completeness Correction Factors')
#plt.plot(SF_midbins,)
plt.errorbar(SF_midbins1,SF_frac1,yerr=SF_frac_error1, fmt='ob',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.8, mfc='none')
plt.errorbar(Q_midbins1,Q_frac1,yerr=Q_frac_error1, fmt='or',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.8, mfc='none')
plt.errorbar(SF_midbins2,SF_frac2,yerr=SF_frac_error2, fmt='db',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.8, mfc='none')
plt.errorbar(Q_midbins2,Q_frac2,yerr=Q_frac_error2, fmt='dr',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.8, mfc='none')
plt.plot(SF_midbins1,SF_frac1,'--b', linewidth=0.5, label='Star-forming')
plt.plot(Q_midbins1,Q_frac1,'--r', linewidth=0.5, label='Quiescent')
plt.plot(SF_midbins2,SF_frac2,'-.b', linewidth=0.5, label='Star-forming')
plt.plot(Q_midbins2,Q_frac2,'-.r', linewidth=0.5, label='Quiescent')
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
m=(SF_frac1[1] - SF_frac1[0])/(SF_midbins1[1] - SF_midbins1[0])
b = SF_frac1[1] - (SF_midbins1[1]*m)
X = np.linspace(7.4,8.1,50)
y=[]
for x in X:
    y.append(m*x + b)
plt.plot(X,y)
## extrapolating Q at lo-mass
m=(Q_frac1[1] - Q_frac1[0])/(Q_midbins1[1] - Q_midbins1[0])
b = Q_frac1[1] - (Q_midbins1[1]*m)
X = np.linspace(7.4,8.1,50)
y=[]
for x in X:
    y.append(m*x + b)
plt.plot(X,y)
## extrapolating SF at hi-mass
m=(SF_frac1[4] - SF_frac1[3])/(SF_midbins1[4] - SF_midbins1[3])
b = SF_frac1[4] - (SF_midbins1[4]*m)
X = np.linspace(11.6,12.1,50)
y=[]
for x in X:
    y.append(m*x + b)
plt.plot(X,y)
## extrapolating Q at hi-mass
m=(Q_frac1[2] - Q_frac1[1])/(Q_midbins1[2] - Q_midbins1[1])
b = Q_frac1[2] - (Q_midbins1[2]*m)
X = np.linspace(11.6,12.1,50)
y=[]
for x in X:
    y.append(m*x + b)
plt.plot(X,y)
##group 2
m=(SF_frac2[1] - SF_frac2[0])/(SF_midbins2[1] - SF_midbins2[0])
b = SF_frac2[1] - (SF_midbins2[1]*m)
X = np.linspace(7.4,8.8,50)
y=[]
for x in X:
    y.append(m*x + b)
plt.plot(X,y)
## extrapolating Q at lo-mass
m=(Q_frac2[1] - Q_frac2[0])/(Q_midbins2[1] - Q_midbins2[0])
b = Q_frac2[1] - (Q_midbins2[1]*m)
X = np.linspace(7.4,8.1,50)
y=[]
for x in X:
    y.append(m*x + b)
plt.plot(X,y)
## extrapolating SF at hi-mass
m=(SF_frac2[1] - SF_frac2[0])/(SF_midbins2[1] - SF_midbins2[0])
b = SF_frac2[1] - (SF_midbins2[1]*m)
X = np.linspace(10.9,12.1,50)
y=[]
for x in X:
    y.append(m*x + b)
plt.plot(X,y)
## extrapolating Q at hi-mass
m=(Q_frac2[3] - Q_frac2[2])/(Q_midbins2[3] - Q_midbins2[2])
b = Q_frac2[3] - (Q_midbins2[3]*m)
X = np.linspace(11.6,12.1,50)
y=[]
for x in X:
    y.append(m*x + b)
plt.plot(X,y)
#
# list of correction factors, determined by inspection interpolating SMF mass 
# bins between the mass correction bins
#SF_correction = [4.19,3.73,3.27,2.81,2.355,1.89,1.58,1.42,1.25,1.08,1,1,1,1,1,1]
#Q_correction = [1,.02,.13,.24,.43,.695,.96,1.226,1.4965,1.6771,1.40625,1.136,1,1,1,1]
#SF_correction = [2.5859,2.51,2.4295,2.353,2.2734,1.9766,1.585937,1.1953125,1.0625,1.1875,1.3125,1.3625,1.2375,1.1125,0.985,.8625]      for cutoff .03 & .07 (SF/Q)
#Q_correction = [1.20723684,.885,.5625,.2865,.6096,.9309,1.25776,1.5859375,1.625,1.375,1.16183375,1.083706,1.00558,.796875,.578125,.359375]
#
#SF_correction = [2.75189,2.69367,2.6152,2.548,2.47931,2.34368,2.00757,1.67038,1.33367,1.19593,1.12456,1.05547,.981921,.910591]      #S5 Q7
#Q_correction = [.54758,.476075,.404639,.330849,.487969,.643279,.802955,.9656,1.02748,.889034,.750098,.714258,.677714,.643279]
SF1_correction = [1.33951,1.51897,1.67088,1.87285,1.89296,1.53608,1.17704,1.26723,1.80445,2.34026,2.12618,1.58668,1.05274,0.51642]
Q1_correction = [0.19195,0.241308,0.287054,.334004,.383361,.431515,.477261,.478465,.434551,.392783,.348816,.307048,.263082,.221314]
SF2_correction = [5.19604,4.79841,4.40077,3.99942,3.59411,3.19647,2.79512,2.39749,1.99756,1.59903,1.20051,.800854,.397619,0.00568955]
Q2_correction = [.459455,.894147,1.32308,1.75116,2.17428,2.50885,2.55917,2.60627,2.65469,2.3642,1.95127,1.55232,1.14951,.742818]
#
#weighting: take relative fraction of members in each group, SF & Q seperately (i.e. multiply correction factors by fraction, e.g. for SF, take # in group 1 / # in group 2, mult. corrections factors by that fraction)
# note: a better way to do this would be by individual mass bin, but no time
SFG1_count = smem[0][0]+smem[0][3]+smem[0][4]+smem[0][5]
QG1_count = qmem[0][0]+qmem[0][3]+qmem[0][4]+qmem[0][5] 
SFG2_count = smem[0][1]+smem[0][2]
QG2_count = qmem[0][1]+qmem[0][2]
SF1_scale = SFG1_count / np.sum(smem)
Q1_scale = QG1_count / np.sum(qmem)
SF2_scale = SFG2_count / np.sum(smem)
Q2_scale = QG2_count / np.sum(qmem)
#apply scaling factors
size = len(SF1_correction)
counter = 0
while counter < size:
    SF1_correction[counter] = SF1_correction[counter]*SF1_scale
    counter +=1
size = len(Q1_correction)
counter = 0
while counter < size:
    Q1_correction[counter] = Q1_correction[counter]*Q1_scale
    counter +=1
size = len(SF2_correction)
counter = 0
while counter < size:
    SF2_correction[counter] = SF2_correction[counter]*SF2_scale
    counter +=1
size = len(Q2_correction)
counter = 0
while counter < size:
    Q2_correction[counter] = Q2_correction[counter]*Q2_scale
    counter +=1
###     NOTE: relative fraction calculated above at end of Seciont (ii), and are
###     stored in arrays "frac_smf" & "frac_field_smf"
#
total_smf = np.empty_like(SF_smf, dtype='float64')
##  apply correction factors to SMFs:
SF_smfg1 = SF_smf1+SF_smf4+SF_smf5+SF_smf6
SF_smfg1 = SF_smfg1/SF1_correction
Q_smfg1 = Q_smf1+Q_smf4+Q_smf5+Q_smf6
Q_smfg1 = Q_smfg1/Q1_correction
SF_smfg2 = SF_smf2+SF_smf3
SF_smfg2 = SF_smfg2/SF2_correction
Q_smfg2 = Q_smf2+Q_smf3
Q_smfg2 = Q_smfg2/Q2_correction
total_smf = SF_smfg1 + Q_smfg1 + SF_smfg2 + Q_smfg2
SF_plot = SF_smfg1 + SF_smfg2
Q_plot = Q_smfg1 + Q_smfg2
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
        frac_smf[0][counter] =  (SF_smfg1[counter] + SF_smfg2[counter]) / total_smf[counter]
        frac_smf[1][counter] =  (Q_smfg1[counter]+Q_smfg2[counter]) / total_smf[counter]
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
SF_error = np.sqrt(SF_plot)        #assign errors by population as sqrt of count in each mass bin
Q_error = np.sqrt(Q_plot)
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
xa = plt.errorbar(SF_midbins,SF_plot,yerr=SF_error, fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)#yerr=SF_error,
xb = plt.errorbar(Q_midbins,Q_plot,yerr=Q_error,fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
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
