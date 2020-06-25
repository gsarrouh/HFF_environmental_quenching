#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 14:28:35 2020

@author: gsarrouh
"""
################## master_smfz_8 ##################
## This program will plot the Stellar Mass Function (SMF) for the master_dadta 
## file of all six clusters (most current version: 'master_data_6_final.py); 
## two plots (SF & Q) segregated between spectroscopic and photometric subsamples
#
## v2 includes the parallel fields data; 
## v3 commented out to produce conference 
##    plot. to re-insert, remove all #s and change plot panels to include || 
##    field in cenre column (positions 1 & 4)
## v4 creates SMFs by cluster, and calculates completeness correction by cluster for false pos/neg
##    and attempts compute the completeness correction on a cluster-by-cluster basis.
##    this approach was abandoned in favour of a single set of correciton factors 
##    for the entire sample, justified by the fact that mass completeness limits are nearly the 
##    same for all clusters. 
## v5 removes individual cluster code for mass correction in section (v), plots
##    correction factors for entire sample treated as single population
## v6 moves all totals and SF/Q fractions until AFTER the correction factors 
##    have been made, as this changes the counts in each hist & associated poissonian error
## 
## v8 COMPLETE OVERHAUL: lists are first sorted by cluster, completeness corrections computed 
##    and applied, then each cluster population (both SF & Q) are normalized by the total 
##    number of galaxies in the cluster (i.e. sum(SF)+sum(Q)). 
##
##
##
#
### Section summary:
#
### PROGRAM START
#
### PLOTS THE STELLAR MASS FUNCTION:
### (1)    collect masses into SORTED arrays for SF & Q;
### (1.1)   add DIAG_FLAG_1: summarize sorted arrays;
### (1.2)   import FIELD data;
### (2)    bin into HISTOGRAMS;
### (2.1)    compute bin mid-points; add DIAG_FLAG_2: sub-samples check;
### (3)    CORRECTIONS to raw counts;
### (3.1)   calculate limiting MASS COMPLETENESS correction for low-mass bins; add DIAG_FLAG_3: 
###         display mass correction factors;
### (3.2)   calculate false pos/neg (i.e. spectroscopic) completeness correction for all bins;
###         includes FIGURE for SPECTROSCOPIC COMPLETENESS; add DIAG_FLAG_4: display 
###         corrections for different bin #s; add DIAG_FLAG_5: # of false pos/neg check 
###         & display spec correction factors;
### (3.3)  NORMALIZE; add DIAG_FLAG_6: check normalization result
### (4)    add ERROR BARS to scatter plot;
### (5)    EMCEE simulation; see emcee_chi2_final.py;
### (6)    build SCHECTER best-fit models;
### (7)    PLOT that shit;
#
### PROGRAM END
#
## NOTE: there are flags for diagnostics and plotting throughout the script. search "MAY NEED TO EDIT" to identify where these flags are
#
#
#
###################     PROGRAM START
#
## TIME_FLAG: START
## superior time_flag which supercedes all others and times the entire program
time_flag = 0     # track & print time to execute current section
#
if time_flag == 1:
    start_time = time.time()
#  
# Import modules
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy
from astropy.table import Table
from astropy.table import Column
from scipy.optimize import curve_fit
#
# this line is specific to jupyter notebook, and allows for Figure editting in a GUI, instead of inline
%matplotlib qt  
#
#
## MASTER DIAGNOSTIC FLAG: allows you to "turn off" all diagnostics at once (if equals 0), otherwise diagnostic flags are set one at a time (search "MAY NEED TO EDIT")
diag_flag_master = 1
#
## SECTION (1): collect objects above limiting mass by cluster into a single array in order to plot 
## SF_*/Q_*, and track sublists of objects which have spec vs those which only have phot, separately for SF/Q; creates list of samples to be binned & plotted as histogram/scatterplot
#
## Cluster sample
#
SF_list = [[],[],[],[],[],[]]    # create an empty list filled with 6 empty lists - 1 for each cluster
Q_list = [[],[],[],[],[],[]]
SF_phot_list = [[],[],[],[],[],[]]
SF_spec_list = [[],[],[],[],[],[]]
Q_phot_list = [[],[],[],[],[],[]]
Q_spec_list = [[],[],[],[],[],[]]
#
limiting_mass = [7.5,7.8,8.0,7.5,7.4,7.3] # clusters 1,2,3,4,5,6, see IDs below
#
SF_pos_lost = np.array([0]*6)        # to track SF/Q false pos/neg objects lost due to their being below the mass limit, by cluster
SF_neg_lost = np.array([0]*6)
Q_pos_lost = np.array([0]*6)
Q_neg_lost = np.array([0]*6)
other_lost = np.array([0]*6)    #objects below limiting mass other than false pos/neg
#
# The following loop searches the master catalogue 'master_cat' and separates all objects by 
# cluster. Then, it looks for all objects above the limiting mass for that cluster. It then 
# creates two lists: one for SF and one for Q (e.g. SF_*/Q_*). It further splits these lists into those objects with spectrscopy, and those without (e.g. SF*_spec/phot)
#
for cluster in range(len(limiting_mass)):          # loop through clusters one at a time; "cluster" takes on values [0,1,2,3,4,5]
    for counter in range(len(master_cat)):
        if master_cat['cluster'][counter] == (cluster+1):    # cluster #
            if master_cat['lmass'][counter] > limiting_mass[cluster]:    # limiting mass of cluster: 7.5
                if master_cat['member'][counter] == 0:    # cluster member = 0
                    if master_cat['type'][counter] == 1:   # SF type = 1
                        SF_list[cluster].append(master_cat['lmass'][counter])
                        if master_cat['sub'][counter]==2:     # sub=2 is objects w/ photometry only
                            SF_phot_list[cluster].append(master_cat['lmass'][counter])
                        elif master_cat['sub'][counter]==1 or master_cat['sub'][counter]==3:  # sub=1 is objects w/ both spec & phot; sub=3 is for spec only
                            SF_spec_list[cluster].append(master_cat['lmass'][counter])
                    elif master_cat['type'][counter] == 2: # Q type = 2
                        Q_list[cluster].append(master_cat['lmass'][counter])
                        if master_cat['sub'][counter]==2:     # sub=2 is objects w/ photometry only
                            Q_phot_list[cluster].append(master_cat['lmass'][counter])
                        elif master_cat['sub'][counter]==1 or master_cat['sub'][counter]==3:  # sub=1 is objects w/ both spec & phot; sub=3 is for spec only
                            Q_spec_list[cluster].append(master_cat[counter]['lmass'])
            elif master_cat['lmass'][counter] < limiting_mass[cluster]:
                if master_cat['member'][counter] == 2:    # member false pos = 2
                    if master_cat['type'][counter] == 1:   # SF type = 1
                        SF_pos_lost[cluster]+=1
                    elif master_cat['type'][counter] == 2: # Q type = 2
                        Q_pos_lost[cluster]+=1
                elif master_cat['member'][counter] ==3:   # member false neg = 3
                    if master_cat['type'][counter] == 1:   # SF type = 1
                        SF_neg_lost[cluster]+=1
                    elif master_cat['type'][counter] == 2: # Q type = 2
                        Q_neg_lost[cluster]+=1
                else: other_lost[cluster]+=1
#
#
## SECTION (1.1): display summary
#
### MAY NEED TO EDIT ### diag_flag_1
##  displays summary of sorted lists above
diag_flag_1 = 1             # 0=off (don't display diagnostic); 1=on (display diagnostic table)
#
if diag_flag_1 == 1 and diag_flag_master == 1:
    ## Summarize initial data stats in table
    #
    # calculate list lengths
    SF_len = np.array([0]*6)
    Q_len = np.array([0]*6)
    SF_phot_len = np.array([0]*6)
    SF_spec_len = np.array([0]*6)
    Q_phot_len = np.array([0]*6)
    Q_spec_len = np.array([0]*6)
    #
    for ii in range(len(SF_list)):
        SF_len[ii] = len(SF_list[ii])
        Q_len[ii] = len(Q_list[ii])
        SF_phot_len[ii] = len(SF_phot_list[ii])
        SF_spec_len[ii] = len(SF_spec_list[ii])
        Q_phot_len[ii] = len(Q_phot_list[ii])
        Q_spec_len[ii] = len(Q_spec_list[ii])
    #
    data_names = Column(['SF','Q','SF_phot','SF_spec','Q_phot','Q_spec','Total'],name='Property')
    # 
    # setup arrays for displaying in table
    col_names=['Property','Total','macs0416','macs1149','macs0717','abell370','abell1063','abell2744']
    SF_tabular = ['SF',np.sum(SF_len),SF_len[0],SF_len[1],SF_len[2],SF_len[3],SF_len[4],SF_len[5]]
    Q_tabular = ['Q',np.sum(Q_len),Q_len[0],Q_len[1],Q_len[2],Q_len[3],Q_len[4],Q_len[5]]
    SF_phot_tabular = ['SF_phot',np.sum(SF_phot_len),SF_phot_len[0],SF_phot_len[1],SF_phot_len[2],SF_phot_len[3],SF_phot_len[4],SF_phot_len[5]]
    SF_spec_tabular = ['SF_spec',np.sum(SF_spec_len),SF_spec_len[0],SF_spec_len[1],SF_spec_len[2],SF_spec_len[3],SF_spec_len[4],SF_spec_len[5]]
    Q_phot_tabular = ['Q_phot',np.sum(Q_phot_len),Q_phot_len[0],Q_phot_len[1],Q_phot_len[2],Q_phot_len[3],Q_phot_len[4],Q_phot_len[5]]
    Q_spec_tabular = ['Q_spec',np.sum(Q_spec_len),Q_spec_len[0],Q_spec_len[1],Q_spec_len[2],Q_spec_len[3],Q_spec_len[4],Q_spec_len[5]]
    # display table
    print('\nSection 1: Cluster MEMBERSHIP stats: ')
    from tabulate import tabulate
    print(tabulate([SF_tabular,Q_tabular,SF_phot_tabular,SF_spec_tabular,Q_phot_tabular,Q_spec_tabular],headers=col_names))
    #
    print('\nObjects lost below limiting mass: ')
    print('SF false pos.: ',str(np.sum(SF_pos_lost)))
    print('SF false neg.: ',str(np.sum(SF_neg_lost)))
    print('Q false pos.: ',str(np.sum(Q_pos_lost)))
    print('Q false neg.: ',str(np.sum(Q_neg_lost)))
    print('Other (objects below limiting mass other than false pos/neg): ',str(np.sum(other_lost)))
    print(other_lost)
#
#
#
#
## SECTION (1.2): Field sample
## 
#
#
##### IMPORT DATA FROM AM's EMAIL 06/19/20
#
#
#
#
#
## SECTION (2): sort objects into HISTOGRAMS bins for both SF & Q populations, then sum 
## for 'total' population. then normalize each cluster SMF by the total cluster mass. compute midbins. use for total pop plot, and for relative fractions
#
## cluster populations arrays: (SF/Q/total)_smf & (SF/Q/total)_field_smf
#
range2 = [7.3,12.3]     #sets range for all histrograms to be computer: cluster,field,false pos/neg
bin_width = 0.2  # in dex
num_points = (round((range2[1]-range2[0])/bin_width))+1       # compute # of data points
num_bins = np.linspace(range2[0],range2[1],num_points)
#
## smf histograms for individual clusters
#
SF_raw_smf = [[],[],[],[],[],[]]       # initialize list of lists to store histograms of SMFs
Q_raw_smf = [[],[],[],[],[],[]]
#
for ii in range(len(SF_list)):
    SF_raw, mass_bins = np.histogram(SF_list[ii], bins=num_bins,range=range2)
    Q_raw, mass_bins = np.histogram(Q_list[ii], bins=num_bins,range=range2)
    SF_raw_smf[ii].append(SF_raw)
    Q_raw_smf[ii].append(Q_raw)
#
total_raw_smf = SF_raw_smf + Q_raw_smf
#
# Display some data for total, SF, Q: 
print('\nSection 2: RAW totals')
print('SF: ',str(np.sum(SF_raw_smf)))
print('Q: ',str(np.sum(Q_raw_smf)))
print('Total: ',str(np.sum(total_raw_smf)),'\n')
#
#
## section (2.1): compute MIDBINS
#
## find midpoint of hist. bins. all populations have been binned identically, so the one 'midbin' will serve for all data arrays to be plotted. for visual clarity when plotting, offset the Q_midpoints by delta_x = 0.05
#
#
## define a function to compute the mid-points of the bins from a histogram
def midbins(bins):
    size = len(bins)-1
    x_midbins = np.empty([size,1],dtype='float64')
    for x in range(size):
        x_midbins[x] = (bins[x] + bins[(x+1)])/2
    return x_midbins
#
# compute midbins        
SF_midbins = midbins(mass_bins)
# offset Q midbins for plotting clarity
Q_midbins = SF_midbins + 0.05
#
#
## SORT the spec/phot subsamples into histograms as well, and confirm that spec + phot = total in each mass bin for each type of galaxy
#
# sort spec/phot subsamples into histograms for each cluster
#
SF_phot_smf = [[],[],[],[],[],[]]       # initialize list of lists to store histograms of SMFs
SF_spec_smf = [[],[],[],[],[],[]]
Q_phot_smf = [[],[],[],[],[],[]]       
Q_spec_smf = [[],[],[],[],[],[]]
#
for ii in range(len(SF_spec_list)):
    SF_spec, mass_bins = np.histogram(SF_spec_list[ii], bins=num_bins,range=range2)
    SF_phot, mass_bins = np.histogram(SF_phot_list[ii], bins=num_bins,range=range2)
    Q_spec, mass_bins = np.histogram(Q_spec_list[ii], bins=num_bins,range=range2)
    Q_phot, mass_bins = np.histogram(Q_phot_list[ii], bins=num_bins,range=range2)
    SF_spec_smf[ii].append(SF_spec)
    SF_phot_smf[ii].append(SF_phot)
    Q_spec_smf[ii].append(Q_spec)
    Q_phot_smf[ii].append(Q_phot)
#
# convert SMF lists to arrays
SF_raw_smf = np.array(SF_raw_smf)
Q_raw_smf = np.array(Q_raw_smf)
total_raw_smf = np.array(total_raw_smf)
SF_phot_smf = np.array(SF_phot_smf)
SF_spec_smf = np.array(SF_spec_smf)
Q_phot_smf = np.array(Q_phot_smf)
Q_spec_smf = np.array(Q_spec_smf)
#
### MAY NEED TO EDIT: diag_flag_2
## DIAGNOSTIC: add spec & phot subsampes together for each cluster, and ensure they equal the total raw count in each mass bin
diag_flag_2 = 0           # 0=off, i.e. don't do diagnostic; 1=on, i.e. perform diagnostic
#
if diag_flag_2 == 1 and diag_flag_master == 1:
    # compute differences, e.g.: SF_smf1 = SF1_spec_smf + SF1_phot_smf for each mass bin. they should be the same
    SF_diff = np.array([[0]*len(SF_midbins)]*6)     # initialize array to store difference between sample & sub-samples, by cluster
    Q_diff = np.array([[0]*len(Q_midbins)]*6)
    #
    print('Section 2.1: Differences between raw cluster count and (spec + phot) subsamples, by cluster')
    for ii in range(len(SF_raw_smf)):
        SF_diff[ii] = SF_raw_smf[ii] - (SF_phot_smf[ii] + SF_spec_smf[ii])
        Q_diff[ii] = Q_raw_smf[ii] - (Q_phot_smf[ii] + Q_spec_smf[ii])
        print('SF',str(ii+1),' difference: ',str(np.sum(SF_diff[ii])))
        print('Q',str(ii+1),' difference: ',str(np.sum(Q_diff[ii])))
        print('Total difference: %s'%np.sum([SF_diff,Q_diff]),'\n')
#    
#
#
#
#
## SECTION (3): calculate corrections to raw counts. There are two separate corrections - one for limiting mass completeness (i.e. correct for the fact that not all clusters are complete down to our lowest mass bin), and one for spectrscopic completeness (i.e. correct for the fact that there are false pos/neg objects in the sample of galaxies which have both spec & phot, and make correction to photometric sample to account for the ratio of false pos/neg)
#
#
## SECTION (3.1): calculate MASS COMPLETENESS corrections. compute correction to low-mass bin points due to varying mass completenesses of each cluster. The correction factor is: (total # of clusters) / (# of clusters complete at that mass bin), and will be multiplied by the raw number count of galaxies in each mass bin. 
#
## an examination of the limiting mass for each cluster (see list "limiting_mass", above) shows that the following bin midpoints have the following corresponding number of clusters complete at that mass: [7.3,7.5,7.7,7.9,8.1] ---> [1,4,4,5,6]. So all 6 clusters are complete at a mass  of 8.1, but only 1 cluster is complete down to 7.3. the corresponding corrections are as follows:
#
#mass_completeness_correction = np.array([6,1.5,1.5,1.2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
#
## if we assume all clusters have the same general composition of SF to Q galaxies, i.e. the same relative fraction in each cluster, then if only one cluster is complete, we "scale it up" by a factor of 6 (# of clusters total / # of clusters complete at that mass). if 4 clusters are complete, we scale it up by (6/4) = 1.5. In this way, it's as if each cluster were complete down to the limiting mass of 7.3
#
## the following loop automates the above explanation, making the code adaptable should i decide to change the number of bins in future
#
mass_completeness_correction = np.zeros_like(SF_midbins)
for ii in range(len(mass_completeness_correction)):
    for jj in range(len(limiting_mass)):
        if limiting_mass[jj] <= SF_midbins[ii]:    # count # of clusters complete at each mass bin
            mass_completeness_correction[ii]+=1
mass_completeness_correction = np.transpose(6/mass_completeness_correction)  # the correction factor is: (total # of clusters) / (# of clusters complete at that mass bin); return as a row vector
#
#
# define a function to compute the incremental difference of applying a correction to the raw count lists
def correction_difference(raw_smf,completeness_correction):
    corrected_smf = raw_smf*completeness_correction
    diff = corrected_smf- raw_smf
    return diff
#
# compute how many objects are added to each mass bin as a result of applying the mass_completeness_correction to the *_raw_smf lists. confirm that (# added to SF) + (# added to Q) = (# added to total)
SF_mass_completeness_diff = correction_difference(SF_raw_smf,mass_completeness_correction)
Q_mass_completeness_diff = correction_difference(Q_raw_smf,mass_completeness_correction)
total_mass_completeness_diff = correction_difference(total_raw_smf,mass_completeness_correction)
#
### MAY NEED TO EDIT: diag_flag_3
diag_flag_3 = 0
#
if diag_flag_3 == 1 and diag_flag_master == 1:
# Display correction factors
    print('/nSection 3.1: Mass completeness correction factors by bin: ',str(mass_completeness_correction),'\n')
    #
    # Display some data for total, SF, Q: 
    print('Section 3.1: Galaxies added due to MASS COMPLETENESS correction')
    print('SF: ',str(SF_mass_completeness_diff),'\nor ',str((np.sum(SF_mass_completeness_diff)/np.sum(SF_raw_smf))*100),'%\n')
    print('Q: ',str(Q_mass_completeness_diff),'\nor ',str((np.sum(Q_mass_completeness_diff)/np.sum(Q_raw_smf))*100),'%\n')
    print('Total: ',str(total_mass_completeness_diff),'\nor ',str((np.sum(total_mass_completeness_diff)/np.sum(total_raw_smf))*100),'%\n')
else:
    print('Section 3.1: Galaxies added due to MASS COMPLETENESS correction\nSF: %s'%np.sum(SF_mass_completeness_diff),'\nQ: %s'%np.sum(SF_mass_completeness_diff))
#
#
#
#
## SECTION (3.2): calculate SPECTROSCOPIC COMPLETENESS correction. basically, look at all the false positives/false negatives, and sort them by type (i.e. SF/Q). then bin them (i.e. make histograms of false pos/neg for each of SF/Q). take their ratio of false pos to false neg, and plot that ratio. it is the correction factor to be applied to the photometric subsample
#
#
SF_pos = []
SF_neg = []
Q_pos = []
Q_neg = []
pos_by_cluster = np.array([[0]*6]*2)    #for tracking false pos/neg by cluster; row_1=SF, row_2=Q
neg_by_cluster = np.array([[0]*6]*2)
objects_below_lim_mass = np.array([0]*6)    # for tracking objects below the limiting mass of each cluster
#
for counter in range(len(master_cat)):
    for ii in range(len(limiting_mass)):
        if master_cat['cluster'][counter] == (ii+1):           # only look at objects above the limiting mass for each cluster
            if master_cat['lmass'][counter] > limiting_mass[ii]:      
                if master_cat['type'][counter] == 1:      # type=1 for SF
                    if master_cat['member'][counter] == 2:     # member=2 for false pos
                        SF_pos.append(master_cat['lmass'][counter])
                        for ii in range(len(pos_by_cluster[0])):
                            if master_cat['cluster'][ii] == (ii+1):
                                pos_by_cluster[0]+=1           # track false pos for SF
                    elif master_cat['member'][counter] == 3:   # member=3 for false neg
                        SF_neg.append(master_cat['lmass'][counter])
                        for ii in range(len(pos_by_cluster[0])):
                            if master_cat['cluster'][ii] == (ii+1):
                                neg_by_cluster[0]+=1           # track false neg for SF
                elif master_cat['type'][counter] == 2:     # type=2 for Q
                    if master_cat['member'][counter] == 2:     # member=2 for false pos
                        Q_pos.append(master_cat['lmass'][counter])
                        for ii in range(len(pos_by_cluster[1])):
                            if master_cat['cluster'][ii] == (ii+1):
                                pos_by_cluster[1]+=1           # track false pos for Q
                    elif master_cat['member'][counter] == 3:   # member=3 for false neg
                        Q_neg.append(master_cat['lmass'][counter])
                        for ii in range(len(neg_by_cluster[1])):
                            if master_cat['cluster'][ii] == (ii+1):
                                neg_by_cluster[1]+=1           # track false neg for Q
            else: 
                objects_below_lim_mass[ii]+=1
#
### bin SF & Q, then compute false pos/neg fractions by mass bin for correction factors. 
### NOTE: that # of bins was determined as the largest number which would ensure that all bins are populated
num_bins_array = [8,7,6,5,4,3,2]   # number of histogram (i.e. SMF) mass bins to try
#
#
### MAY NEED TO EDIT: diag_flag_4
# SPEC. BINNING: iterate through different number of histogram bins to see which yields a set of corrections closest in general to ~1
diag_flag_4 = 0
#
if diag_flag_4 == 1 and diag_flag_master == 1:
    # write a loop that interatively uses a different number of bins in the histrogram, to isolate the largest number for which all bins have at least one entry; NOTE: the lines that stop the loop have been commented out, to investigate the relative fraction of false pos/neg for each different # of bins
    print('Section 3.2: Spec. completeness correction binning\n')
    #SF
    for number in range(len(num_bins_array)):
        print('# of bins to try (SF): %s'%num_bins_array[number])
        # make histograms
        SF_pos_hist, bins_SF = np.histogram(SF_pos, bins=num_bins_array[number], range=range2)
        SF_neg_hist, bins_SF = np.histogram(SF_neg, bins=num_bins_array[number], range=range2)
        print('SF_pos: %s'%SF_pos_hist)
        print('SF_neg: %s'%SF_neg_hist)
        print('SF Bins: %s'%bins_SF)
        for jj in range(len(SF_pos_hist)):
            if SF_pos_hist[jj] == 0 and SF_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
                SF_pos_hist[jj] = 1
                SF_neg_hist[jj] = 1
        SF_fraction = SF_pos_hist / SF_neg_hist
        print('SF_fraction: %s'%SF_fraction)
        total = np.sum(SF_pos_hist==0) + np.sum(SF_neg_hist==0)
        #if total == 0:
        #    num_binsSF = num_bins_array[number]         # set number of SF bins as greatest number for which each bin is populated
        #    print('# of SF bins: %s'%num_bins_array[number])
        #    break
    # Q
    for number in range(len(num_bins_array)):
        print('# of bins to try (Q): %s'%num_bins_array[number])
        # make histograms
        Q_pos_hist, bins_Q = np.histogram(Q_pos, bins=num_bins_array[number], range=range2)
        Q_neg_hist, bins_Q = np.histogram(Q_neg, bins=num_bins_array[number], range=range2)
        print('Q_pos: %s'%Q_pos_hist)
        print('Q_neg: %s'%Q_neg_hist)
        print('Q Bins: %s'%bins_Q)
        for jj in range(len(Q_pos_hist)):
            if Q_pos_hist[jj] == 0 and Q_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
                Q_pos_hist[jj] = 1
                Q_neg_hist[jj] = 1
        Q_fraction = Q_pos_hist / Q_neg_hist
        print('Q_fraction: %s'%Q_fraction,'\n')
        total = np.sum(Q_pos_hist==0) + np.sum(Q_neg_hist==0)
        #if total == 0:
        #    num_binsQ = num_bins_array[number]           # set number of Q bins as greatest number for which each bin is populated
        #    print('# of Q bins: %s'%num_bins_array[number])
    #    break
#
###### The following few lines are for if you want to choose the # of bins independent of the criteria described above (i.e. if you want to try fewer bins than the largest # for which each bin is populated, or if you want to specify non-symmetric bin widths); marked by 5 hashtags #####
###### recall: range2 = [7.3,12.3]
num_binsSF = [7.3,8.55,9.8,11.14,12.3]#3   ### ASYMMETRIC BINNING: even bins would be [ 7.3   8.55  9.8  11.05 12.3 ]; i made the last bin a bit smaller so that it is empty for SF false pos, given there are no hi-mass false neg. (max SF_neg = 10.63, max SF_pos = 11.13)
#####num_binsQ = [7.3,8.85,10.1,11.05,12.3]#4
#
###num_binsSF = 5
num_binsQ = 5
#
SF_pos_hist, bins_SF = np.histogram(SF_pos, bins=num_binsSF, range=range2)
SF_neg_hist, bins_SF = np.histogram(SF_neg, bins=num_binsSF, range=range2)
Q_pos_hist, bins_Q = np.histogram(Q_pos, bins=num_binsQ, range=range2)
Q_neg_hist, bins_Q = np.histogram(Q_neg, bins=num_binsQ, range=range2)
#
#####print('SF: %s'%bins_SF)
#####print('Q: %s'%bins_Q)
#
### MAY NEED TO EDIT: diag_flag_5
# display diagnostics
diag_flag_5 = 0              # 0=off, 1=on
#
if diag_flag_5 == 1 and diag_flag_master == 1:
    # sum into total list, to compare with totals reported in "spec_stats1" table from "master_data_*.py"
#    total_pos_hist = SF_pos_hist + Q_pos_hist     
#    total_neg_hist = SF_neg_hist + Q_neg_hist
    # print
    print('Section 3.2: spec.completeness correction:\nThe data preparation file reports:')
    print('# of SF false pos: ',str(np.sum(pos[0])))
    print('# of SF false neg: ',str(np.sum(neg[0])))
    print('# of Q false pos: ',str(np.sum(pos[1])))
    print('# of Q false neg: ',str(np.sum(neg[1])))
    print('This SMF preparation file finds:')
    print('# of SF false pos: ',str(np.sum(SF_pos_hist)))
    print('# of SF false neg: ',str(np.sum(SF_neg_hist)))
    print('# of Q false pos: ',str(np.sum(Q_pos_hist)))
    print('# of Q false neg: ',str(np.sum(Q_neg_hist)))
    print('\nFalse pos/neg lost below limiting mass: ')
    print('SF false pos.: ',str(np.sum(SF_pos_lost)))
    print('SF false neg.: ',str(np.sum(SF_neg_lost)))
    print('Q false pos.: ',str(np.sum(Q_pos_lost)))
    print('Q false neg.: ',str(np.sum(Q_neg_lost)),'\n')
    print('# objects lost below limiting mass (by cluster): ',objects_below_lim_mass)
#
## compute false pos/ false neg ratio; there is a diagnostic built in for error handling - since we require that all mass bins be populated by at least one false pos. and one false neg. (so that we may compute their ratio), the program BREAKS when an empty mass bin is encountered, and you are prompted to try a new number of bins 
SF_frac = np.empty_like(SF_pos_hist, dtype='float32')
Q_frac = np.empty_like(Q_pos_hist, dtype='float32')
# compute fractions for SF, exiting loop if a bin value of zero is encountered; somewhat deprecated given the loop added above, which tests different # of bins for SF/Q. 
for ii in range(len(SF_pos_hist)):               
    if SF_pos_hist[ii] == 0 and SF_neg_hist[ii] == 0:
#        print('Zero false pos. AND zero false neg in bin',str(ii+1))
#        print('Adjust the value of "num_binsSF"\n')
#        break
        SF_frac[ii] = 1
        pass
    elif SF_pos_hist[ii] == 0:
        print('Zero false pos. in bin',str(ii+1))
        print('Adjust the value of "num_binsSF"\n')
        break
    elif SF_neg_hist[ii] == 0:
        print('Zero false neg. in bin',str(ii+1))
        print('Adjust the value of "num_binsSF"\n')
        break
    else:
        SF_frac[ii] = SF_pos_hist[ii] / SF_neg_hist[ii]
# compute fractions for Q, exiting loop if a bin value of zero is encountered 
for ii in range(len(Q_pos_hist)):               
    if Q_pos_hist[ii] == 0 and Q_neg_hist[ii] == 0:
#        print('Zero false pos. AND zero false neg in bin',str(ii+1))
#        print('Adjust the value of "num_binsQ"\n')
#        break
        Q_frac[ii] = 1
        pass
    elif Q_pos_hist[ii] == 0:
        print('Zero false pos. in bin',str(ii+1))
        print('Adjust the value of "num_binsQ"\n')
        break
    elif Q_neg_hist[ii] == 0:
        print('Zero false neg. in bin',str(ii+1))
        print('Adjust the value of "num_binsQ"\n')
        break
    else:
        Q_frac[ii] = Q_pos_hist[ii] / Q_neg_hist[ii]        
#
# compute midbins for spec. mass completeness plot (i.e. plot of false pos/false neg ratios)
SF_frac_midbins = midbins(bins_SF)
Q_frac_midbins = midbins(bins_Q)
#
#
## now compute the errors for the spec. completeness plot, which is simply sqrt(N) since the spectroscopic uncertainty is Poissonian in nature. do so by computing the relative error for the false pos & false neg histograms for each of SF/Q, and then sum in quadrature to determine relative error of fractions
#              
SF_relerr_pos = (np.sqrt(SF_pos_hist))/SF_pos_hist
SF_relerr_neg = (np.sqrt(SF_neg_hist))/SF_neg_hist
Q_relerr_pos = (np.sqrt(Q_pos_hist))/Q_pos_hist
Q_relerr_neg = (np.sqrt(Q_neg_hist))/Q_neg_hist
#              
SF_frac_err = np.sqrt((SF_relerr_pos**2) + (SF_relerr_neg**2))*SF_frac              
Q_frac_err = np.sqrt((Q_relerr_pos**2) + (Q_relerr_neg**2))*Q_frac              
#              
#    
#
## FIGURE ##
#
### MAY NEED TO EDIT: plot_flag
plot_flag = 0        # 0=off (i.e. don't make plot), 1=on (i.e. make plot)
if plot_flag == 1:
    # plot Spectroscopic completion correction factors 
    plt.close()
    MC = plt.figure(num=2)
    #MC.suptitle('Spectroscopic Completeness Correction Factors')
    plt.errorbar(SF_frac_midbins,SF_frac,yerr=SF_frac_err, fmt='ob',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.8, mfc='none')
    plt.errorbar(Q_frac_midbins,Q_frac,yerr=Q_frac_err, fmt='or',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.8, mfc='none')
    plt.plot(SF_frac_midbins,SF_frac,'-b', linewidth=1.0, label='Star-forming')
    plt.plot(Q_frac_midbins,Q_frac,'-r', linewidth=1.0, label='Quiescent')
    plt.plot([0,13],[1,1],'--k',linewidth = 0.5)
    plt.legend(loc='upper right', frameon=False)
    plt.xlim=(7.1,12.5)
    plt.xlabel('$log(M/M_{\odot})$')
    plt.ylim=(-0.5,4.1)
    plt.ylabel('Correction factor\n(false pos / false neg)')
    plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on', labelleft='on')
    plt.minorticks_on()
    plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
    MC.xlim=(7.1,12.5)
    #
## Now interpolate/extrapolate between these data points
#
# initialize arrays to store slopes/intercepts for extrapolation/interpolation of spec mass completeness correction factors
m_SF = np.zeros((len(SF_frac_midbins)-1))     
b_SF = np.zeros((len(SF_frac_midbins)-1))
m_Q = np.zeros((len(Q_frac_midbins)-1))     
b_Q = np.zeros((len(Q_frac_midbins)-1))
SF_spec_completeness_correction = np.zeros_like(SF_midbins,dtype='float32')
Q_spec_completeness_correction = np.zeros_like(Q_midbins,dtype='float32')
#
## SF
for ii in range(len(SF_frac_midbins)-1):
    m_SF[ii] = (SF_frac[ii+1] - SF_frac[ii]) / (SF_frac_midbins[ii+1] - SF_frac_midbins[ii]) # calc slope
    b_SF[ii] = SF_frac[ii] - (SF_frac_midbins[ii]*m_SF[ii])   # calc intercept
#
for ii in range(len(SF_midbins)):
    if SF_spec_completeness_correction[ii] == 0:     # don't overwrite cell once correction factor is computed
        if SF_midbins[ii] < SF_frac_midbins[0]:      # extrapolate below lowest mass bin
            SF_spec_completeness_correction[ii] = m_SF[0]*SF_midbins[ii] + b_SF[0]    
        elif SF_midbins[ii] > SF_frac_midbins[-1]:    # extrapolate above highest mass bin
            SF_spec_completeness_correction[ii] = m_SF[-1]*SF_midbins[ii] + b_SF[-1]    
        elif SF_midbins[ii] > SF_frac_midbins[0] and SF_midbins[ii] < SF_frac_midbins[-1]:    # interpolate in between all other points
            for jj in range(len(SF_frac_midbins)-1):
                if SF_midbins[ii] > SF_frac_midbins[jj] and SF_midbins[ii] < SF_frac_midbins[jj+1]:
                    SF_spec_completeness_correction[ii] = m_SF[jj]*SF_midbins[ii] + b_SF[jj]
        else:
            print('Error in SF spec completeness correction computation. ABORT')
            break   
#

## Q
for ii in range(len(Q_frac_midbins)-1):
    m_Q[ii] = (Q_frac[ii+1] - Q_frac[ii]) / (Q_frac_midbins[ii+1] - Q_frac_midbins[ii]) # calc slope
    b_Q[ii] = Q_frac[ii] - (Q_frac_midbins[ii]*m_Q[ii])   # calc intercept
#
for ii in range(len(Q_midbins)):
    if Q_spec_completeness_correction[ii] == 0:     # don't overwrite cell once correction factor is computed
        if Q_midbins[ii] < Q_frac_midbins[0]:      # extrapolate below lowest mass bin
            Q_spec_completeness_correction[ii] = m_Q[0]*Q_midbins[ii] + b_Q[0]    
        elif Q_midbins[ii] > Q_frac_midbins[-1]:    # extrapolate above highest mass bin
            Q_spec_completeness_correction[ii] = m_Q[-1]*Q_midbins[ii] + b_Q[-1]    
        elif Q_midbins[ii] > Q_frac_midbins[0] and Q_midbins[ii] < Q_frac_midbins[-1]:    # interpolate in between all other points
            for jj in range(len(Q_frac_midbins)-1):
                if Q_midbins[ii] > Q_frac_midbins[jj] and Q_midbins[ii] < Q_frac_midbins[jj+1]:
                    Q_spec_completeness_correction[ii] = m_Q[jj]*Q_midbins[ii] + b_Q[jj]
        else:
            print('Error in Q spec completeness correction computation. ABORT')
            break   
#    
if plot_flag == 1:                       # plot interpolated/extrapolated points on top of computed correction fractions
    plt.scatter(SF_midbins,SF_spec_completeness_correction,c='b', marker='+', linewidths = 0)
    plt.scatter(Q_midbins,Q_spec_completeness_correction,c='r', marker='x', linewidths = 0)
    MC.xlim=(7.25,12.5)
#
# apply correction; NOTE: need to divide raw_SMF by the spec_completeness_correction, not multiply, hence taking the inverse
SF_spec_completeness_correction = (1/SF_spec_completeness_correction)   
Q_spec_completeness_correction = (1/Q_spec_completeness_correction)
# compute how many objects are added to each mass bin as a result of applying the spec_completeness_correction to the *_raw_smf lists. confirm that (# added to SF) + (# added to Q) = (# added to total)
SF_spec_completeness_diff = correction_difference(SF_phot_smf,np.transpose(SF_spec_completeness_correction))  
Q_spec_completeness_diff = correction_difference(Q_phot_smf,np.transpose(Q_spec_completeness_correction))
total_spec_completeness_diff = SF_spec_completeness_diff + Q_spec_completeness_diff
#
if diag_flag_5 == 1 and diag_flag_master == 1:
    # Display correction factors
    print('\nSection 3.2: Spectroscopic completeness correction factors by bin (multiplicative): ')
    print('SF: ',str(np.transpose(SF_spec_completeness_correction)))
    print('Q: ',str(np.transpose(Q_spec_completeness_correction)),'\n')
    # Display some data for total, SF, Q: 
    print('Galaxies added due to SPECTROSCOPIC COMPLETENESS correction')
    print('SF: ',str(SF_spec_completeness_diff),'\nor ',str((np.sum(SF_spec_completeness_diff)/np.sum(SF_raw_smf))*100),'%\n')
    print('Q: ',str(Q_spec_completeness_diff),'\nor ',str((np.sum(Q_spec_completeness_diff)/np.sum(Q_raw_smf))*100),'%\n')
    print('Total: ',str(total_spec_completeness_diff),'\n')
else:
    print('\nSection 3.2: Galaxies added due to SPEC COMPLETENESS correction\nSF: %s'%np.sum(SF_spec_completeness_diff),'\nQ: %s'%np.sum(Q_spec_completeness_diff))
    #
#
#
#
#
## SECTION (3.3): NORMALIZE the SMF lists for each cluster by the TOTAL MASS in that cluster (i.e. integral under the SMF - so cluster*_smf x *_midbins); begin by adding the corrections just computed to the raw totals
#
## Add mass & spec corrections to *_raw_smf lists
SF_raw_smf_corrected = SF_raw_smf + SF_mass_completeness_diff + SF_spec_completeness_diff
Q_raw_smf_corrected = Q_raw_smf + Q_mass_completeness_diff + Q_spec_completeness_diff
#
## fix the shape of these arrays to be [6 clusters ,# of midbins in SMF], i.e. an array with 6 cells, each cell containing an array with (#of midbins) data points
SF_raw_smf_corrected = SF_raw_smf_corrected.reshape((6,25))
Q_raw_smf_corrected = Q_raw_smf_corrected.reshape((6,25))
#
N,M = SF_raw_smf_corrected.shape      # store dimensions of raw_smf lists
#
## compute the total mass (i.e. sum((# of galaxies in each bin)*(bin mass))
total_mass = np.array([0]*6,dtype='float32')
for ii in range(N):          # go through clusters 1 at a time
    total_mass[ii] = np.sum((SF_raw_smf_corrected[ii]+Q_raw_smf_corrected[ii])*SF_midbins)
#
## now compute the raw relative mass fraction in each bin, for SF/Q by cluster; that is, create an array (row_1=SF, row_2=Q)
SF_rel_mass = np.empty_like(SF_raw_smf_corrected)
Q_rel_mass = np.empty_like(SF_raw_smf_corrected)
for ii in range(N):          # go through clusters 1 at a time
    for jj in range(M):      # go through each mass bin one at a time
        SF_rel_mass[ii][jj] = ((SF_raw_smf_corrected[ii][jj]*SF_midbins[0][jj]) / total_mass[ii])   # mass in the jj'th bin of cluster ii, divided by total mass of cluster ii
        Q_rel_mass[ii][jj] = (Q_raw_smf_corrected[ii][jj]*SF_midbins[0][jj]) / total_mass[ii] # NOT AN ERROR: Q_midbins are offset by 0.05 for plotting purposes. the true value of the point is store in SF_midbins
#
## compute relative mass fractions of SF/Q by cluster after normalization
mass_fraction_by_cluster = np.array([[0]*6]*2,dtype='float32')
#
mass_fraction_by_cluster[0] = np.sum(SF_rel_mass, axis=1)
mass_fraction_by_cluster[1] = np.sum(Q_rel_mass, axis=1)
#
# Mass in a bin is equal to: (# count in that bin) * (mass value of that bin). The above normalizes by mass, so the arrays "*_normalied_mass" contain the normalized amount of MASS IN EACH MASS BIN. To get the normalized # count, you need to divide (the normalized amount of mass in each bin) by (the mass value of that bin)
## Normalize SMFs by cluster. 
SF_smf_by_cluster = np.empty_like(SF_raw_smf_corrected)     #initialize arrays for FINAL SMF
Q_smf_by_cluster = np.empty_like(Q_raw_smf_corrected)
## normalize each cluster individually, for both SF & Q
#
for ii in range(N):          # go through clusters 1 at a time
    for jj in range(M):      # go through each mass bin (within each cluster) 1 at a time
        SF_smf_by_cluster[ii][jj] = (SF_rel_mass[ii][jj]) / SF_midbins[0][jj]  # normalize each cluster by the TOTAL MASS IN EACH CLUSTER (hence the midbins multiplication); 
        Q_smf_by_cluster[ii][jj] = Q_rel_mass[ii][jj] / SF_midbins[0][jj] 
#
## for ease of future calling, combine all clusters into a sinlge SMF
SF_smf = np.sum(SF_smf_by_cluster,axis=0)
Q_smf = np.sum(Q_smf_by_cluster,axis=0)
total_smf = SF_smf + Q_smf
#
#
## now we check that the mass fractions are still the same
## now compute the normalized amount of mass in each bin, for SF/Q by cluster; that is, create an array (row_1=SF, row_2=Q)
## compute the total mass (i.e. sum((# of galaxies in each bin)*(bin mass))
total_norm_mass = np.array([0]*6,dtype='float32')
for ii in range(N):          # go through clusters 1 at a time
    total_norm_mass[ii] = np.sum((SF_smf_by_cluster[ii]+Q_smf_by_cluster[ii])*SF_midbins)
#
SF_norm_mass = np.empty_like(SF_smf_by_cluster)
Q_norm_mass = np.empty_like(Q_smf_by_cluster)
for ii in range(N):      # go through each mass bin one at a time
    for jj in range(M):      # go through each mass bin (within each cluster) 1 at a time
        SF_norm_mass[ii][jj] = ((SF_smf_by_cluster[ii][jj]*SF_midbins[0][jj]) / total_norm_mass[ii])   # mass in the ii'th bin
        Q_norm_mass[ii][jj] = ((Q_smf_by_cluster[ii][jj]*SF_midbins[0][jj]) / total_norm_mass[ii]) # NOT AN ERROR: Q_midbins are offset by 0.05 for plotting purposes. the true value of the point is store in SF_midbins
#
## compute relative mass fractions of SF/Q by cluster after normalization
mass_fraction_by_cluster_norm = np.array([[0]*6]*2,dtype='float32')
#
mass_fraction_by_cluster_norm[0] = np.sum(SF_norm_mass, axis=1)
mass_fraction_by_cluster_norm[1] = np.sum(Q_norm_mass, axis=1)
#






### MAY NEED TO EDIT: diag_flag_6
# display diagnostics before/after normalization
diag_flag_6 = 0              # 0=off, 1=on
#
if diag_flag_6 == 1 and diag_flag_master == 1:
    # the above should make the area under each of the SF/Q curves equal to 1 (by cluster), so the total area under the curve (sum of all clusters) should be 6. 
    print('\nSection 3.3: Normalization diagnostic:\nMass fractions below are per cluster. The mass in each bin was calculated as: \n(count in mass bin)*(value of mass bin) = mass in that bin')
    print('\nTotal corrected raw samples\nSF: %s'%np.sum(SF_raw_smf_corrected,axis=1),'\nQ: %s'%np.sum(Q_raw_smf_corrected,axis=1))
    
#    print('\nTotal SF fraction by mass in each cluster - raw: \n%s'%     $$$   )
    
    print('\nPre-Normalization:\nTotal SF fraction by mass in each cluster: \n%s'%mass_fraction_by_cluster[0])
    print('Total Q fraction by mass in each cluster: \n%s'%mass_fraction_by_cluster[1])
    print('Total relative mass per cluster (check by summing SF + Q from above): \n%s'%np.sum(mass_fraction_by_cluster,axis=0))
    print('\nPost-Normalization:\nTotal SF fraction by mass in each cluster: \n%s'%mass_fraction_by_cluster_norm[0])
    print('Total Q fraction by mass in each cluster: \n%s'%mass_fraction_by_cluster_norm[1])
    print('Total relative mass per cluster (check by summing SF + Q from above): \n%s'%np.sum(mass_fraction_by_cluster_norm,axis=0))
    #
    ## so we normalized the area under each cluster to 1, then summed all 6 clusters. so the total area under the curve should be 6. CHECK THAT
    SF_area = np.sum(SF_smf*SF_midbins)
    Q_area = np.sum(Q_smf*SF_midbins)
    print('\nCheck on final SMF - total area under curve should be 6\nArea under SF_smf: %s'%SF_area,'\nArea under Q_smf: %s'%Q_area,'\nArea under TOTAL curve: %s'%(SF_area+Q_area))
    print('\nTotal NORMALIZED SF fraction by mass in each cluster: \n%s'%mass_fraction_by_cluster_norm[0])
    print('Total Q fraction by mass in each cluster: \n%s'%mass_fraction_by_cluster_norm[1])
    print('Total relative mass per cluster (check by summing SF + Q from above): \n%s'%np.sum(mass_fraction_by_cluster_norm,axis=0))
else:
    print('\nSection 3.3 NORMALIZATION\nTotal RAW SF: %s'%np.sum(SF_raw_smf_corrected),'\nTotal RAW Q: %s'%np.sum(Q_raw_smf_corrected))
    print('Rel. MASS fractions - pre-normalization:   SF: %s'%(np.sum(mass_fraction_by_cluster[0])/6),'    Q: %s'%(np.sum(mass_fraction_by_cluster[1])/6))
    print('\nNOTE: each cluster was normalized by mass to 1, so total normalized curve has area of 6\nTotal NORMALIZED SF: %s'%np.sum(mass_fraction_by_cluster_norm[0]))    
    print('Total NORMALIZED Q: %s'%np.sum(mass_fraction_by_cluster_norm[1]))
    print('Rel. MASS fractions - post-normalization:   SF: %s'%(np.sum(mass_fraction_by_cluster_norm[0])/6),'    Q: %s'%(np.sum(mass_fraction_by_cluster_norm[1])/6))

    
    
    
    
    
                                                                      
# sum cluster array into a single SMF array for the entire sample
SF_smf = np.sum(SF_smf_by_cluster, axis=0)
Q_smf = np.sum(Q_smf_by_cluster, axis=0)

  





#
#
## combine all clusters into full SMF list
SF_smf = SF_smf1 + SF_smf2 + SF_smf3 + SF_smf4 + SF_smf5 + SF_smf6
Q_smf = Q_smf1 + Q_smf2 + Q_smf3 + Q_smf4 + Q_smf5 + Q_smf6
total_smf = SF_smf + Q_smf
#
#
## Now the total area under the total_smf curve is some arbitrary number (~6.8 for 26 bins). Re-normalize all three curves s.t. the area under the total_smf curve is equal to 1. Recall that by dividing each of SF/Q by the total area under the curve (i.e. np.sum(total_smf)) we preserve the relative fraction of SF-to-Q in each mass bin.
#
#diagnostic to confirm relative fraction of SF-to-Q is preserved
rel_SF_1 = SF_smf / total_smf
rel_Q_1 = Q_smf / total_smf
#
# display diagnostic initial status
print('Pre-normalization')
print('Area under SF_smf: ',str(np.sum(SF_smf)))
print('Area under Q_smf: ',str(np.sum(Q_smf)))
print('Area under total_smf: ',str(np.sum(total_smf)))
#
# re-normalize the curves
SF_smf = SF_smf / np.sum(total_smf)
Q_smf = Q_smf / np.sum(total_smf)
total_smf = total_smf / np.sum(total_smf)
#
# re-calculate relative fraction diagnostic post-normalization
rel_SF_2 = SF_smf / total_smf
rel_Q_2 = Q_smf / total_smf
#
# display final diagnostic status
print('\nPost-normalization')
print('Area under SF_smf: ',str(np.sum(SF_smf)))
print('Area under Q_smf: ',str(np.sum(Q_smf)))
print('Area under total_smf: ',str(np.sum(total_smf)))
diff_SF = rel_SF_1 - rel_SF_2
diff_Q = rel_Q_1 - rel_Q_2
print('\nDifference in relative fractions before & after normalization:')
print('SF: ',str(np.sum(diff_SF)))
print('Q: ',str(np.sum(diff_Q)))
#
#
#
#
#
######### THE FOLLOWING CODE IS DEPRECATED #############
#
### I left it in because it was replaced by the few lines of code above in section (iii).2. This just shows what I've learned in the two years since I wrote the below code, as it has been replaced by a few lines above.
#
#
### Mass corrections
### compute total number of SF/Q galaxies in all 6 clusters, then for individual 
### cluster down to its limiting mass. take the ratio and divide the number count 
### for that galaxy in that mass bin by the ratio to correct for mass comleteness 
### by cluster
#
#
##limiting_mass = [7.5,7.8,8.0,7.5,7.4,7.3]  # limiting mass by cluster as determined by magnitude-to-mass plot
### cluster 1: macs1149; limiting mass 7.5
#gal_count1 = np.array([0,0],dtype = float)  # [SF, Q]
#total_count1 = np.array([0,0],dtype = float)
#counter = 0
#size = len(master_cat)
#while counter < size:
#    if master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==1: #secure member, SF
#        if master_cat[counter]['lmass'] > limiting_mass[0]:
#            total_count1[0] +=1
#            if master_cat[counter]['cluster'] == 1:
#                gal_count1[0] +=1
#    elif master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==2: #secure member, Q
#        if master_cat[counter]['lmass'] > limiting_mass[0]:
#            total_count1[1] +=1
#            if master_cat[counter]['cluster'] == 1:
#                gal_count1[1] +=1
#    counter +=1
#mass_correction1 = gal_count1/total_count1
#SF_smf1[0] = SF_smf1[0]/mass_correction1[0]
#Q_smf1[0] = Q_smf1[0]/mass_correction1[1]
##
#### cluster 2: macs1149; limiting mass 7.8
#gal_count2 = np.array([0,0],dtype = float)  # [SF, Q]
#total_count2 = np.array([0,0],dtype = float)
#counter = 0
#size = len(master_cat)
#while counter < size:
#    if master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==1: #secure member, SF
##        if master_cat[counter]['lmass'] > limiting_mass[1]:
#            total_count2[0] +=1
#            if master_cat[counter]['cluster'] == 2:
#                gal_count2[0] +=1
##    elif master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==2: #secure member, SF
#        if master_cat[counter]['lmass'] > limiting_mass[1]:
#            total_count2[1] +=1
#            if master_cat[counter]['cluster'] == 2:
#                gal_count2[1] +=1
#    counter +=1
#mass_correction2 = gal_count2/total_count2
#SF_smf2[0] = SF_smf2[0]/mass_correction2[0]
#SF_smf2[1] = SF_smf2[0]/mass_correction2[0]
#Q_smf2[0] = Q_smf2[1]/mass_correction2[1]
#Q_smf2[1] = Q_smf2[1]/mass_correction2[1]
##
### cluster 3: macs 0717; limiting mass 8.0
#gal_count3 = np.array([0,0],dtype = float)  # [SF, Q]
#total_count3 = np.array([0,0],dtype = float)
#counter = 0
#size = len(master_cat)
#while counter < size:
#    if master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==1: #secure member, SF
#        if master_cat[counter]['lmass'] > limiting_mass[2]:
#            total_count3[0] +=1
#            if master_cat[counter]['cluster'] == 3:
#                gal_count3[0] +=1
#    elif master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==2: #secure member, Q
#        if master_cat[counter]['lmass'] > limiting_mass[2]:
#            total_count3[1] +=1
#            if master_cat[counter]['cluster'] == 3:
#                gal_count3[1] +=1
#    counter +=1
#mass_correction3 = gal_count3/total_count3
#SF_smf3[0] = SF_smf3[0]/mass_correction3[0]
#SF_smf3[1] = SF_smf3[1]/mass_correction3[1]
#Q_smf3[0] = Q_smf3[1]/mass_correction3[1]
#Q_smf3[1] = Q_smf3[1]/mass_correction3[1]
##
### cluster 4: abell370; limiting mass 7.5
#gal_count4 = np.array([0,0],dtype = float)  # [SF, Q]
#total_count4 = np.array([0,0],dtype = float)
#counter = 0
#size = len(master_cat)
#while counter < size:
#    if master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==1: #secure member, SF
#        if master_cat[counter]['lmass'] > limiting_mass[3]:
#            total_count4[0] +=1
#            if master_cat[counter]['cluster'] == 4:
##                gal_count4[0] +=1
#    elif master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==2: #secure member, SF
#        if master_cat[counter]['lmass'] > limiting_mass[3]:
#            total_count4[1] +=1
#            if master_cat[counter]['cluster'] == 4:
#                gal_count4[1] +=1
#    counter +=1
#mass_correction4 = gal_count4/total_count4
#SF_smf4[0] = SF_smf4[0]/mass_correction4[0]
#Q_smf4[0] = Q_smf4[1]/mass_correction4[1]
##
### cluster 5: abell1063; limiting mass 7.4
#gal_count5 = np.array([0,0],dtype = float)  # [SF, Q]
#total_count5 = np.array([0,0],dtype = float)
#counter = 0
#size = len(master_cat)
#while counter < size:
#    if master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==1: #secure member, SF
##        if master_cat[counter]['lmass'] > limiting_mass[4]:
#            total_count5[0] +=1
#            if master_cat[counter]['cluster'] == 5:
#                gal_count5[0] +=1
#    elif master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==2: #secure member, SF
#        if master_cat[counter]['lmass'] > limiting_mass[4]:
#            total_count5[1] +=1
#            if master_cat[counter]['cluster'] == 5:
#                gal_count5[1] +=1
#    counter +=1
#mass_correction5 = gal_count5/total_count5
#SF_smf5[0] = SF_smf5[0]/mass_correction5[0]
#Q_smf5[0] = Q_smf5[1]/mass_correction5[1]
##
### cluster 6: abell2744; limiting mass 7.3
#gal_count6 = np.array([0,0],dtype = float)  # [SF, Q]
#total_count6 = np.array([0,0],dtype = float)
#counter = 0
#size = len(master_cat)
##while counter < size:
#    if master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==1: #secure member, SF
##        if master_cat[counter]['lmass'] > limiting_mass[5]:
#            total_count6[0] +=1
#            if master_cat[counter]['cluster'] == 6:
#                gal_count6[0] +=1
#    elif master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==2: #secure member, SF
#        if master_cat[counter]['lmass'] > limiting_mass[5]:
#            total_count6[1] +=1
#            if master_cat[counter]['cluster'] == 6:
#                gal_count6[1] +=1
#    counter +=1
#mass_correction6 = gal_count6/total_count6
#SF_smf6[0] = SF_smf6[0]/mass_correction6[0]
#Q_smf6[0] = Q_smf6[1]/mass_correction6[1]
#
#
#
#
#d######## END OF DEPRECATED CODE #############
#
#
#
#
## (iv) error bars & relative fractions
## Poissonian error bars will be added to the spec. sample only; phot sample requires MCMC error to account for cluster membership correction
#
## Method: treat each bin as its own Poisson distribution, w/ expectation equal to count in each respective bin. the error is then the sqrt of that count
#
## *******this is a test code to plot error bars, which will be applied to the entire field sample. adjustment to be made later to isolate spec sample alone
#
#
## Sort the raw SMFs into spec & phot sub-samples. confirm the total number of raw clusters members are accounted for. compute error for spectroscopy (i.e. sqrt(N)). The phot sub-sample will be used later to make the false pos/neg correction. display how many galaxies of each type there are as this information will be put in the paper. 
#
#
#
#
#
#
#
#
#
## field populations
SF_field_smf, SF_field_bins = np.histogram(SF_field_hist, bins=SF_bins,range=range2)
Q_field_smf, Q_field_bins = np.histogram(Q_field_hist, bins=SF_bins,range=range2)      #use same bins for both populations to facilitate comparison
frac_field_smf = np.empty(shape=(2,len(SF_bins)-1),dtype='float64')       #1st row: SF fraction; 2nd row Q fraction
#apply correction factors to field
#SF_field_smf = SF_field_smf/SF_correction   #spec mass correction only applies to cluster, not fli
#Q_field_smf = Q_field_smf/Q_correction
total_field_smf = SF_field_smf + Q_field_smf
#
## I. error values by sample type, cluster & field
SF_error = np.sqrt(SF_smf)        #assign errors by population as sqrt of count in each mass bin
Q_error = np.sqrt(Q_smf)
total_error = np.sqrt(total_smf)
SF_field_error = np.sqrt(SF_field_smf)        #assign errors by population as sqrt of count in each mass bin
Q_field_error = np.sqrt(Q_field_smf)
total_field_error = np.sqrt(total_field_smf)
#

### Progagating errors
### make a correction for the spectroscopic completeness errors: add relative errors in quadrature with # count error defined just above in SF_error, Q_error, 
SF_frac_rel = np.empty_like(SF_frac_error, dtype='float64')       #store relative error of spec correction
Q_frac_rel = np.empty_like(Q_frac_error, dtype='float64')
SF_spec_bin_error = np.empty_like(SF_error, dtype='float64') #interpolate rel error into bins
Q_spec_bin_error = np.empty_like(Q_error, dtype='float64')
#
SF_frac_rel = SF_frac_error/SF_frac
Q_frac_rel = Q_frac_error/Q_frac
# interpolate SF relative error for SMF mass bins
i = 0               #main counting index
j = 1               #secondary index for interpolation b/w correction data points
k = 1
for x in SF_spec_bin_error:
    if i < 2:
        SF_spec_bin_error[i] = SF_frac_rel[0]
        i +=1
    elif i >= 2 and i < 8:
        SF_spec_bin_error[i] = ((5-j)/4)*SF_frac_rel[0] + ((j-1)/4)*SF_frac_rel[1]
        i +=1
        j +=1
    elif i >= 8 and i <= 12:
        SF_spec_bin_error[i] = ((5-k)/4)*SF_frac_rel[1] + ((k-1)/4)*SF_frac_rel[2]
        i +=1
        k +=1
    elif i > 12:
        SF_spec_bin_error[i] = SF_frac_rel[2]
#
# interpolate Q relative error for SMF mass bins
i = 0               #main counting index
j = 1               #secondary indices for interpolation b/w correction data points
k = 1
l = 1
m = 1
for x in Q_spec_bin_error:
    if i < 1:
        Q_spec_bin_error[i] = Q_frac_rel[0]
        i +=1
    elif i >= 1 and i < 4:
        Q_spec_bin_error[i] = ((3-j)/2)*Q_frac_rel[0] + ((j-1)/2)*Q_frac_rel[1]
        i +=1
        j +=1
    elif i >= 4 and i <= 6:
        Q_spec_bin_error[i] = ((3-k)/2)*Q_frac_rel[1] + ((k-1)/2)*Q_frac_rel[2]
        i +=1
        k +=1
    elif i > 6 and i <= 9:
        Q_spec_bin_error[i] = ((3-l)/2)*Q_frac_rel[2] + ((l-1)/2)*Q_frac_rel[3]
        i +=1
        l +=1
    elif i > 9 and i <= 12:
        Q_spec_bin_error[i] = ((3-m)/2)*Q_frac_rel[3] + ((m-1)/2)*Q_frac_rel[4]
        i +=1
        m +=1
    elif i > 12:
        Q_spec_bin_error[i] = Q_frac_rel[4]
#       
###   now compute relative error from SMF, then sum in quadrature
SF_rel = SF_error/SF_smf
Q_rel = Q_error/Q_smf
#SF_error = np.sqrt(SF_rel**2 + SF_spec_bin_error**2)*SF_smf
SF_error = np.sqrt(SF_smf)
#Q_error = np.sqrt(Q_rel**2 + Q_spec_bin_error**2)*Q_smf
Q_error = np.sqrt(Q_smf)
#
#
# initialize fractional error arrays
frac_error = [[0 for x in range(len(SF_smf))] for y in range(2)]
frac_field_error = [[0 for x in range(len(SF_field_smf))] for y in range(2)]
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
        frac_field_error [1][counter] = 0
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
##
#
#
#######################
#######################
## this next step is a bit of a necessary diversion: NORMALIZING THE CURVE
#
#
#
#
total_rel_error = total_error/total_smf
# normalize each array to 1
SF_sum = np.sum(SF_smf)
Q_sum = np.sum(Q_smf)
#total_sum = np.sum(total_smf)
SF_smf = SF_smf / SF_sum
Q_smf = Q_smf / Q_sum
total_smf = SF_smf + Q_smf
# normalize error by preserving the relative error same as before
SF_error = SF_rel*SF_smf
Q_error = Q_rel*Q_smf
total_error = total_rel_error*total_smf
#
# do the same for field pop
SF_field_rel = SF_field_error / SF_field_smf
Q_field_rel = Q_field_error / Q_field_smf
total_field_rel_error = total_field_error/total_field_smf
# normalize each array to 1
SF_field_sum = np.sum(SF_field_smf)
Q_field_sum = np.sum(Q_field_smf)
#total_sum = np.sum(total_smf)
SF_field_smf = SF_field_smf / SF_field_sum
Q_field_smf = Q_field_smf / Q_field_sum
total_field_smf = SF_field_smf + Q_field_smf
# normalize error by preserving the relative error same as before
SF_field_error = SF_field_rel*SF_field_smf
Q_field_error = Q_field_rel*Q_field_smf
total_field_error = total_field_rel_error*total_field_smf
#
#
#
frac_smf = np.empty(shape=(2,len(SF_bins)-1),dtype='float64')       #1st row: SF fraction; 2nd row Q fraction
## calculate bin fractions for lower panels of SMF plot
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
#
#
########################
########################
#
# BREAK
#
########################
########################
#
#
#
#
#
#
#
#    
#
## (vi) fit a SCHECHTER function to scatter plot
#
######## This section has been broken out into its own program, called "emcee_initial" & "emcee_chi2", the latter of which uses chi-squared as the cost function and is the code to be implemented in the final run for this project 
#
## The following summarizes the result of the MCMC simulation and sets up the appropriate arrays for plotting
#
## define x array
## generate points to plot Schechter fit
x = np.linspace(SF_midbins[0],SF_midbins[len(SF_midbins)-1],num=1000)#
#
SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [10.14499598,0.02096564,-1.2720757]     # 100 walkers, 500,000 steps
SFM_star_sigma, SFphi_sigma, SFalpha_sigma = [0.57568334,0.00589964,0.03997277]
#single-schechter Q pop
QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [10.76047713,0.05031532,-0.94764447]     # 100 walkers, 200,000 steps
QM_star_sigma, Qphi_sigma, Qalpha_sigma = [0.04034471,0.00262553,0.01046056]
#double-schechter Q pop
#QM_star_mcmc, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [7.4296233,-15.09718791,-19.94457933,2.81944241,-2.9740897]     # 100 walkers, 200,000 steps
#QM_star_sigma, Qphi1_sigma, Qalpha1_sigma, Qphi2_sigma, Qalpha2_sigma = [0.00010148,0.00010006,0.00010002,0.00010189,0.00010069]
#
TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [10.88441285,0.03862147,-1.17262005]
TM_star_sigma, Tphi_sigma, Talpha_sigma = [0.04220442,0.00251507,0.00934645]
#
## build a model for plotting
SF_model_mcmc = np.log(10)*SFphi_mcmc*(10**((x-SFM_star_mcmc)*(1+SFalpha_mcmc)))*np.exp(-10**(x-SFM_star_mcmc))
#single-schechter Q pop
Q_model_mcmc = np.log(10)*Qphi_mcmc*(10**((x-QM_star_mcmc)*(1+Qalpha_mcmc)))*np.exp(-10**(x-QM_star_mcmc))
#double-schechter Q pop
#Q_model_mcmc = np.log(10)*np.exp(-10**(x-QM_star_mcmc))*((Qphi1_mcmc*(10**((x-QM_star_mcmc)*(1+Qalpha1_mcmc))))+(Qphi2_mcmc*(10**((x-QM_star_mcmc)*(1+Qalpha2_mcmc)))))
T_model_mcmc = np.log(10)*Tphi_mcmc*(10**((x-TM_star_mcmc)*(1+Talpha_mcmc)))*np.exp(-10**(x-TM_star_mcmc))
#
# create arrays for plotting purposes to display values down to a y_min = 1
SF_model_ml_plot = []
SF_model_mcmc_plot = []
x_plot_SF = []
for ii in range(len(SF_model_mcmc)):
    if SF_model_mcmc[ii] > 1e-4:
        SF_model_mcmc_plot.append(SF_model_mcmc[ii])
        x_plot_SF.append(x[ii])
# do same for Q pop
Q_model_ml_plot = []
Q_model_mcmc_plot = []
x_plot_Q = []
for ii in range(len(Q_model_mcmc)):
    if Q_model_mcmc[ii] > 1e-4:
        Q_model_mcmc_plot.append(Q_model_mcmc[ii])
        x_plot_Q.append(x[ii])
# do same for total pop
T_model_ml_plot = []
T_model_mcmc_plot = []
x_plot_T = []
for ii in range(len(Q_model_mcmc)):
    if T_model_mcmc[ii] > 1e-4:
        T_model_mcmc_plot.append(T_model_mcmc[ii])
        x_plot_T.append(x[ii])
#
#
## SECTION (vii): create plot
## upper: SMF for cluster, field
## lower: fractions of SF/Q in cluster, field
### SMF
## TOTAL (all clusters)
plt.close()
SMF = plt.figure(num=1)
gs = gridspec.GridSpec(2,2, wspace=0, hspace=0, width_ratios=[1,1], height_ratios=[2,1])   #make a tiled-plot like vdB2013 w/ fractions below, this line sets the proporitons of plots in the figure
#gs = gridspec.GridSpec(2,3, width_ratios=[1,1,1], height_ratios=[2,1])   #make a tiled-plot like vdB2013 w/ fractions below, this line sets the proporitons of plots in the figure
# Cluster
cluster = plt.subplot(gs[0])      
xa = plt.errorbar(SF_midbins,SF_smf,yerr=SF_error, fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)#yerr=SF_error,
xb = plt.errorbar(Q_midbins,Q_smf,yerr=Q_error,fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
xc = plt.errorbar(SF_midbins,total_smf,yerr=total_error,fmt='.k',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
#Plot Schechter fits:
#plt.plot(x_plot_Q,Q_model_ml_plot, ':r')
plt.plot(x_plot_Q,Q_model_mcmc_plot, '--r')
#plt.plot(x_plot_SF,SF_model_ml_plot, ':c', label = 'Max. Likelihood', linewidth = 0.5)
plt.plot(x_plot_SF,SF_model_mcmc_plot, '--b', label = 'MCMC', linewidth = 0.5)
plt.plot(x_plot_T,T_model_mcmc_plot, 'k')

#plt.plot(SF_midbins,Q_smf,'.r',linewidth=0.5)#Qx_new, Qy_new,'-r',linewidth=0.5)
#plt.plot(SF_midbins,total_smf,'.k',linewidth=0.5)#totalx_new, totaly_new,'-k',linewidth=0.5)
#plt.errorbar(SF_midbins,SF_smf,fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)#yerr=SF_error,
plt.xscale('linear')
#plt.xlabel('log(M/M_sol)')
plt.xlim=(8,12.25)
plt.yscale('log')
plt.ylim=(1,1001)
cluster.minorticks_on()
cluster.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off',labelbottom='off')
cluster.yaxis.set_label_position("left")
plt.ylabel('# count')
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
#
#Plot Schechter fits:
#plt.plot(x,FQ_model, 'r')
#plt.plot(x,FT_model, 'k')
###plt.plot(x,FSF_model, 'b')
#
plt.xscale('linear')
#plt.xlabel('log(M/M_sol)')
plt.xlim=(8,12.25)
plt.yscale('log')
plt.ylabel('# count')
plt.ylim=(1,1001)
field.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=True,labelleft=False,labelbottom=False)
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
########
########
#
#
#
#
#
#
########
########
# Plot: compare environments by populatin (i.e. plot Cluster vs Field for Total, SF, & Q)
plt.close()
SMF = plt.figure(num=1)
gs = gridspec.GridSpec(1,3, wspace=0, hspace=0, width_ratios=[1,1,1])   #make a tiled-plot 
# Total population
cluster = plt.subplot(gs[0])      
xa = plt.errorbar(SF_midbins,total_smf,yerr=total_error,fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
xb = plt.errorbar(SF_midbins,total_field_smf,yerr=total_field_error, fmt='.g',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
#Plot Schechter fits:
plt.plot(x_plot_T,T_model_mcmc_plot, '--b')
plt.xscale('linear')
plt.xlabel('log(M/$M_{\odot}$)')
plt.xlim=(8,12.25)
plt.yscale('log')
plt.ylim=(1,1001)
cluster.minorticks_on()
cluster.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=False,labelbottom=True)
cluster.yaxis.set_label_position("left")
plt.ylabel('???')
plt.title('Total')
plt.legend((xa,xb),('Cluster','Field'),scatterpoints=1,loc='lower left', frameon=False, fontsize = 'x-small')
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')
#
# SF population
cluster = plt.subplot(gs[1])      
xa = plt.errorbar(SF_midbins,SF_smf,yerr=SF_error,fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
xb = plt.errorbar(SF_midbins,SF_field_smf,yerr=SF_field_error, fmt='.g',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
#Plot Schechter fits:
plt.plot(x_plot_SF,SF_model_mcmc_plot, '--b', label = 'MCMC', linewidth = 0.5)
plt.xscale('linear')
plt.xlabel('log(M/$M_{\odot}$)')
plt.xlim=(8,12.25)
plt.yscale('log')
plt.ylim=(1,1001)
cluster.minorticks_on()
cluster.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelleft=False,labelright=False,labelbottom=True)
#cluster.yaxis.set_label_position("left")
#plt.ylabel('???')
plt.title('Star-Forming')
#plt.legend((xa,xb),('Cluster','Field'),scatterpoints=1,loc='lower left', frameon=False, fontsize = 'x-small')
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')
#
# Q population
cluster = plt.subplot(gs[2])      
xa = plt.errorbar(SF_midbins,Q_smf,yerr=Q_error,fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
xb = plt.errorbar(SF_midbins,Q_field_smf,yerr=Q_field_error, fmt='.g',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
#Plot Schechter fits:
plt.plot(x_plot_Q,Q_model_mcmc_plot, '--b', label = 'MCMC', linewidth = 0.5)
plt.xscale('linear')
plt.xlabel('log(M/$M_{\odot}$)')
plt.xlim=(8,12.25)
plt.yscale('log')
plt.ylim=(1,1001)
cluster.minorticks_on()
cluster.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelleft=False,labelright=True,labelbottom=True)
cluster.yaxis.set_label_position("right")
plt.ylabel('???')
plt.title('Quiescent')
#plt.legend((xa,xb),('Cluster','Field'),scatterpoints=1,loc='lower left', frameon=False, fontsize = 'x-small')
plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')
#
#
#
#
#
#
#    
#################
#################  END  #################
#################
#
#
#
#
#
#

