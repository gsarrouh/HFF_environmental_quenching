# Created on Fri Jun 26 13:17:46 2020
#
################## spec_completeness_binning.py ##################
#
#
### WHAT THIS PROGRAM DOES:
### This program used to be SECTION (3.2) from "master_smfz*.py". Now it stands alone, and is CALLED BY "master_data*.py". For a given spectroscopic & photometric redshift cut, a population of false positives & false negatives results. The following program uses 3 different binning techniques to try various numbers of bins, evaluating the resulting spectroscopic correction factors, and computing a "mertic of merit", by which we may judge various combinations of redshift cuts, bins, and binning techniques.
#
## This is the old description of SECTION (3.2):
### SECTION (3.2): calculate SPECTROSCOPIC COMPLETENESS correction. basically, look at all the false positives/false negatives, and sort them by type (i.e. SF/Q). then bin them (i.e. make histograms of false pos/neg for each of SF/Q). take their ratio of false pos to false neg, and plot that ratio. this is the correction factor to be applied to the photometric subsample
#
## SPEC. BINNING: iterate through different number of histogram bins to see which yields a set of corrections closest in general to ~1; this has now been expanded to test: 
#
## method = 1, SYMMETRIC binning; 
## method = 2, ASYMMETRIC binning, EQUAL NUMBER of objects in each bin; 
## method = 3, ASYMMETRIC binning, EQUAL MASS in each bin; 
## method = 4; ASYMMETRIC binning using built-in "np.histogram_bin_edges" function, using black-box binning method "auto", which is the Maximum of the "STURGES" & "FREEDMAN DIACONIS" estimators
#
#
#
### Section summary:
#
### PROGRAM START
#
### EVALUATE DIFFERENT BIN #s & BINNING TECHNIQUES:
### (0)    modules & definitions
### (1)    sort data into false pos/neg lists for SF/Q;
### (1.1)   add Summary Table comparing to false pos/neg found in "master_data*.py"
### (2)    Method 1: SYMMETRIC bins; 
### (3)    Method 2: ASYMMERTIC bins, EQUAL NUMBER distribution among bins;
### (4)    Method 3: ASYMMERTIC bins, EQUAL MASS distribution among bins;
### (5)    compute BIN EDGES and METRIC OF MERIT;
### (6)    prepare strings & to WRITE to file;
#
### PROGRAM END
#
#
#
### NOTE: To find a section, search "SECTION (*)" in all caps, where * is the section number you want. To find fields which require user input, search "MAY NEED TO EDIT" in all caps. Some of these fields may direct you to manually enter input into a sub-program. 
#
#
#
###################     PROGRAM START
#
#
## SECTION (0): modules & definitions
#
## MODULES
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import astropy
from astropy.table import Table
from astropy.table import Column
import time
#from master_data_7_final import master_cat
#
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
#
#
## TIME_FLAG: START
## superior time_flag which supercedes all others and times the entire program
time_flag = 1     # track & print time to execute current section
#
if time_flag == 1:
    if project_time_flag == 1:
        pass
    else:
        start_time = time.time()
#
#
#
## MASTER DIAGNOSTIC FLAG: this will enable (1=on) or suppress (0=off) the diagnostic output of this file, for de-bugging purposes
diag_flag = 1
#
## Section summary table flags
summary_flag_1 = 1                   # initial loop to pick up all the false pos/neg
#
#
## SECTION (1): setup LISTS, i.e. compile lists of SF/Q false pos/neg populations; define "midbins" function
# 
## ignore 'divide by zero' errors, since they are handled with exceptions below
np.seterr(divide='ignore')
#
## define some key information about the study: limiting mass, mass range, SMF bin midpoints

if 'limiting_mass' in locals():
    pass
else:
    limiting_mass = [7.14,7.63,8.2,7.5,7.34,7.07]        #NOTE: limiting masses are z_cutoff dependent
##### OLD LIMITING MASSES   limiting_mass = [7.5,7.8,8.0,7.5,7.4,7.3] # clusters 1,2,3,4,5,6, see IDs below
range2 = [min(limiting_mass),12.3]
cluster_names = ['M0416','M1149','M0717','A370','A1063','A2744']   # in the order corresponding to limiting mags above
## now setup of the bin midpoints for the SMF. no need to actually create the SMF, just define an array with the same bin midpoints
#
## MAY NEED TO EDIT!!! -- UPDATE: bin_width now set in main_project_file.py
## if you change the bin width in the "master_smf*.py" code, YOU MUST ALSO CHANGE IT HERE. just search search "MAY NEED TO EDIT" in that file.
num_points = int((round((range2[1]-range2[0])/bin_width))+1)       # compute # of data points
num_bins = np.linspace(range2[0],range2[1],num_points)
#
SF_midbins = midbins(num_bins)                                # define 
#
SF_mem = []
SF_pos = []
SF_neg = []
Q_mem = []
Q_pos = []
Q_neg = []
mem_by_cluster = np.array([[0]*6]*2)
pos_by_cluster = np.array([[0]*6]*2)    #for tracking false pos/neg by cluster; row_1=SF, row_2=Q
neg_by_cluster = np.array([[0]*6]*2)
mem_below_lim_mass = np.array([[0]*6]*2)
objects_below_lim_mass = np.array([[0]*6]*2)    # for tracking objects below the limiting mass of each cluster
#
for counter in range(len(master_cat)):
    for ii in range(len(limiting_mass)):
        if master_cat['cluster'][counter] == (ii+1):           # only look at objects above the limiting mass for each cluster
            if master_cat['lmass'][counter] > limiting_mass[ii]:      
                if master_cat['type'][counter] == 1:      # type=1 for SF
                    if master_cat['member'][counter] == 0:     # member=0 for cluster member
                        SF_mem.append(master_cat['lmass'][counter])
                        for ii in range(len(mem_by_cluster[0])):
                            if master_cat['cluster'][counter] == (ii+1):
                                mem_by_cluster[0][ii]+=1
                    elif master_cat['member'][counter] == 2:     # member=2 for false pos
                        SF_pos.append(master_cat['lmass'][counter])
                        for ii in range(len(pos_by_cluster[0])):
                            if master_cat['cluster'][counter] == (ii+1):
                                pos_by_cluster[0][ii]+=1           # track false pos for SF
                    elif master_cat['member'][counter] == 3:   # member=3 for false neg
                        SF_neg.append(master_cat['lmass'][counter])
                        for ii in range(len(pos_by_cluster[0])):
                            if master_cat['cluster'][counter] == (ii+1):
                                neg_by_cluster[0][ii]+=1           # track false neg for SF
                elif master_cat['type'][counter] == 2:     # type=2 for Q
                    if master_cat['member'][counter] == 0:     # member=0 for cluster member
                        Q_mem.append(master_cat['lmass'][counter])
                        for ii in range(len(mem_by_cluster[1])):
                            if master_cat['cluster'][counter] == (ii+1):
                                mem_by_cluster[1][ii]+=1
                    elif master_cat['member'][counter] == 2:     # member=2 for false pos
                        Q_pos.append(master_cat['lmass'][counter])
                        for ii in range(len(pos_by_cluster[1])):
                            if master_cat['cluster'][counter] == (ii+1):
                                pos_by_cluster[1][ii]+=1           # track false pos for Q
                    elif master_cat['member'][counter] == 3:   # member=3 for false neg
                        Q_neg.append(master_cat['lmass'][counter])
                        for ii in range(len(neg_by_cluster[1])):
                            if master_cat['cluster'][counter] == (ii+1):
                                neg_by_cluster[1][ii]+=1           # track false neg for Q
            else: 
                if master_cat['member'][counter] == 0:
                    if master_cat['type'][counter] == 1:           # SF  false pos/neg below limiting mass
                        for ii in range(len(mem_below_lim_mass[0])):
                            if master_cat['cluster'][counter] == (ii+1):
                                mem_below_lim_mass[0][ii]+=1
                    elif master_cat['type'][counter] == 2:         # Q  false pos/neg below limiting mass
                        for ii in range(len(mem_below_lim_mass[1])):
                            if master_cat['cluster'][counter] == (ii+1):
                                mem_below_lim_mass[1][ii]+=1
                elif master_cat['member'][counter] == 2 or master_cat['member'][counter] == 3:
                    if master_cat['type'][counter] == 1:           # SF  false pos/neg below limiting mass
                        for ii in range(len(objects_below_lim_mass[0])):
                            if master_cat['cluster'][counter] == (ii+1):
                                objects_below_lim_mass[0][ii]+=1
                    elif master_cat['type'][counter] == 2:         # Q  false pos/neg below limiting mass
                        for ii in range(len(objects_below_lim_mass[1])):
                            if master_cat['cluster'][counter] == (ii+1):
                                objects_below_lim_mass[1][ii]+=1
#
#
## SECTION (1.1) - Summary Table 1: Did we pick up all the false pos/neg as reported in "master_data*.py"?
#
if summary_flag_1 == 1 or adams_flag == 1:
    ## Summarize initial data stats in table
    pos_neg_names = Column(['TOTAL Members','TOTAL False Pos.','TOTAL False Neg.','SF - Members','SF - False Pos.','SF - False Neg.','SF (mem)- LBLM','SF (pos/neg)- LBLM','Q - Members','Q - False Pos.','Q - False Neg.','Q (mem) - LBLM','Q (pos/neg) - LBLM','SUM Member','SUM False Pos./Neg.'],name='Property')
    col_names = cluster_names
    # SF table
    pos_neg0 = Column([np.sum([mem_phot,mem_spec]),np.sum(pos_spec),np.sum(neg_spec),np.sum(mem_by_cluster[0]),np.sum(pos_by_cluster[0]),np.sum(neg_by_cluster[0]),np.sum(mem_below_lim_mass[0]),np.sum(objects_below_lim_mass[0]),np.sum(mem_by_cluster[1]),np.sum(pos_by_cluster[1]),np.sum(neg_by_cluster[1]),np.sum(mem_below_lim_mass[1]),np.sum(objects_below_lim_mass[1]),np.sum([mem_by_cluster,mem_below_lim_mass]),np.sum([pos_by_cluster,neg_by_cluster,objects_below_lim_mass])],name='Total')  # total column
    pos_neg_stats = Table([pos_neg_names,pos_neg0])
    for ii in range(len(mem_spec[0])):
        col = Column([np.sum([mem_phot[0][ii],mem_spec[0][ii],mem_phot[1][ii],mem_spec[1][ii]]),np.sum([pos_spec[0][ii],pos_spec[1][ii]]),np.sum([neg_spec[0][ii],neg_spec[1][ii]]),mem_by_cluster[0][ii],pos_by_cluster[0][ii],neg_by_cluster[0][ii],mem_below_lim_mass[0][ii],objects_below_lim_mass[0][ii],mem_by_cluster[1][ii],pos_by_cluster[1][ii],neg_by_cluster[1][ii],mem_below_lim_mass[1][ii],objects_below_lim_mass[1][ii],np.sum([mem_by_cluster[0][ii],mem_by_cluster[1][ii],mem_below_lim_mass[0][ii],mem_below_lim_mass[1][ii]]),np.sum([pos_by_cluster[0][ii],neg_by_cluster[0][ii],objects_below_lim_mass[0][ii],pos_by_cluster[1][ii],neg_by_cluster[1][ii],objects_below_lim_mass[1][ii]])],name=col_names[ii])
        pos_neg_stats.add_column(col)  # add columns to table one cluster at a time
    #
    print('\nSummary Table: False Pos./Neg.\n%s'%pos_neg_stats)
    print('NOTE: TOTALs reported in first three rows are from Summary Tables 4/6 in "master_data*.py".\nNOTE: LBLM = galaxies Lost Below Limiting Mass for their cluster.\n')
#
#
## sort false pos/neg lists in ascending order
SF_pos = np.sort(SF_pos)
SF_neg = np.sort(SF_neg)
Q_pos = np.sort(Q_pos)
Q_neg = np.sort(Q_neg)
#
#
#
## open a file to print to
f = open('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/diagnostic_outputs/spec_binning/z_spec_%.3f'%z_cutoff[0]+'/binning_%.3f'%z_cutoff[0]+'_spec_cutoff_%.3f'%z_cutoff[1]+'_phot_cutoff.txt','w+')
#
#
## to be used in building strings throughout program
delim = ','   
## write a header for the file, start with hashtag to identify comment
header1 = 'z_spec_cutoff'+delim+'z_phot_cutoff'+delim+'type'+delim+'[limiting_mass]'+delim+'MoM'+delim+'TOTAL_MoM'+delim+'cluster_members'+delim+'false pos.'+delim+'false neg.'+delim+'binning_method'+delim+'#_of_bins'+delim+'[bin_edges]\n'
#
f.write(header1)
#
#
#########  Beginning of binning methods
#
### SECTION (2): METHOD 1
### METHOD 1: bin SF & Q in SYMMETRIC BINS, then compute false pos/neg fractions by mass bin for correction factors. 
#
## number of histogram mass bins to try for spec. completeness correction
num_bins_to_try = np.arange(3,6,1)   # np.arange([...]) = [array_start,array_end_plus_increment,increment]
# write a loop that interatively uses a different number of bins in the histrogram, to isolate the largest number for which all bins have at least one entry; NOTE: the lines that stop the loop have been commented out, to investigate the relative fraction of false pos/neg for each different # of bins
#
# 
method = 1
for number in range(len(num_bins_to_try)):
    #
    # make histograms
    SF_hist, bins_SF = np.histogram(SF_mem, bins=num_bins_to_try[number], range=range2)
    SF_pos_hist, bins_SF = np.histogram(SF_pos, bins=num_bins_to_try[number], range=range2)
    SF_neg_hist, bins_SF = np.histogram(SF_neg, bins=num_bins_to_try[number], range=range2)
    Q_hist, bins_Q = np.histogram(Q_mem, bins=num_bins_to_try[number], range=range2)
    Q_pos_hist, bins_Q = np.histogram(Q_pos, bins=num_bins_to_try[number], range=range2)
    Q_neg_hist, bins_Q = np.histogram(Q_neg, bins=num_bins_to_try[number], range=range2)
    #
    bins_SF = np.round(bins_SF,decimals=3)
    bins_Q = np.round(bins_Q,decimals=3)
    #
    ## compute (member + false pos) / (member + false neg) RATIO
    SF_ratio = np.round(((SF_hist + SF_pos_hist)/(SF_hist + SF_neg_hist)),decimals=3)
    Q_ratio = np.round(((Q_hist + Q_pos_hist)/(Q_hist + Q_neg_hist)),decimals=3)
    #
    # compute midbins for spec. mass completeness plot (i.e. plot of false pos/false neg ratios)
    SF_ratio_midbins = midbins(bins_SF)
    Q_ratio_midbins = midbins(bins_Q)
    #
    #
    #    
    ## call a new file: "spec_correction_factors.py" to interpolate/extrapolate the correction factors to the SMF, and return the metric of merit
    #
    exec(open('spec_correction_factors.py').read())      #opens and executes the script
    #
    #
    ## compute variance of SF/Q ratios from 1
    if np.isnan(SF_metric) == 1:
        SF_var = float('NaN')
    else:
        SF_var = SF_metric
    if np.isnan(Q_metric) == 1:
        Q_var = float('NaN')
    else:
        Q_var = Q_metric
    if np.sum([np.isnan(SF_var),np.isnan(Q_var)]) == 0:
        tot_var = SF_var + Q_var
    else:
        tot_var = float('NaN')
    #
    if diag_flag == 1:
        print('#bins: %s'%num_bins_to_try[number],'; method: %s'%method,'; spec cut: %s'%z_cutoff[0],'; phot cut: %s'%z_cutoff[1])
        print('SF metric: %s'%SF_var,';  Q metric: %s'%Q_var)
    #
    #if np.isnan(SF_var) == 1 or np.isnan(Q_var) == 1:
    #    total_var = float('NaN')
    #else: total_var = np.round((SF_var+Q_var),decimals=3)
    #
    #
    ## SECIONT (6): prepare what to WRITE to file - (this is copied&pasted from below for formatting consistency)
    bin_entry1 = str(np.round(z_cutoff[0],decimals=3))+delim+str(np.round(z_cutoff[1],decimals=3))+delim+'SF'+delim+str(limiting_mass)+delim+'%.5f'%SF_var+delim+str(np.sum([mem_spec[0],mem_phot[0]]))+delim+'%s'%len(SF_pos)+delim+'%s'%len(SF_neg)+delim+str(method)+delim+str(num_bins_to_try[number])+delim+str(bins_SF)
    bin_entry2 = str(np.round(z_cutoff[0],decimals=3))+delim+str(np.round(z_cutoff[1],decimals=3))+delim+'Q'+delim+str(limiting_mass)+delim+'%.5f'%Q_var+delim+str(np.sum([mem_spec[1],mem_phot[1]]))+delim+'%s'%len(Q_pos)+delim+'%s'%len(Q_neg)+delim+str(method)+delim+str(num_bins_to_try[number])+delim+str(bins_Q)
    writer = '%s'%str(bin_entry1)+'\n'+'%s'%str(bin_entry2)+'\n'
    f.write(writer)
    #
#
#
## ASYMMETRIC BINNING START
#
method_designations = [2]#[2,3]
#
##
for m in range(len(method_designations)):
    method = method_designations[m]
        #
        #
        #
### SECTION (3): METHOD 2
### METHOD 2: ASYMMETRIC BINS for bins with ~EQUAL NUMBER OF OBJECTS in each bin (for each of the SF/Q false pos/neg pairs), return the asymmetric bin edges and compute the midpoints b/w the *_pos/*_neg bin edges. Use these midpoints as the bin edges for a new histogram for the *_pos/*_neg lists, and re-compute the false pos/neg ratio for each of SF/Q. Then print relevant data to a file.
        #
    if method == 2:
        #    
        for number in range(len(num_bins_to_try)):    #number of bin EDGES (i.e. # of bins + 1)
            ## compute the index corresponding to the first evenly-space bin edge
            SF_pos_index = int(np.ceil(len(SF_pos))/(num_bins_to_try[number]))
            SF_neg_index = int(np.ceil(len(SF_neg))/(num_bins_to_try[number]))
            Q_pos_index = int(np.ceil(len(Q_pos))/(num_bins_to_try[number]))
            Q_neg_index = int(np.ceil(len(Q_neg))/(num_bins_to_try[number]))
            ## build arrays over range [7.3,12.3] with "num_bins_to_try[number]" bins, where each bin has ~equal # of objects in it; then compute the midpoints between the false pos bin edges and the false neg bin edges
            num_bins_SF = [[],[]]     # [pos, neg]  i.e. bin edges of mass bins
            num_bins_Q = [[],[]]     # [pos, neg]
            #a=0
            #b=0
            for ii in range(num_bins_to_try[number]):   # -1 b/c there are 1 fewer bins than bin edges
                num_bins_SF[0].append(SF_pos[(ii*SF_pos_index)])
                num_bins_SF[1].append(SF_neg[(ii*SF_neg_index)])
                num_bins_Q[0].append(Q_pos[(ii*Q_pos_index)])
                num_bins_Q[1].append(Q_neg[(ii*Q_neg_index)])
            #
            ## fix the upper limit to be the highest-mass object in the list
            num_bins_SF[0].append(SF_pos[-1])
            num_bins_SF[1].append(SF_neg[-1])
            num_bins_Q[0].append(Q_pos[-1])
            num_bins_Q[1].append(Q_neg[-1])
            #
            ## check that the number of bins is correct and that you didn't add an extra element at the end of the array in this last step (could happen if the length of the SF_pos list is exactly divisible by SF_pos_index, for example)
            #
            for ii in range(2):
                while len(num_bins_SF[ii]) > (num_bins_to_try[number]+1):
                    del num_bins_SF[ii][-2]
                while len(num_bins_Q[ii]) > (num_bins_to_try[number]+1):
                    del num_bins_Q[ii][-2]
            #
            num_bins_SF = np.array(num_bins_SF)
            num_bins_Q = np.array(num_bins_Q)
            #print('# bins: %s'%(num_bins_to_try[number]))
            #print('length of num_bins_* array (which stores bin edges): %s'%len(num_bins_SF[0]))
            #print(('num_bins_SF: %s'%num_bins_SF)
            #
            ## set first/last entry as limits of mass range for smf
            num_bins_SF[0][0] = range2[0]         # reminder: row [0] = false pos; row [1] = false neg
            num_bins_SF[0][-1] = range2[-1]
            num_bins_SF[1][0] = range2[0] 
            num_bins_SF[1][-1] = range2[-1]
            num_bins_Q[0][0] = range2[0]         # reminder: row [0] = false pos; row [1] = false neg
            num_bins_Q[0][-1] = range2[-1]
            num_bins_Q[1][0] = range2[0] 
            num_bins_Q[1][-1] = range2[-1]    #
            ## compute midpoints between false pos/neg bin edges as starting point of search for ideal bin edges
            bin_edge_means_SF = np.mean(num_bins_SF,axis=0)
            bin_edge_means_Q = np.mean(num_bins_Q,axis=0)
            #
            #
            #
            ## Call the file that now stores Section (5) - compute midbins, correction factors, & metric of merit
            #
            #
            exec(open('spec_asymmetric_binning_metric.py').read())
            #
            ######
            #####if np.sum([np.isnan(SF_var),np.isnan(Q_var)]) == 0:
            #####    tot_var = SF_var + Q_var
            #####else:
            #####    tot_var = float('NaN')
            ######
            #
            ## SECIONT (3.1): prepare what to WRITE to file
            bin_entry1 = str(np.round(z_cutoff[0],decimals=3))+delim+str(np.round(z_cutoff[1],decimals=3))+delim+'SF'+delim+str(limiting_mass)+delim+'%.5f'%SF_var+delim+'%.5f'%tot_var+delim+str(np.sum([mem_spec[0],mem_phot[0]]))+delim+'%s'%len(SF_pos)+delim+'%s'%len(SF_neg)+delim+str(method)+delim+str(num_bins_to_try[number])+delim+str(bins_SF)
            bin_entry2 = str(np.round(z_cutoff[0],decimals=3))+delim+str(np.round(z_cutoff[1],decimals=3))+delim+' Q'+delim+str(limiting_mass)+delim+'%.5f'%Q_var+delim+'%.5f'%tot_var+delim+str(np.sum([mem_spec[1],mem_phot[1]]))+delim+'%s'%len(Q_pos)+delim+'%s'%len(Q_neg)+delim+str(method)+delim+str(num_bins_to_try[number])+delim+str(bins_Q)
            writer = '%s'%str(bin_entry1)+'\n'+'%s'%str(bin_entry2)+'\n'
            f.write(writer)
            #
    #      
    elif method == 3:
        #
        #
        #
        #
        #
### SECTION (4): METHOD 3
### METHOD 3: ASYMMETRIC BINS for bins with ~EQUAL AMOUNTS OF MASS in each bin (for each of the SF/Q false pos/neg pairs), return the asymmetric bin edges and compute the midpoints b/w the *_pos/*_neg bin edges. Use these midpoints as the bin edges for a new histogram for the *_pos/*_neg lists, and re-compute the false pos/neg ratio for each of SF/Q. Then print relevant data to a file.
#
        #
        for number in range(len(num_bins_to_try)):
            SF_pos_sum_index = np.ceil(np.sum(SF_pos)/num_bins_to_try[number]) 
            SF_neg_sum_index = np.ceil(np.sum(SF_neg)/num_bins_to_try[number])
            Q_pos_sum_index = np.ceil(np.sum(Q_pos)/num_bins_to_try[number])
            Q_neg_sum_index = np.ceil(np.sum(Q_neg)/num_bins_to_try[number])
            ## build arrays over range [7.3,12.3] with num_bins_to_try[number] bins, where each bin has ~equal amounts of mass in it; then compute the midpoints between the false pos bin edges and the false neg bin edges
            num_bins_SF_index = [[],[]]     # [pos, neg]
            num_bins_Q_index = [[],[]]     # [pos, neg]
            # look through SF/Q false pos/neg lists one galaxy at a time, summing their mass until you reach the amount stored in *_sum_index (e.g. SF_pos_sum_index for SF_pos list), then record that index as a bin edge
            #num_bins_SF[0].append(SF_pos[0])
            ## SF_posfor ii in range(num_bins_to_try[number]):
            for ii in range(num_bins_to_try[number]):
                mass_sum = 0
                for jj in range(len(SF_pos)):
                    mass_sum = mass_sum + SF_pos[jj]
                    if mass_sum >= ii*SF_pos_sum_index:
                        #print('\nSF_pos\nmass sum : %s'%mass_sum,'\nmass limit index equals: %s'%(ii*SF_pos_sum_index),'jj: %s'%jj)
                        num_bins_SF_index[0].append(jj)
                        break
            ## SF_neg
            for ii in range(num_bins_to_try[number]):
                mass_sum = 0
                for jj in range(len(SF_neg)):
                    mass_sum = mass_sum + SF_neg[jj]
                    if mass_sum >= ii*SF_neg_sum_index:
                        #print('\nSF_neg\nmass sum : %s'%mass_sum,'\nmass limit index equals: %s'%(ii*SF_neg_sum_index),'jj: %s'%jj)
                        num_bins_SF_index[1].append(jj)
                        break
            ## Q_pos
            for ii in range(num_bins_to_try[number]):
                mass_sum = 0
                for jj in range(len(Q_pos)):
                    mass_sum = mass_sum + Q_pos[jj]
                    if mass_sum >= ii*Q_pos_sum_index:
                        #print('\nQ_pos\nmass sum : %s'%mass_sum,'\nmass limit index equals: %s'%(ii*Q_pos_sum_index),'jj: %s'%jj)
                        num_bins_Q_index[0].append(jj)
                        break
            ## Q_neg
            for ii in range(num_bins_to_try[number]):
                mass_sum = 0
                for jj in range(len(Q_neg)):
                    mass_sum = mass_sum + Q_neg[jj]
                    if mass_sum >= ii*Q_neg_sum_index:
                        #print('\nQ_neg\nmass sum : %s'%mass_sum,'\nmass limit index equals: %s'%(ii*Q_neg_sum_index),'jj: %s'%jj)
                        num_bins_Q_index[1].append(jj)
                        break
            #
            ## add last (highest-mass) bin edge
            num_bins_SF_index[0].append(len(SF_pos)-1)
            num_bins_SF_index[1].append(len(SF_neg)-1)
            num_bins_Q_index[0].append(len(Q_pos)-1)
            num_bins_Q_index[1].append(len(Q_neg)-1)
            #
            ## check that the length of the false pos list matches the length of the false neg list before converting the lists to arrays (i.e. check that the length matches the # of bin edges). if one is short, append a value of [0] to its list. this is just a patch to make the code run. introducing a zero at the end of the index list will mean the lowest mass is both the first and last bin edge. the algorithm will penalize this selection when computing the metric of merit, assigning it a value of "NaN", thus it will not contaminate our analysis
            #
            ## record length of each list
            len_SF_pos = len(num_bins_SF_index[0])
            len_SF_neg = len(num_bins_SF_index[1])
            len_Q_pos = len(num_bins_Q_index[0])
            len_Q_neg = len(num_bins_Q_index[1])
            #
            ## now compare to number of bin edges (or # of bins + 1, per below. same difference)
            ## check SF false pos
            while len_SF_pos < (num_bins_to_try[number]+1): 
                num_bins_SF_index[0].append(0)
                len_SF_pos = len(num_bins_SF_index[0])
            while len_SF_pos > (num_bins_to_try[number]+1): 
                del num_bins_SF_index[0][-1]
                len_SF_pos = len(num_bins_SF_index[0])
            ## check SF false neg
            while len_SF_neg < (num_bins_to_try[number]+1): 
                num_bins_SF_index[1].append(0)
                len_SF_neg = len(num_bins_SF_index[1])
            while len_SF_neg > (num_bins_to_try[number]+1): 
                del num_bins_SF_index[1][-1]
                len_SF_neg = len(num_bins_SF_index[1])
            ## check Q false pos    
            while len_Q_pos < (num_bins_to_try[number]+1): 
                num_bins_Q_index[0].append(0)
                len_Q_pos = len(num_bins_Q_index[0])
            while len_Q_pos > (num_bins_to_try[number]+1): 
                del num_bins_Q_index[0][-1]
                len_Q_pos = len(num_bins_Q_index[0])
            ## check Q false neg    
            while len_Q_neg < (num_bins_to_try[number]+1): 
                num_bins_Q_index[1].append(0)
                len_Q_neg = len(num_bins_Q_index[1])
            while len_Q_neg > (num_bins_to_try[number]+1): 
                del num_bins_Q_index[1][-1]
                len_Q_neg = len(num_bins_Q_index[1])
            #
            ## convert to arrays
            num_bins_SF_index = np.array(num_bins_SF_index)
            num_bins_Q_index = np.array(num_bins_Q_index)
            #
            ## enable diagnostic output
            if diag_flag == 1 or project_diagnostic_flag == 1:
                if project_diagnostic_flag == 0:
                    pass
                else:
                    print('# bins: %s'%(num_bins_to_try[number]))
                    print('length of num_bins_* array: %s'%len(num_bins_SF_index[0]),'\n# of SF false pos.: %s'%len(SF_pos),'\n# of SF false neg.: %s'%len(SF_neg),'\n# of Q false pos.: %s'%len(Q_pos),'\n# of Q false neg.: %s'%len(Q_neg))
                    print('indices corresponding to bin edges\nSF [pos,neg]: %s'%num_bins_SF_index,'\nQ: %s'%num_bins_Q_index)
            #
            ## convert to arrays
            num_bins_SF = np.empty_like(num_bins_SF_index,dtype='float32')
            num_bins_Q = np.empty_like(num_bins_Q_index,dtype='float32')
            #
            for ii in range(len(num_bins_SF[0])):
                num_bins_SF[0][ii] = SF_pos[num_bins_SF_index[0][ii]]
                num_bins_SF[1][ii] = SF_neg[num_bins_SF_index[1][ii]]
                num_bins_Q[0][ii] = Q_pos[num_bins_Q_index[0][ii]]
                num_bins_Q[1][ii] = Q_neg[num_bins_Q_index[1][ii]]
            #
            ## set first/last entry as limits of mass range for smf
            num_bins_SF[0][0] = range2[0]         # reminder: row [0] = false pos; row [1] = false neg
            num_bins_SF[1][0] = range2[0]
            num_bins_SF[0][-1] = range2[-1] 
            num_bins_SF[1][-1] = range2[-1]
            num_bins_Q[0][0] = range2[0]         # reminder: row [0] = false pos; row [1] = false neg
            num_bins_Q[1][0] = range2[0]
            num_bins_Q[0][-1] = range2[-1]
            num_bins_Q[1][-1] = range2[-1]
            #
            if diag_flag == 1 or project_diagnostic_flag == 1:
                if project_diagnostic_flag == 0:
                    pass
                else:
                    print('Bin edges\nSF [pos,neg]: %s'%num_bins_SF,'\nQ: %s'%num_bins_Q)
            #
            ## compute midpoints between false pos/neg bin edges
            bin_edge_means_SF = np.mean(num_bins_SF,axis=0)
            bin_edge_means_Q = np.mean(num_bins_Q,axis=0)
            #
            #
            ## Call the file that now stores Section (5) - compute midbins, correction factors, & metric of merit
            #
            #
            exec(open('spec_asymmetric_binning_metric.py').read())
            #
            ######
            #####if np.sum([np.isnan(SF_var),np.isnan(Q_var)]) == 0:
            #####    tot_var = SF_var + Q_var
            #####else:
            #####    tot_var = float('NaN')
            ######
            #
            ## SECIONT (4.1): prepare what to WRITE to file
            bin_entry1 = str(np.round(z_cutoff[0],decimals=3))+delim+str(np.round(z_cutoff[1],decimals=3))+delim+'SF'+delim+str(limiting_mass)+delim+'%.5f'%SF_var+delim+str(np.sum([mem_spec[0],mem_phot[0]]))+delim+'%s'%len(SF_pos)+delim+'%s'%len(SF_neg)+delim+str(method)+delim+str(num_bins_to_try[number])+delim+str(bins_SF)
            bin_entry2 = str(np.round(z_cutoff[0],decimals=3))+delim+str(np.round(z_cutoff[1],decimals=3))+delim+'Q'+delim+str(limiting_mass)+delim+'%.5f'%Q_var+delim+str(np.sum([mem_spec[1],mem_phot[1]]))+delim+'%s'%len(Q_pos)+delim+'%s'%len(Q_neg)+delim+str(method)+delim+str(num_bins_to_try[number])+delim+str(bins_Q)
            writer = '%s'%str(bin_entry1)+'\n'+'%s'%str(bin_entry2)+'\n'
            f.write(writer)
            #
        #
## use built-in "Freedman Diaconis Estimator" for np.histrogram_bin_edges ASYMMETRIC binning method 
#method = 4
##
## make histograms
#SF_pos_hist, bins_SF = np.histogram_bin_edges(SF_pos, bins='auto', range=range2)
#SF_neg_hist, bins_SF = np.histogram_bin_edges(SF_neg, bins='auto', range=range2)
#Q_pos_hist, bins_Q = np.histogram_bin_edges(Q_pos, bins='auto', range=range2)
#Q_neg_hist, bins_Q = np.histogram_bin_edges(Q_neg, bins='auto', range=range2)
##
# bins_SF = np.round(bins_SF,decimals=3)
# bins_Q = np.round(bins_Q,decimals=3)
##
#SF_ratio = np.round((SF_pos_hist/SF_neg_hist),decimals=3)
#Q_ratio = np.round((Q_pos_hist/Q_neg_hist),decimals=3)
#for jj in range(len(SF_pos_hist)):
#    if SF_pos_hist[jj] == 0 and SF_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
#        SF_ratio[jj] = float('NaN')
#for jj in range(len(Q_pos_hist)):
#    if Q_pos_hist[jj] == 0 and Q_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
#        Q_ratio[jj] = float('NaN')
##
### compute variance of SF/Q ratios from 1
#SF_var = np.round(np.sum((1 - SF_ratio)**2),decimals=3)
#Q_var = np.round(np.sum((1 - Q_ratio)**2),decimals=3)
#if np.isnan(SF_var) == 1 or np.isnan(Q_var) == 1:
#   total_var = float('NaN')
#else: total_var = np.round((SF_var+Q_var),decimals=3)
### prepare what to write to file
#bin_entry1 = str(np.round(z_cutoff[0],decimals=3))+delim+str(np.round(z_cutoff[1],decimals=3))+delim+'SF'+delim+str(method)+delim+str(len(bins_SF))+delim+str(bins_SF)+delim+str(np.sum(mem[0]))+delim+str(np.round(SF_var,decimals=3))+delim+str(np.round(total_var,decimals=3))
#bin_entry2 = str(np.round(z_cutoff[0],decimals=3))+delim+str(np.round(z_cutoff[1],decimals=3))+delim+'Q'+delim+str(method)+delim+str(len(bins_Q))+delim+str(bins_Q)+delim+str(np.sum(mem[1]))+delim+str(np.round(Q_var,decimals=3))+delim+str(np.round(total_var,decimals=3))
#writer = '%s'%str(bin_entry1)+'\n'+'%s'%str(bin_entry2)+'\n'
#f.write(writer)
#
#
#
## close the file       
f.close()
#
#
#
## TIME_FLAG END
#
if time_flag == 1:
    if project_time_flag == 1:
        pass
    else:
        print('Program "spec_completeness_binning.py" for z_spec < %.3f'%z_cutoff[0],' and z_phot < %.3f'%z_cutoff[1],' took: %s seconds to run.\n\n' % (time.time() - start_time))
#
#
#
#
## enable diagnostic output
if diag_flag == 1 or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:    
        print('\n\n"spec_completeness_binning.py" terminated successfully.\n')
#
#                      
###### PROGRAM END ######