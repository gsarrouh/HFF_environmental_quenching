tester = 'tester_tester_tester'
#Created on Fri Jun 27 13:17:46 2020
#
################## spec_mass_binning.py ##################
#
#
#
### This program used to be SECTION (3.2) from "master_smfz*.py". Now it stands alone, and is called by "master_data*.py". For a given spectroscopic & photometric redshift cut, a population of false positives & false negatives results. The following program uses 3 different binning techniques to try various numbers of bins, evaluating the resulting spectroscopic correction factors, and computing a "mertic of merit", by which we may judge various combinations of redshift cuts, bins, and binning techniques.
#
## This is the old description of SECTION (3.2):
### SECTION (3.2): calculate SPECTROSCOPIC COMPLETENESS correction. basically, look at all the false positives/false negatives, and sort them by type (i.e. SF/Q). then bin them (i.e. make histograms of false pos/neg for each of SF/Q). take their ratio of false pos to false neg, and plot that ratio. it is the correction factor to be applied to the photometric subsample
#
## SPEC. BINNING: iterate through different number of histogram bins to see which yields a set of corrections closest in general to ~1; this has now been expanded to test: 
#
## method = 1, SYMMETRIC binning; 
## method = 2, ASYMMETRIC binning, EQUAL NUMBER of objects in each bin; 
## method = 3, ASYMMETRIC binning, EQUAL MASS in each bin; 
## method = 4; ASYMMETRIC binning using built-in "np.histogram_bin_edges" function, using black-box binning method "auto", which is the Maximum of the "STURGES" & "FREEDMAN DIACONIS" estimators
#
#
### Section summary:
#
### PROGRAM START
#
### EVALUATE DIFFERENT BIN #s & BINNING TECHNIQUES:
### (0)    modules & definitions
### (1)    sort data into false pos/neg lists for SF/Q;
### (2)    Method 1: SYMMETRIC bins; 
### (3)    Method 2: ASYMMERTIC bins, EQUAL NUMBER distribution among bins;
### (4)    Method 3: ASYMMERTIC bins, EQUAL MASS distribution among bins;
### (5)    compute BIN EDGES and METRIC OF MERIT;
### (6)    prepare strings & to WRITE to file;
#
### PROGRAM END
#
## NOTE: there are flags for diagnostics and plotting throughout the script. search "MAY NEED TO EDIT" to identify where these flags are
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
    start_time = time.time()
#
#
#
## SECTION (1): setup LISTS, i.e. compile lists of SF/Q false pos/neg populations; define "midbins" function
# 
## ignore 'divide by zero' errors, since they are handled with exceptions below
np.seterr(divide='ignore')
#
limiting_mass = [7.5,7.8,8.0,7.5,7.4,7.3] # clusters 1,2,3,4,5,6, see IDs below
range2 = [7.3,12.3]
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
header1 = 'z_spec_cutoff'+delim+'z_phot_cutoff'+delim+'type'+delim+'binning_method'+delim+'#_of_bins'+delim+'bin_edges'+delim+'cluster_members'+delim+'var'+delim+'tot_var\n'
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
num_bins_to_try = np.arange(2,5,1)   # np.arange([...]) = [array_start,array_end_plus_increment,increment]
# write a loop that interatively uses a different number of bins in the histrogram, to isolate the largest number for which all bins have at least one entry; NOTE: the lines that stop the loop have been commented out, to investigate the relative fraction of false pos/neg for each different # of bins
#
# 
method = 1
for number in range(len(num_bins_to_try)):
    #
    # make histograms
    SF_pos_hist, bins_SF = np.histogram(SF_pos, bins=num_bins_to_try[number], range=range2)
    SF_neg_hist, bins_SF = np.histogram(SF_neg, bins=num_bins_to_try[number], range=range2)
    Q_pos_hist, bins_Q = np.histogram(Q_pos, bins=num_bins_to_try[number], range=range2)
    Q_neg_hist, bins_Q = np.histogram(Q_neg, bins=num_bins_to_try[number], range=range2)
    #
    bins_SF = np.round(bins_SF,decimals=3)
    bins_Q = np.round(bins_Q,decimals=3)
    #
    SF_ratio = np.round((SF_pos_hist/SF_neg_hist),decimals=3)
    Q_ratio = np.round((Q_pos_hist/Q_neg_hist),decimals=3)
    for jj in range(len(SF_pos_hist)):
        if SF_pos_hist[jj] == 0 and SF_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
            SF_ratio[jj] = float('NaN')
    for jj in range(len(Q_pos_hist)):
        if Q_pos_hist[jj] == 0 and Q_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
            Q_ratio[jj] = float('NaN')
    # compute midbins for spec. mass completeness plot (i.e. plot of false pos/false neg ratios)
    SF_frac_midbins = midbins(bins_SF)
    Q_frac_midbins = midbins(bins_Q)
    #
    #
    #    
    ## call a new file: "correction_factors.py" to interpolate/extrapolate the correction factors to the SMF, and return the metric of merit
    #
    exec(open('correction factor.py').read())      #opens and executes the script
    #
    #
    ## compute variance of SF/Q ratios from 1
    SF_var = SF_metric
    Q_var = Q_metric
    #
    if np.isnan(SF_var) == 1 or np.isnan(Q_var) == 1:
        total_var = float('NaN')
    else: total_var = np.round((SF_var+Q_var),decimals=3)
    #
    #
    ## prepare what to write to file
    bin_entry1 = str(np.round(z_cutoff[0],decimals=3))+delim+str(np.round(z_cutoff[1],decimals=3))+delim+'SF'+delim+str(method)+delim+str(num_bins_to_try[number])+delim+str(bins_SF)+delim+str(np.sum(mem[0]))+delim+str(np.round(SF_var,decimals=3))+delim+str(np.round(total_var,decimals=3))
    bin_entry2 = str(np.round(z_cutoff[0],decimals=3))+delim+str(np.round(z_cutoff[1],decimals=3))+delim+'Q'+delim+str(method)+delim+str(num_bins_to_try[number])+delim+str(bins_Q)+delim+str(np.sum(mem[1]))+delim+str(np.round(Q_var,decimals=3))+delim+str(np.round(total_var,decimals=3))
    writer = '%s'%str(bin_entry1)+'\n'+'%s'%str(bin_entry2)+'\n'
    f.write(writer)
#
#
## ASYMMETRIC BINNING START
#
method_designations = [2,3]
#
##
for m in range(len(method_designations)):
    method = method_designations[m]
    if method == 2:
        #
        #
        #
        #
        #
### SECTION (3): METHOD 2
### METHOD 2: ASYMMETRIC BINS for bins with ~EQUAL NUMBER OF OBJECTS in each bin (for each of the SF/Q false pos/neg pairs), return the asymmetric bin edges and compute the midpoints b/w the *_pos/*_neg bin edges. Use these midpoints as the bin edges for a new histogram for the *_pos/*_neg lists, and re-compute the false pos/neg ratio for each of SF/Q. Then print relevant data to a file.
#
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
            ## convert to arrays
            num_bins_SF_index = np.array(num_bins_SF_index)
            num_bins_Q_index = np.array(num_bins_Q_index)
            #print('# bins: %s'%(num_bins_to_try[number]-1))
            #print('length of num_bins_* array: %s'%len(num_bins_SF_index[0]),'\nlength of SF array: %s'%len(SF_pos),'\nlength of Q array: %s'%len(Q_pos))
            #print('indices corresponding to bin edges\nSF: %s'%num_bins_SF_index,'\nQ: %s'%num_bins_Q_index)
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
            ## compute midpoints between false pos/neg bin edges
            bin_edge_means_SF = np.mean(num_bins_SF,axis=0)
            bin_edge_means_Q = np.mean(num_bins_Q,axis=0)
    else:
        pass
    #
    #
    #
    #
    ### SECTION (5): compute BIN EDGES, CORRECTION FACTORS, & METRIC OF MERIT
    ## The goal now is to produce a list of bin edges that will be applied to both the false pos/neg lists. for the ii'th bin edge from each of the false pos/neg lists, test all possible values for the ii'th bin between false_pos_bin_edge[ii] & false_neg_bin_edge[ii] in increments of 0.01. choose the bin edge which yields the lowest average deviation from 1. then call (i.e. execute the file) "correction_factors.py", which is based on code from "master_smfz*" S3.2. It interpolates/extrapolates SMF bin correction factors from the false pos/neg ratios computed below, and returns the "metric of merit", i.e. (the sum of squared deviations from 1 for each pos/neg ratio bin) divided by (the number of bins). 
    #
    ## SF
    for ii in range(len(bin_edge_means_SF)):
        SF_var_list = []
        if num_bins_SF[0][ii] < num_bins_SF[1][ii]:    # this if statement won't apply to the first/last bin, since they have been set to equal "range2", the mass range for the study
            bin_edge_SF_start = num_bins_SF[0][ii]     # store the lower bound of the ii'th bin edge b/w the false pos/neg bin edges
            bin_edge_SF_end = num_bins_SF[1][ii]
            equal_flag = 0
        elif num_bins_SF[0][ii] > num_bins_SF[1][ii]:
            bin_edge_SF_start = num_bins_SF[1][ii]
            bin_edge_SF_end = num_bins_SF[0][ii]
            equal_flag = 0
        else:                                          # for when the i'th bin edge is the same for both false pos/neg lists
            equal_flag = 1
            bin_edge_SF_start = num_bins_SF[0][ii]
            bin_edge_SF_end = num_bins_SF[0][ii]
            mid_pt_to_use = num_bins_SF[0][ii]
        print('Start of : %.2f'%bin_edge_SF_start)
        #print(bin_edge_SF_start)
        #print(bin_edge_SF_end)
        #
        ## 
        while bin_edge_SF_start < bin_edge_SF_end:
            if ii == 0:
                bin_edge_means_SF[ii] = range2[0]    
            elif ii == (len(bin_edge_means_SF)-1):
                bin_edge_means_SF[ii] = range2[1]  
            else:
                bin_edge_means_SF[ii] = bin_edge_SF_start
            try:     # try making histograms. if it doesn't work because the bin_edge* array isn't monotonically increasing, raise an error and record a variance of 'NaN'
                SF_pos_hist, bins_SF = np.histogram(SF_pos, bins=bin_edge_means_SF, range=range2)
                SF_neg_hist, bins_SF = np.histogram(SF_neg, bins=bin_edge_means_SF, range=range2)
                bins_SF = np.round(bins_SF,decimals=3)
                ## find ratios, set ratio==1 for bins with no false pos or false neg
                SF_ratio = np.round((SF_pos_hist/SF_neg_hist),decimals=3)
                for jj in range(len(SF_pos_hist)):      # go through each bin 
                    if np.isnan(bins_SF[jj]) == 1:
                        SF_ratio[jj] = float('NaN')
                    elif SF_pos_hist[jj] == 0 and SF_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
                        SF_ratio[jj] = 1
                    elif SF_pos_hist[jj] == 0 or SF_neg_hist[jj] == 0:
                        SF_ratio[jj] = float('NaN') 
                #
                ## compute variance of SF/Q ratios from 1
                SF_var = np.sum((1 - SF_ratio)**2)
                #
                ## store variance in array for comparison of different values of bin edges
                SF_var_list.append([SF_var,bin_edge_means_SF[ii],num_bins_to_try[number]])
                #print('No error. appended to list for # of bins: %s'%num_bins_to_try[number])
                #break
            except:
                print("ERROR 1: SF false pos/neg bins overlap s.t. bin edges don't increase\nmonotonically for cutoffs - spec: %s"%z_cutoff[0],";  phot: %s"%z_cutoff[1],";   for %s"%num_bins_to_try[number]," bins, method = %s"%method,'\nUSE FEWER BINS\n"NaN" appended to list.')
                bins_SF = np.array([float('NaN')]*len(num_bins_SF[0]))
                SF_var = float('NaN')
                SF_var_list.append([SF_var,bin_edge_means_SF[ii],num_bins_to_try[number]])
                #print('error_flag=1 appended list for # of bins: %s'%num_bins_to_try[number])
                #break            
            bin_edge_SF_start = np.round(bin_edge_SF_start+0.01,decimals=2)
            print(bin_edge_SF_start)
        print(bin_edge_means_SF)
        if equal_flag == 0:
            #print('SF_var_list: %s'%SF_var_list)
            #    
            SF_var_list.sort(key=lambda x: x[0])
            mid_pt_to_use = SF_var_list[0][1]        # assign best midpoint to bin_edge* array 
        if ii == 0 or ii == (len(bin_edge_means_SF)-1):
            pass
        else:
            bin_edge_means_SF[ii] = mid_pt_to_use 
        #print('bin_edge_means_SF[%s]'%ii,' = %s'%bin_edge_means_SF[ii])
        if ii != (len(bin_edge_means_SF)-1):
            if bin_edge_means_SF[ii] > bin_edge_means_SF[(ii+1)]:
                print("ERROR 2: the next bin value is lower than the current bin value for cutoffs - spec: %s"%z_cutoff[0],";  phot: %s"%z_cutoff[1],";   for %s"%num_bins_to_try[number]," bins, method = %s"%method,'\nUSE FEWER BINS')
                bin_edge_means_SF[ii] = bin_edge_means_SF[(ii+1)]            
    #
    ##
    ### insert code here calling new correction factors file    
        
        
        
#        **************************
        
    
    
        
    
    
    
    
    #       
    ## Q
    for ii in range(len(bin_edge_means_Q)):
        Q_var_list = []
        if num_bins_Q[0][ii] < num_bins_Q[1][ii]:    # this if statement won't apply to the first/last bin, since they have been set to equal "range2", the mass range for the study
            bin_edge_Q_start = num_bins_Q[0][ii]     # store the lower bound of the ii'th bin edge b/w the false pos/neg bin edges
            bin_edge_Q_end = num_bins_Q[1][ii]
            equal_flag = 0
        elif num_bins_Q[0][ii] > num_bins_Q[1][ii]:
            bin_edge_Q_start = num_bins_Q[1][ii]
            bin_edge_Q_end = num_bins_Q[0][ii]
            equal_flag = 0
        else:                                          # for when the i'th bin edge is the same for both false pos/neg lists
            equal_flag = 1
            bin_edge_Q_start = num_bins_Q[0][ii]
            bin_edge_Q_end = num_bins_Q[0][ii]
            mid_pt_to_use = num_bins_Q[0][1]
        print('Start of : %.2f'%bin_edge_Q_start)
        print(bin_edge_Q_start)
        print(bin_edge_Q_end)
        #
        ## 
        while bin_edge_Q_start < bin_edge_Q_end:
            if ii == 0:
                bin_edge_means_Q[ii] = range2[0]    
            elif ii == (len(bin_edge_means_Q)-1):
                bin_edge_means_Q[ii] = range2[1]  
            else:
                bin_edge_means_Q[ii] = bin_edge_Q_start
            try:     # try making histograms. if it doesn't work because the bin_edge* array isn't monotonically increasing, raise an error and record a variance of 'NaN'
                Q_pos_hist, bins_Q = np.histogram(Q_pos, bins=bin_edge_means_Q, range=range2)
                Q_neg_hist, bins_Q = np.histogram(Q_neg, bins=bin_edge_means_Q, range=range2)
                bins_Q = np.round(bins_Q,decimals=3)
                ## find ratios, set ratio==1 for bins with no false pos or false neg
                Q_ratio = np.round((Q_pos_hist/Q_neg_hist),decimals=3)
                for jj in range(len(Q_pos_hist)):      # go through each bin 
                    if np.isnan(bins_Q[jj]) == 1:
                        Q_ratio[jj] = float('NaN')
                    elif Q_pos_hist[jj] == 0 and Q_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
                        Q_ratio[jj] = 1
                    elif Q_pos_hist[jj] == 0 or Q_neg_hist[jj] == 0:
                        Q_ratio[jj] = float('NaN') 
                #
                ## compute variance of SF/Q ratios from 1
                Q_var = np.sum((1 - Q_ratio)**2)
                #
                ## store variance in array for comparison of different values of bin edges
                Q_var_list.append([Q_var,bin_edge_means_Q[ii],num_bins_to_try[number]])
                #print('No error. appended to list for # of bins: %s'%num_bins_to_try[number])
                #break
            except:
                print("ERROR 1: Q false pos/neg bins overlap s.t. bin edges don't increase\nmonotonically for cutoffs - spec: %s"%z_cutoff[0],";  phot: %s"%z_cutoff[1],";   for %s"%num_bins_to_try[number]," bins, method = %s"%method,'\nUSE FEWER BINS\n"NaN" appended to list.')
                bins_Q = np.array([float('NaN')]*len(num_bins_Q[0]))
                Q_var = float('NaN')
                Q_var_list.append([Q_var,bin_edge_means_Q[ii],num_bins_to_try[number]])
                #print('error_flag=1 appended list for # of bins: %s'%num_bins_to_try[number])
                #break            
            bin_edge_Q_start = np.round(bin_edge_Q_start+0.01,decimals=2)
            print(bin_edge_Q_start)
        print(bin_edge_means_Q)
        if equal_flag == 0:
            #print('SF_var_list: %s'%Q_var_list)
            #    
            Q_var_list.sort(key=lambda x: x[0])
            mid_pt_to_use = Q_var_list[0][1]        # assign best midpoint to bin_edge* array 
        if ii == 0 or ii == (len(bin_edge_means_Q)-1):
            pass
        else:
            bin_edge_means_Q[ii] = mid_pt_to_use 
        #print('bin_edge_means_Q[%s]'%ii,' = %s'%bin_edge_means_Q[ii])
        if ii != (len(bin_edge_means_Q)-1):
            if bin_edge_means_Q[ii] > bin_edge_means_Q[(ii+1)]:
                print("ERROR 2: the next bin value is lower than the current bin value for cutoffs - spec: %s"%z_cutoff[0],";  phot: %s"%z_cutoff[1],";   for %s"%num_bins_to_try[number]," bins, method = %s"%method,'\nUSE FEWER BINS')
                bin_edge_means_Q[ii] = bin_edge_means_Q[(ii+1)]
    #
    ##
    ### insert code here calling new correction factors file    
        
        
        
#        **************************
        
    
    
        
    #
    #
    #
    if np.isnan(SF_var) == 1 or np.isnan(Q_var) == 1:
        total_var = float('NaN')
    else: 
        total_var = np.round((SF_var + Q_var),decimals=3)
    #
    #
    #
    #
    #
    ## SECIONT (6): prepare what to WRITE to file
    bin_entry1 = str(np.round(z_cutoff[0],decimals=3))+delim+str(np.round(z_cutoff[1],decimals=3))+delim+'SF'+delim+str(method)+delim+str(num_bins_to_try[number])+delim+str(bins_SF)+delim+str(np.sum(mem[0]))+delim+str(np.round(SF_var,decimals=3))+delim+str(np.round(total_var,decimals=3))
    bin_entry2 = str(np.round(z_cutoff[0],decimals=3))+delim+str(np.round(z_cutoff[1],decimals=3))+delim+'Q'+delim+str(method)+delim+str(num_bins_to_try[number])+delim+str(bins_Q)+delim+str(np.sum(mem[1]))+delim+str(np.round(Q_var,decimals=3))+delim+str(np.round(total_var,decimals=3))
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
#bins_SF = np.round(bins_SF,decimals=3)
#bins_Q = np.round(bins_Q,decimals=3)
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
    print('Program "spec_completeness_binning.py" for z_spec < %.3f'%z_cutoff[0],' and z_phot < %.3f'%z_cutoff[1],' took: %s seconds to run.\n\n' % (time.time() - start_time))
#
#
#
#
print('\n\nProgram terminated successfully.')
#
#                      
###### PROGRAM END ######