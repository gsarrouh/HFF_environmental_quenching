#Created on Tue Jul 07 20:29:15 2020
#
#
################## master_smfz_9.py ##################
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
## v9 FINAL version, done after the kinks in variational analysis have been ironed out. This file
##    now just executes a given redshift cut and bin edges, as defined in "main_project_file.py"
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
### (3.4)  compute QUENCHED FRACTION (i.e. # Q galaxies / (total # SF+Q galaxies) in each mass bin
### (4)    add ERROR BARS for scatter plot;
### (5)    EMCEE simulation; see emcee_chi2_final.py;
### (6)    build SCHECTER best-fit models;
### (7)    PLOT that shit - SMF (CLUSTER v FIELD);
### (7.1)   PLOT that shit - SMF by POPULATION;
#
### PROGRAM END
#
## NOTE: there are flags for diagnostics and plotting throughout the script. search "MAY NEED TO EDIT" to identify where these flags are
#
## NOTE: search "MAY NEED TO EDIT" to find where user-input is required
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
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
from astropy.table import Table
from astropy.table import Column
#from scipy.optimize import curve_fit
#
#
#
# define a function to compute the incremental difference of applying a correction to the raw count lists
def correction_difference(raw_smf,completeness_correction):
    corrected_smf = raw_smf*completeness_correction
    diff = corrected_smf- raw_smf
    return diff
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
## define a function to mind the max value in a list of lists of various lengths
def nested_list_max(list0):
    max_0 = 0
    #
    for ii in range(len(list0)):       # for each list in list of lists
        for jj in range(len(list0[ii])):    # go through each list individually 
            if list0[ii][jj] > max_0:
                max_0 = list0[ii][jj]
    return max_0
#
## define a function to compute multiplicative limiting mass corrections by mass bin
#
def mass_completeness_correction_function(mass_bin_edges,limiting_masses):
    mass_completeness_correction_factors = np.zeros([(len(mass_bin_edges)-1),1],dtype='float64')
    for ii in range(len(mass_completeness_correction_factors)):
        for jj in range(len(limiting_masses)):
            if limiting_masses[jj] <= mass_bin_edges[ii]:    # count # of clusters complete at each mass bin
                mass_completeness_correction_factors[ii]+=1
    mass_completeness_correction_factors = np.transpose(6/mass_completeness_correction_factors) # take recipricol for multiplicative factors
    return mass_completeness_correction_factors
#
## define a function to take in a list of 6 SMFs, normalize them all by MASS to a total mass of 1 per SMF, and then combines them into a single SMF and divides by 6, effectively taking the average SMF, whose area-under-the-curve (i.e. total mass) has been normalized to 1; see Section (3.3)
#
def normalize_smf_mass(SF_raw_smf,Q_raw_smf,midbins):
    #
    ## compute the total mass (i.e. sum((# of galaxies in each bin)*(bin mass))
    total_mass = np.array([0]*6,dtype='float64')
    for ii in range(N):          # go through clusters 1 at a time
        total_mass[ii] = np.sum((SF_raw_smf[ii]+Q_raw_smf[ii])*np.transpose(midbins))
    #
    ## now compute the raw relative mass fraction in each bin, for SF/Q by cluster; that is, create an array (row_1=SF, row_2=Q)
    SF_rel_mass = np.empty_like(SF_raw_smf,dtype='float64')
    Q_rel_mass = np.empty_like(Q_raw_smf,dtype='float64')
    for ii in range(N):          # go through clusters 1 at a time
        for jj in range(M):      # go through each mass bin one at a time
            SF_rel_mass[ii][jj] = ((SF_raw_smf[ii][jj]*midbins[jj]) / total_mass[ii])   # mass in the jj'th bin of cluster ii, divided by total mass of cluster ii
            Q_rel_mass[ii][jj] = ((Q_raw_smf[ii][jj]*midbins[jj]) / total_mass[ii]) 
    #
    ## Mass in a bin is equal to: (# count in that bin) * (mass value of that bin). The above normalizes by mass, so the arrays "*_normalied_mass" contain the normalized amount of MASS IN EACH MASS BIN. To get the normalized # count, you need to divide (the normalized amount of mass in each bin) by (the mass value of that bin)
    ## Normalize SMFs by cluster. 
    SF_smf_by_cluster = np.empty_like(SF_raw_smf,dtype='float64')     #initialize arrays for FINAL SMF
    Q_smf_by_cluster = np.empty_like(Q_raw_smf,dtype='float64')
    ## normalize each cluster individually, for both SF & Q
    #
    for ii in range(N):          # go through clusters 1 at a time
        for jj in range(M):      # go through each mass bin (within each cluster) 1 at a time
            SF_smf_by_cluster[ii][jj] = (SF_rel_mass[ii][jj] / midbins[jj])  # normalize each cluster by the TOTAL MASS IN EACH CLUSTER (hence the midbins multiplication); 
            Q_smf_by_cluster[ii][jj] = (Q_rel_mass[ii][jj] / midbins[jj]) 
    #
    ## for ease of future calling, combine all clusters into a sinlge SMF
    SF_smf = (np.sum(SF_smf_by_cluster,axis=0)/6)
    Q_smf = (np.sum(Q_smf_by_cluster,axis=0)/6)
    total_smf = SF_smf + Q_smf
    return SF_smf, Q_smf, total_smf
    #
#####
#
## define a function to take in a list of 6 SMFs, normalize them all by NUMBER COUNT of 1 galaxy per SMF, and then combines them into a single SMF and divides by 6, effectively taking the average SMF, whose area-under-the-curve (i.e. total # of galaxies) has been normalized to 1; see Section (3.3)
#
def normalize_smf_count(SF_raw_smf,Q_raw_smf,midbins):
## Add mass & spec corrections to *_raw_smf lists 
    #
    ## compute the total number count (i.e. sum((# of galaxies in each bin))
    total_count = np.array([0.0]*N)
    for ii in range(N):
        for jj in range(M):
            total_count[ii] = total_count[ii] + (SF_raw_smf[ii][jj]+Q_raw_smf[ii][jj])
    #
    ## now compute the raw relative mass fraction in each bin, for SF/Q by cluster; that is, create an array (row_1=SF, row_2=Q)
    SF_rel_smf = np.empty_like(SF_raw_smf)
    Q_rel_smf = np.empty_like(Q_raw_smf)
    for ii in range(N):          # go through clusters 1 at a time
        for jj in range(M):      # go through each mass bin one at a time
            SF_rel_smf[ii][jj] = (SF_raw_smf[ii][jj]/ total_count[ii])   # count in the jj'th bin of cluster ii, divided by total mass of cluster ii
            Q_rel_smf[ii][jj] = (Q_raw_smf[ii][jj]/ total_count[ii]) 
    #
    ## compute relative mass fractions of SF/Q by cluster after normalization
    count_fraction_by_cluster = np.array([[0]*6]*2,dtype='float32')
    #
    count_fraction_by_cluster[0] = np.sum(SF_rel_smf, axis=1)
    count_fraction_by_cluster[1] = np.sum(Q_rel_smf, axis=1)
    #
    ## Now combine the 6 SMFs into one, and divide by 6
    SF_smf = (np.sum(SF_rel_smf,axis=0)/6)
    Q_smf = (np.sum(Q_rel_smf,axis=0)/6)
    total_smf = SF_smf + Q_smf
    return SF_smf, Q_smf, total_smf, count_fraction_by_cluster
    #
#####
#
#
## define a function to interpolate the relative errors of the Spec. Completeness Correction ratios (n bins) to the SMF mass bin correction factors (N bins), where n < N. see Section (4); 
#
def interpolate_errors(ratios,ratios_err,correction_factors,ratio_bins,mass_bins):
    correction_factors_err = np.empty_like(correction_factors)
    for ii in range(len(correction_factors)): 
        for jj in range(len(ratio_bins)):
            if mass_bins[ii] < ratio_bins[0]:
                correction_factors_err[ii] = (ratios_err[0] / ratios[0]) * correction_factors[ii]
            elif mass_bins[ii] > ratio_bins[-1]:
                correction_factors_err[ii] = (ratios_err[-1] / ratios[-1]) * correction_factors[ii]
            elif mass_bins[ii] >= ratio_bins[jj] and mass_bins[ii] < ratio_bins[jj+1]:
                D = ratio_bins[jj+1] - ratio_bins[jj]
                d = mass_bins[ii] - ratio_bins[jj]
                correction_factors_err[ii] = np.sqrt( ((D-d)/D)*(ratios_err[jj]/ratios[jj])**2 + (d/D)*(ratios_err[jj+1]/ratios[jj+1])**2 )*correction_factors[ii]
    return correction_factors_err
#
#
#
#
#
### TEMPORARY WRITING FLAG
where_im_at_flag = 0
#
#
## MAY NEED TO EDIT: NORMALIZATION FLAG
## Choose your normalization:   1=by mass;  2=by volume;   3=by number count (i.e. # of galaxies in each cluster)
normalization_flag = 3
#
#
## MAY NEED TO EDIT: choose the method to CALCULATE ERRORS
smf_error_method = 3             #   1 = propagate errors;   2 = bootstrap resampling
spec_completeness_error_flag = 2     # 1 = sum of rel. errors; 2 = (rel. error of numerator)^2 + (rel. error of denom.)^2
#
#
## MASTER DIAGNOSTIC FLAG: allows you to "turn off" all diagnostics at once (if equals 0), "turn on" all flags (if equal to 1), or set diagnostic flags individually (equal to 2 - search "MAY NEED TO EDIT" to find flags)
diag_flag_master = 2       # 0= all flags turned off;     1= all flags turned on;     2= may turn on flags individually
#
## diagnostic flags
diag_flag_1 = 1            # counting array through initial loop, tracking all object categories
diag_flag_2 = 1            # add spec & phot subsampes together for each cluster, and ensure they equal the total raw count in each mass bin
diag_flag_3 = 1            # limiting mass completeness correction factors
diag_flag_4 = 1            # display bootstrap sampling results
diag_flag_5 = 1            # spec completeness correction factors
diag_flag_6 = 1            # display diagnostics before/after normalization
#
summary_flag_1 = 1         # initial Summary Table: ensure lists agree w/ "master_data*.py"
summary_flag_2 = 1
#
plot_flag_1 = 1            # spec completeness correction factors
plot_flag_2 = 1           # SMF
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
SF_pos_list = [[],[],[],[],[],[]]
SF_neg_list = [[],[],[],[],[],[]]
Q_pos_list = [[],[],[],[],[],[]]
Q_neg_list = [[],[],[],[],[],[]]
SF_lost = [[],[],[],[],[],[]]
Q_lost = [[],[],[],[],[],[]]
#
#
if 'limiting_mass' in locals():
    pass
else:
    limiting_mass = [7.62,7.63,8.2,7.5,6.75,6.64] # clusters 1,2,3,4,5,6, see IDs below; 
#
SF_pos_lost = np.array([0]*6)        # to track SF/Q false pos/neg objects lost due to their being below the mass limit, by cluster
SF_neg_lost = np.array([0]*6)
Q_pos_lost = np.array([0]*6)
Q_neg_lost = np.array([0]*6)
other_lost = np.array([0]*6)    #objects below limiting mass other than false pos/neg
#
counting_array = np.array([0]*13)
#
# The following loop searches the master catalogue 'master_cat' and separates all objects by 
# cluster. Then, it looks for all objects above the limiting mass for that cluster. It then 
# creates two lists: one for SF and one for Q (e.g. SF_*/Q_*). It further splits these lists into those objects with spectrscopy, and those without (e.g. SF*_spec/phot)
#
for cluster in range(len(limiting_mass)):          # loop through clusters one at a time; "cluster" takes on values [0,1,2,3,4,5]
    for counter in range(len(master_cat)):
        if master_cat['cluster'][counter] == (cluster+1):    # cluster #
            counting_array[0]+=1                                    # all objects in all clusters
            if master_cat['lmass'][counter] > limiting_mass[cluster]:    # limiting mass of cluster: 7.5
                counting_array[1]+=1                                # all objects above limiting mass
                if master_cat['member'][counter] == 0:    # cluster member = 0
                    counting_array[2]+=1                            # all cluster members
                    if master_cat['type'][counter] == 1:   # SF type = 1
                        counting_array[3]+=1                        # SF cluster members
                        SF_list[cluster].append(master_cat['lmass'][counter])               # SF cluster members
                        if master_cat['sub'][counter]==2:     # sub=2 is objects w/ photometry only
                            counting_array[4]+=1                    # SF PHOT cluster members
                            SF_phot_list[cluster].append(master_cat['lmass'][counter])
                        elif master_cat['sub'][counter]==1:  # sub=1 is objects w/ both spec & phot
                            SF_spec_list[cluster].append(master_cat['lmass'][counter])
                            counting_array[5]+=1                    # SF SPEC cluster members
                    elif master_cat['type'][counter] == 2: # Q type = 2
                        Q_list[cluster].append(master_cat['lmass'][counter])               # Q cluster members
                        counting_array[6]+=1                    # Q cluster members
                        if master_cat['sub'][counter]==2:     # sub=2 is objects w/ photometry only
                            Q_phot_list[cluster].append(master_cat['lmass'][counter])
                            counting_array[7]+=1                    # Q PHOT cluster members
                        elif master_cat['sub'][counter]==1:  # sub=1 is objects w/ both spec & phot
                            Q_spec_list[cluster].append(master_cat[counter]['lmass'])
                            counting_array[8]+=1                    # Q SPEC cluster members
                elif master_cat['member'][counter] == 2:  
                    if master_cat['type'][counter] == 1:
                        SF_pos_list[cluster].append(master_cat[counter]['lmass'])
                        counting_array[9]+=1                    # SF pos 
                    elif master_cat['type'][counter] ==2:
                        Q_pos_list[cluster].append(master_cat[counter]['lmass'])
                        counting_array[10]+=1                    # Q pos 
                elif master_cat['member'][counter] == 3:
                    if master_cat['type'][counter] == 1:
                        SF_neg_list[cluster].append(master_cat[counter]['lmass'])
                        counting_array[11]+=1                    # SF neg 
                    elif master_cat['type'][counter] ==2:
                        Q_neg_list[cluster].append(master_cat[counter]['lmass'])
                        counting_array[12]+=1                    # Q neg 
                elif master_cat['member'][counter] == 1:                   ### FIELD!
                    if master_cat['type'][counter] == 1:
                        SF_field_list[cluster].append(master_cat[counter]['lmass'])
                        counting_array[11]+=1                    # SF neg 
                    elif master_cat['type'][counter] ==2:
                        Q_field_list[cluster].append(master_cat[counter]['lmass'])
                        counting_array[12]+=1 
            elif master_cat['lmass'][counter] < limiting_mass[cluster]:
                if master_cat['member'][counter] == 6:    # member = 0
                    if master_cat['type'][counter] == 1:   # SF type = 1
                        SF_lost[cluster].append(master_cat['lmass'][counter])               # SF <lim. mass
                    elif master_cat['type'][counter] == 2: # Q type = 2
                        Q_lost[cluster].append(master_cat['lmass'][counter])               # Q <lim. mass
                elif master_cat['member'][counter] == 2:    # member false pos = 2
                    if master_cat['type'][counter] == 1:   # SF type = 1
                        SF_pos_lost[cluster]+=1
                    elif master_cat['type'][counter] == 2: # Q type = 2
                        Q_pos_lost[cluster]+=1
                elif master_cat['member'][counter] ==3:   # member false neg = 3
                    if master_cat['type'][counter] == 1:   # SF type = 1
                        SF_neg_lost[cluster]+=1
                    elif master_cat['type'][counter] == 2: # Q type = 2
                        Q_neg_lost[cluster]+=1
                else: other_lost[cluster]+=1                                               # catchall for everything else
#
#
if (diag_flag_1 == 1 and diag_flag_master == 2) or diag_flag_master == 1:
    print('\n[all,above_lim_mass,all_cluster_members,SF_members,SF_phot,SF_spec,Q_members,Q_phot,Q_spec,SF_pos,Q_pos,SF_neg,Q_neg]\n%s'%counting_array)
#
## SECTION (1.1) - Summary Table 1: Did we pick up all the false pos/neg as reported in "master_data*.py"?
#
if summary_flag_1 == 1 or adams_flag == 1:
    ## Summarize initial data stats in table
    num_SF_pos = np.sum([len(SF_pos_list[0]),len(SF_pos_list[1]),len(SF_pos_list[2]),len(SF_pos_list[3]),len(SF_pos_list[4]),len(SF_pos_list[5])])
    num_SF_neg = np.sum([len(SF_neg_list[0]),len(SF_neg_list[1]),len(SF_neg_list[2]),len(SF_neg_list[3]),len(SF_neg_list[4]),len(SF_neg_list[5])])
    num_Q_pos = np.sum([len(Q_pos_list[0]),len(Q_pos_list[1]),len(Q_pos_list[2]),len(Q_pos_list[3]),len(Q_pos_list[4]),len(Q_pos_list[5])])
    num_Q_neg = np.sum([len(Q_neg_list[0]),len(Q_neg_list[1]),len(Q_neg_list[2]),len(Q_neg_list[3]),len(Q_neg_list[4]),len(Q_neg_list[5])])
    #
    pos_neg_names = Column(['TOTAL False Pos.','TOTAL False Neg.','SF - False Pos.','SF - False Neg.','SF <lim. mass','Q - False Pos.','Q - False Neg.','Q <lim. mass','SUM'],name='Property')
    col_names = cluster_names
    # SF table
    pos_neg0 = Column([np.sum(pos_spec),np.sum(neg_spec),num_SF_pos,num_SF_neg,np.sum([SF_pos_lost,SF_neg_lost]),num_Q_pos,num_Q_neg,np.sum([Q_pos_lost,Q_neg_lost]),np.sum([num_SF_pos,num_SF_neg,np.sum(SF_pos_lost),np.sum(SF_neg_lost),num_Q_pos,num_Q_neg,np.sum(Q_pos_lost),np.sum(Q_neg_lost)])],name='Total')  # total column
    pos_neg_stats = Table([pos_neg_names,pos_neg0])
    for ii in range(len(mem_spec[0])):
        col = Column([np.sum([pos_spec[0][ii],pos_spec[1][ii]]),np.sum([neg_spec[0][ii],neg_spec[1][ii]]),len(SF_pos_list[ii]),len(SF_neg_list[ii]),np.sum([SF_pos_lost[ii],SF_neg_lost[ii]]),len(Q_pos_list[ii]),len(Q_neg_list[ii]),np.sum([Q_pos_lost[ii],Q_neg_lost[ii]]),np.sum([len(SF_pos_list[ii]),len(SF_neg_list[ii]),SF_pos_lost[ii],SF_neg_lost[ii],len(Q_pos_list[ii]),len(Q_neg_list[ii]),Q_pos_lost[ii],Q_neg_lost[ii]])],name=col_names[ii])
        pos_neg_stats.add_column(col)  # add columns to table one cluster at a time
    #
    #
    ## Now prepare a summary table for the cluster MEMBERS
    #
    num_SF = np.sum([len(SF_list[0]),len(SF_list[1]),len(SF_list[2]),len(SF_list[3]),len(SF_list[4]),len(SF_list[5])])
    num_SF_phot = np.sum([len(SF_phot_list[0]),len(SF_phot_list[1]),len(SF_phot_list[2]),len(SF_phot_list[3]),len(SF_phot_list[4]),len(SF_phot_list[5])])
    num_SF_spec = np.sum([len(SF_spec_list[0]),len(SF_spec_list[1]),len(SF_spec_list[2]),len(SF_spec_list[3]),len(SF_spec_list[4]),len(SF_spec_list[5])])
    num_SF_lost = np.sum([len(SF_lost[0]),len(SF_lost[1]),len(SF_lost[2]),len(SF_lost[3]),len(SF_lost[4]),len(SF_lost[5])])
    num_Q = np.sum([len(Q_list[0]),len(Q_list[1]),len(Q_list[2]),len(Q_list[3]),len(Q_list[4]),len(Q_list[5])])
    num_Q_phot = np.sum([len(Q_phot_list[0]),len(Q_phot_list[1]),len(Q_phot_list[2]),len(Q_phot_list[3]),len(Q_phot_list[4]),len(Q_phot_list[5])])
    num_Q_spec = np.sum([len(Q_spec_list[0]),len(Q_spec_list[1]),len(Q_spec_list[2]),len(Q_spec_list[3]),len(Q_spec_list[4]),len(Q_spec_list[5])])
    num_Q_lost = np.sum([len(Q_lost[0]),len(Q_lost[1]),len(Q_lost[2]),len(Q_lost[3]),len(Q_lost[4]),len(Q_lost[5])])
    #
    member_smf_names = Column(['TOTAL Members (master_data*.py)','Total SF >lim. mass','SF - Phot.','SF - Spec.','SF <lim. mass','Total Q >lim. mass','Q - Phot.','Q - Spec.','Q <lim. mass','SUM'],name='Property')
    col_names = cluster_names
    # SF table
    member_smf0 = Column([np.sum([mem_phot,mem_spec]),num_SF,num_SF_phot,num_SF_spec,num_SF_lost,num_Q,num_Q_phot,num_Q_spec,num_Q_lost,np.sum([num_SF_phot,num_SF_spec,num_Q_phot,num_Q_spec,num_SF_lost,num_Q_lost])],name='Total')  # total column
    member_smf_stats = Table([member_smf_names,member_smf0])
    for ii in range(len(mem_spec[0])):
        col = Column([np.sum([mem_phot[0][ii],mem_phot[1][ii],mem_spec[0][ii],mem_spec[1][ii]]),len(SF_list[ii]),len(SF_phot_list[ii]),len(SF_spec_list[ii]),len(SF_lost[ii]),len(Q_list[ii]),len(Q_phot_list[ii]),len(Q_spec_list[ii]),len(Q_lost[ii]),np.sum([len(SF_phot_list[ii]),len(SF_spec_list[ii]),len(SF_lost[ii]),len(Q_phot_list[ii]),len(Q_spec_list[ii]),len(Q_lost[ii])])],name=col_names[ii])
        member_smf_stats.add_column(col)  # add columns to table one cluster at a time
    #
    print('\nSummary Table 1A: False Pos./Neg.\n%s'%pos_neg_stats)
    print('NOTE: TOTALs reported in first two rows are from Summary Table 4 in "master_data*.py".\n')
    print('\nSummary Table 1B: MEMBERS\n%s'%member_smf_stats)
    #
    #
#
#
#
#
#
## SECTION (1.2): Field sample
# 
## utilize lists for field samples created in "master_data*.py" and "master_parallel*.py"
#
##### RECALL THE NAMES OF KEY LISTS:
#
##### use "SF_field_par_list/Q_field_par_list" from "master_parallel*.py"
#
##### use "SF_field_list/Q_field_list" from "master_data*.py"
#
#
## Procedure:
##### correct for limiting mass completeness separately for all 12 frames (6 cluster, 6 parallel)
#
##### normalize all 12 frames by (mass/volume/#count) to 1, take average
#
##### calculate quenched fraction
#
##### add to SMF plot
#
#
## SECTION (2): sort objects into HISTOGRAMS bins for both SF & Q populations, then sum 
## for 'total' population. then normalize each cluster SMF by the total cluster mass. compute midbins. use for total pop plot, and for relative fractions
#
## cluster populations arrays: (SF/Q/total)_smf & (SF/Q/total)_field_smf
#
## MAY NEED TO EDIT  - 
#
## DEFINE "range2"
#
max_SF = nested_list_max(SF_list)
max_Q = nested_list_max(Q_list)
range2 = [range2[0],max(max_SF,max_Q)]
#
num_points = int((round((range2[1]-range2[0])/bin_width))+1)       # compute # of data points;  bin_width set in "main_project_file.py"
num_bins = np.linspace(range2[0],range2[1],num_points)
#
## smf histograms for individual clusters
## cluster
SF_raw_smf = [[],[],[],[],[],[]]       # initialize list of lists to store histograms of SMFs
Q_raw_smf = [[],[],[],[],[],[]]
## field
SF_field_raw_smf = [[],[],[],[],[],[]]       
SF_field_par_raw_smf = [[],[],[],[],[],[]]       
Q_field_raw_smf = [[],[],[],[],[],[]]
Q_field_par_raw_smf = [[],[],[],[],[],[]]
#
for ii in range(len(SF_list)):
    SF_field_raw, mass_bins = np.histogram(SF_field_list[ii], bins=num_bins,range=range2)
    SF_field_par_raw, mass_bins = np.histogram(SF_field_par_list[ii], bins=num_bins,range=range2)
    Q_field_raw, mass_bins = np.histogram(Q_field_list[ii], bins=num_bins,range=range2)
    Q_field_par_raw, mass_bins = np.histogram(Q_field_par_list[ii], bins=num_bins,range=range2)
    SF_raw, mass_bins = np.histogram(SF_list[ii], bins=num_bins,range=range2)
    Q_raw, mass_bins = np.histogram(Q_list[ii], bins=num_bins,range=range2)
    ## cluster
    SF_raw_smf[ii].append(SF_raw)
    Q_raw_smf[ii].append(Q_raw)
    ## field
    SF_field_raw_smf[ii].append(SF_field_raw)
    SF_field_par_raw_smf[ii].append(SF_field_par_raw)
    Q_field_raw_smf[ii].append(Q_field_raw)
    Q_field_par_raw_smf[ii].append(Q_field_par_raw)
#
## convert lists to arrays so we can do math operations on them
## cluster
SF_raw_smf = np.array(SF_raw_smf)
Q_raw_smf = np.array(Q_raw_smf)
## field
SF_field_raw_smf = np.array(SF_field_raw_smf)
SF_field_par_raw_smf = np.array(SF_field_par_raw_smf)
Q_field_raw_smf = np.array(Q_field_raw_smf)
Q_field_par_raw_smf = np.array(Q_field_par_raw_smf)
## totals
total_raw_smf = SF_raw_smf + Q_raw_smf
total_field_raw_smf = SF_field_raw_smf + Q_field_raw_smf
total_field_par_raw_smf = SF_field_par_raw_smf + Q_field_par_raw_smf
#
# Display some data for total, SF, Q: 
print('\nSection 2: RAW totals - Members')
print('SF: ',str(np.sum(SF_raw_smf)))
print('Q: ',str(np.sum(Q_raw_smf)))
print('Total: ',str(np.sum(total_raw_smf)),'\n')
print('\nRAW totals - Field\nSF field (clu): %s'%np.sum(SF_field_raw_smf),'\nSF field (par): %s'%np.sum(SF_field_par_raw_smf))
print('Q field (clu): %s'%np.sum(Q_field_raw_smf),'\nQ field (par): %s'%np.sum(Q_field_par_raw_smf))
print('Total (clu): %s'%np.sum(total_field_raw_smf),'\nTotal (par): %s'%np.sum(total_field_par_raw_smf))
#
#
## section (2.1): compute MIDBINS
#
## find midpoint of hist. bins. all populations have been binned identically, so the one 'midbin' will serve for all data arrays to be plotted. for visual clarity when plotting, offset the Q_midpoints by delta_x = 0.05
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
## define mass_bins midpoints
SF_midbins = midbins(mass_bins)
Q_midbins = SF_midbins + 0.05
#
## convert lists to arrays so we can do math operations on them
SF_phot_smf = np.array(SF_phot_smf)
SF_spec_smf = np.array(SF_spec_smf)
Q_phot_smf = np.array(Q_phot_smf)
Q_spec_smf = np.array(Q_spec_smf)
#
#
### MAY NEED TO EDIT: diag_flag_2
## DIAGNOSTIC: add spec & phot subsampes together for each cluster, and ensure they equal the total raw count in each mass bin
#
if (diag_flag_2 == 1 and diag_flag_master == 2) or diag_flag_master == 1:
    # compute differences, e.g.: SF_smf1 = SF1_spec_smf + SF1_phot_smf for each mass bin. they should be the same
    SF_diff = np.array([[0]*len(SF_midbins)]*6)        # initialize array to store difference between sample & sub-samples, by cluster
    Q_diff = np.array([[0]*len(SF_midbins)]*6)
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
## cluster
#mass_completeness_correction = np.zeros_like(SF_midbins)
#for ii in range(len(mass_completeness_correction)):
#    for jj in range(len(limiting_mass)):
#        if limiting_mass[jj] <= mass_bins[ii]:    # count # of clusters complete at each mass bin
#            mass_completeness_correction[ii]+=1
#mass_completeness_correction = np.transpose(6/mass_completeness_correction)  # the correction factor is: (total # of clusters) / (# of clusters complete at that mass bin); return as a row vector
#
## field
## cluster & "cluster field sample" limiting mass corrections (b/c they are from the same fields, have the same completeness limits, so get the same completeness correction factors
mass_completeness_correction = mass_completeness_correction_function(mass_bins,limiting_mass)   # function is the same as steps taken above, but now in function form
#
## Parallel "field sample"
mass_completeness_correction_par = mass_completeness_correction_function(mass_bins,limiting_mass_par)
#
#
#
#
# compute how many objects are added to each mass bin as a result of applying the mass_completeness_correction to the *_raw_smf lists. confirm that (# added to SF) + (# added to Q) = (# added to total)
SF_mass_completeness_diff = correction_difference(SF_raw_smf,mass_completeness_correction)
Q_mass_completeness_diff = correction_difference(Q_raw_smf,mass_completeness_correction)
total_mass_completeness_diff = correction_difference(total_raw_smf,mass_completeness_correction)
#
## fields
## cluster field
SF_field_mass_completeness_diff = correction_difference(SF_field_raw_smf,mass_completeness_correction)
Q_field_mass_completeness_diff = correction_difference(Q_field_raw_smf,mass_completeness_correction)
## parallel field
SF_field_par_mass_completeness_diff = correction_difference(SF_field_par_raw_smf,mass_completeness_correction_par)
Q_field_par_mass_completeness_diff = correction_difference(Q_field_par_raw_smf,mass_completeness_correction_par)
#
### MAY NEED TO EDIT: diag_flag_3
#
if (diag_flag_3 == 1 and diag_flag_master == 2) or diag_flag_master == 1:
# Display correction factors
    print('/nSection 3.1: Mass completeness correction factors by bin\nCluster: ',str(mass_completeness_correction),'\nParallel field: %s'%mass_completeness_correction_par)
    #
    # Display some data for total, SF, Q: 
    print('Section 3.1: Galaxies added due to MASS COMPLETENESS correction')
    print('SF: ',str(np.sum(SF_mass_completeness_diff,axis=0)),'\nTotal: %s'%np.sum(SF_mass_completeness_diff),', or ',str((np.sum(SF_mass_completeness_diff)/np.sum(SF_raw_smf))*100),'%\n')
    print('Q: ',str(np.sum(Q_mass_completeness_diff,axis=0)),'\nTotal: %s'%np.sum(Q_mass_completeness_diff),', or ',str((np.sum(Q_mass_completeness_diff)/np.sum(Q_raw_smf))*100),'%\n')
    print('Total: ',str(np.sum(total_mass_completeness_diff,axis=0)),'\nTotal: %s'%np.sum(total_mass_completeness_diff),', or ',str((np.sum(total_mass_completeness_diff)/np.sum(total_raw_smf))*100),'%\n')
else:
    print('Section 3.1: Galaxies added due to MASS COMPLETENESS correction\nSF: %s'%np.sum(SF_mass_completeness_diff),'\nQ: %s'%np.sum(Q_mass_completeness_diff),'\nTotal: %s'%np.sum(total_mass_completeness_diff))
#
#
#
#
## SECTION (3.2): calculate SPECTROSCOPIC COMPLETENESS correction (cluster members only). basically, look at all the false positives/false negatives, and sort them by type (i.e. SF/Q). then bin them (i.e. make histograms of false pos/neg for each of SF/Q). take their ratio of false pos to false neg, and plot that ratio. it is the correction factor to be applied to the photometric subsample
#
#
## RECALL: spectroscopic member SMFs are stored in the list of lists: SF_spec_smf & Q_spec_smf
SF_pos = []     # track ALL false pos for plotting
SF_neg = []     # track ALL false neg for plotting
Q_pos = []
Q_neg = []
pos_by_cluster = np.array([[0]*6]*2)    #for tracking false pos/neg by cluster; row_1=SF, row_2=Q
neg_by_cluster = np.array([[0]*6]*2)
objects_below_lim_mass = np.array([0]*6)    # for tracking objects below the limiting mass of each cluster
#
## set up lists for plotting false pos/neg ratios
for ii in range(len(SF_pos_list)):
    ## plotting lists
    for jj in range(len(SF_pos_list[ii])):
        SF_pos.append(SF_pos_list[ii][jj])
    for jj in range(len(SF_neg_list[ii])):
        SF_neg.append(SF_neg_list[ii][jj])
    for jj in range(len(Q_pos_list[ii])):
        Q_pos.append(Q_pos_list[ii][jj])
    for jj in range(len(Q_neg_list[ii])):
        Q_neg.append(Q_neg_list[ii][jj])
    ## counting lists
    pos_by_cluster[0][ii] = len(SF_pos_list[ii])
    pos_by_cluster[1][ii] = len(Q_pos_list[ii])
    neg_by_cluster[0][ii] = len(SF_neg_list[ii])
    neg_by_cluster[1][ii] = len(Q_neg_list[ii])
    objects_below_lim_mass[ii] = np.sum([SF_pos_lost[ii],SF_neg_lost[ii],Q_pos_lost[ii],Q_neg_lost[ii]])
#
## setup list for spec members (i.e. sum each bin across all clusters) for both SF/Q
#
SF_spec_mem_flat = [item for sublist in SF_spec_list for item in sublist]     # flatten lists
Q_spec_mem_flat = [item for sublist in Q_spec_list for item in sublist]
#             
## Set up false pos/neg histograms
#
if (diag_flag_4 == -99 and diag_flag_master == 2) or diag_flag_master == -99:   # SYMMETRIC BINNING
    #
    ## this section - the variational analysis testing different binning methods for a varying number of bins - has been broken out into its own program, called "spec_completeness_binning.py". it is called by the "master_data_*" program when the appropriate diagnostic flag (variational_anaylsis_master_flag or project_master_variational_flag) is turned on. the result of that analysis are presented below the 'else' statement (i.e. the bin numbers chosen based on the variational analysis).
    pass
#
else:
###### 
#    RECALL: bin edges set in "main_project_file.py"; range2 set in "data_mass_completeness*.py"
######
    #
    #
    SF_mem, bins_SF = np.histogram(SF_spec_mem_flat, bins=num_bins_SF_pos_neg, range=range2)
    SF_pos_hist, bins_SF = np.histogram(SF_pos, bins=num_bins_SF_pos_neg, range=range2)
    SF_neg_hist, bins_SF = np.histogram(SF_neg, bins=num_bins_SF_pos_neg, range=range2)
    Q_mem, bins_Q = np.histogram(Q_spec_mem_flat, bins=num_bins_Q_pos_neg, range=range2)
    Q_pos_hist, bins_Q = np.histogram(Q_pos, bins=num_bins_Q_pos_neg, range=range2)
    Q_neg_hist, bins_Q = np.histogram(Q_neg, bins=num_bins_Q_pos_neg, range=range2)
    #
    #####print('SF: %s'%bins_SF)
    #####print('Q: %s'%bins_Q)
#
### MAY NEED TO EDIT: diag_flag_5; continues below
# display diagnostics
#
if (diag_flag_5 == 1 and diag_flag_master == 2) or diag_flag_master == 1:
    # sum into total list, to compare with totals reported in "spec_stats1" table from "master_data_*.py"
#    total_pos_hist = SF_pos_hist + Q_pos_hist     
#    total_neg_hist = SF_neg_hist + Q_neg_hist
    # print
    print('Section 3.2: spec.completeness correction:\n\nThe data preparation file reports:')
    print('# of SF false pos: ',str(np.sum(pos_spec[0])))
    print('# of SF false neg: ',str(np.sum(neg_spec[0])))
    print('# of Q false pos: ',str(np.sum(pos_spec[1])))
    print('# of Q false neg: ',str(np.sum(neg_spec[1])))
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
SF_ratio = np.empty_like(SF_pos_hist, dtype='float32')
Q_ratio = np.empty_like(Q_pos_hist, dtype='float32')
# compute fractions for SF/Q
#
SF_ratio = ((SF_mem + SF_neg_hist) / (SF_mem + SF_pos_hist))
Q_ratio = ((Q_mem + Q_neg_hist) / (Q_mem + Q_pos_hist))
#
# compute midbins for spec. mass completeness plot (i.e. plot of false pos/false neg ratios)
SF_ratio_midbins = midbins(bins_SF)
Q_ratio_midbins = midbins(bins_Q)
#
#
## now compute the errors for the spec. completeness plot, which is simply sqrt(N) since the spectroscopic uncertainty is Poissonian in nature. do so by computing the relative error for the false pos & false neg histograms for each of SF/Q, and then sum in quadrature to determine relative error of fractions; 
#
if spec_completeness_error_flag == 1:
    ## compute individual rel errors
    SF_relerr_mem = (np.sqrt(SF_mem))/SF_mem
    SF_relerr_pos = (np.sqrt(SF_pos_hist))/SF_pos_hist
    SF_relerr_neg = (np.sqrt(SF_neg_hist))/SF_neg_hist
    Q_relerr_mem = (np.sqrt(Q_mem))/Q_mem
    Q_relerr_pos = (np.sqrt(Q_pos_hist))/Q_pos_hist
    Q_relerr_neg = (np.sqrt(Q_neg_hist))/Q_neg_hist
    ## compute error of derived quantity              
    SF_ratio_err = np.sqrt((SF_relerr_pos**2) + (SF_relerr_neg**2) + (SF_relerr_mem**2))*SF_ratio              
    Q_ratio_err = np.sqrt((Q_relerr_pos**2) + (Q_relerr_neg**2) + (Q_relerr_mem**2))*Q_ratio         
elif spec_completeness_error_flag == 2: 
    ## compute individual rel errors
    SF_err_mem = np.sqrt(SF_mem)
    SF_err_pos = np.sqrt(SF_pos_hist)
    SF_err_neg = np.sqrt(SF_neg_hist)
    SF_num_relerr = (SF_err_mem + SF_err_neg) / (SF_mem + SF_neg_hist)
    SF_den_relerr = (SF_err_mem + SF_err_pos) / (SF_mem + SF_pos_hist)
    Q_err_mem = np.sqrt(Q_mem)
    Q_err_pos = np.sqrt(Q_pos_hist)
    Q_err_neg = np.sqrt(Q_neg_hist)
    Q_num_relerr = (Q_err_mem + Q_err_neg) / (Q_mem + Q_neg_hist)
    Q_den_relerr = (Q_err_mem + Q_err_pos) / (Q_mem + Q_pos_hist)
    ## compute error of derived quantity              
    SF_ratio_err = np.sqrt((SF_num_relerr**2) + (SF_den_relerr**2))*SF_ratio              
    Q_ratio_err = np.sqrt((Q_num_relerr**2) + (Q_den_relerr**2))*Q_ratio 
#              
#    
## Now interpolate/extrapolate between these data points
#
# initialize arrays to store slopes/intercepts for extrapolation/interpolation of spec mass completeness correction factors
m_SF = np.zeros((len(SF_ratio_midbins)-1))     
b_SF = np.zeros((len(SF_ratio_midbins)-1))
m_Q = np.zeros((len(Q_ratio_midbins)-1))     
b_Q = np.zeros((len(Q_ratio_midbins)-1))
SF_spec_completeness_correction = np.zeros_like(SF_midbins,dtype='float32')
Q_spec_completeness_correction = np.zeros_like(Q_midbins,dtype='float32')
#
## SF
for ii in range(len(SF_ratio_midbins)-1):
    m_SF[ii] = (SF_ratio[ii+1] - SF_ratio[ii]) / (SF_ratio_midbins[ii+1] - SF_ratio_midbins[ii]) # calc slope
    b_SF[ii] = SF_ratio[ii] - (SF_ratio_midbins[ii]*m_SF[ii])   # calc intercept
#
for ii in range(len(SF_midbins)):
    if SF_spec_completeness_correction[ii] == 0:     # don't overwrite cell once correction factor is computed
        if SF_midbins[ii] < SF_ratio_midbins[0]:      # extrapolate below lowest mass bin
            SF_spec_completeness_correction[ii] = m_SF[0]*SF_midbins[ii] + b_SF[0]    
        elif SF_midbins[ii] > SF_ratio_midbins[-1]:    # extrapolate above highest mass bin
            SF_spec_completeness_correction[ii] = m_SF[-1]*SF_midbins[ii] + b_SF[-1]    
        elif SF_midbins[ii] > SF_ratio_midbins[0] and SF_midbins[ii] < SF_ratio_midbins[-1]:    # interpolate in between all other points
            for jj in range(len(SF_ratio_midbins)-1):
                if SF_midbins[ii] > SF_ratio_midbins[jj] and SF_midbins[ii] < SF_ratio_midbins[jj+1]:
                    SF_spec_completeness_correction[ii] = m_SF[jj]*SF_midbins[ii] + b_SF[jj]
        else:
            print('Error in SF spec completeness correction computation. ABORT')
            break   
#

## Q
for ii in range(len(Q_ratio_midbins)-1):
    m_Q[ii] = (Q_ratio[ii+1] - Q_ratio[ii]) / (Q_ratio_midbins[ii+1] - Q_ratio_midbins[ii]) # calc slope
    b_Q[ii] = Q_ratio[ii] - (Q_ratio_midbins[ii]*m_Q[ii])   # calc intercept
#
for ii in range(len(SF_midbins)):
    if Q_spec_completeness_correction[ii] == 0:     # don't overwrite cell once correction factor is computed
        if SF_midbins[ii] < Q_ratio_midbins[0]:      # extrapolate below lowest mass bin
            Q_spec_completeness_correction[ii] = m_Q[0]*SF_midbins[ii] + b_Q[0]    
        elif SF_midbins[ii] > Q_ratio_midbins[-1]:    # extrapolate above highest mass bin
            Q_spec_completeness_correction[ii] = m_Q[-1]*SF_midbins[ii] + b_Q[-1]    
        elif SF_midbins[ii] > Q_ratio_midbins[0] and SF_midbins[ii] < Q_ratio_midbins[-1]:    # interpolate in between all other points
            for jj in range(len(Q_ratio_midbins)-1):
                if SF_midbins[ii] > Q_ratio_midbins[jj] and SF_midbins[ii] < Q_ratio_midbins[jj+1]:
                    Q_spec_completeness_correction[ii] = m_Q[jj]*SF_midbins[ii] + b_Q[jj]
        else:
            print('Error in Q spec completeness correction computation. ABORT')
            break   
#
## Now compute errors on the interpolated correction factors
SF_spec_completeness_correction_err = interpolate_errors(SF_ratio,SF_ratio_err,SF_spec_completeness_correction,SF_ratio_midbins,SF_midbins)
Q_spec_completeness_correction_err = interpolate_errors(Q_ratio,Q_ratio_err,Q_spec_completeness_correction,Q_ratio_midbins,SF_midbins)
#
## FIGURE ##
#
if (plot_flag_1 == 1 and project_plot_flag ==2) or project_plot_flag == 1:
    if project_plot_flag == 0:
        pass
    else:
        # plot Spectroscopic completion correction factors 
        #plt.close()
        fig = plt.figure()
        string = 'Spec: %s'%z_cutoff[0]+'  Phot: %s'%z_cutoff[1]
        fig.suptitle(string, fontsize=30)
        ax = fig.add_subplot(1, 1, 1)
        ax.errorbar(SF_ratio_midbins,SF_ratio,yerr=SF_ratio_err, fmt='sb',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.8, mfc='none')
        ax.errorbar(Q_ratio_midbins,Q_ratio,yerr=Q_ratio_err, fmt='sr',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.8, mfc='none')
        ax.plot(SF_ratio_midbins,SF_ratio,'-b', linewidth=1.0, label='Star-forming')
        ax.plot(Q_ratio_midbins,Q_ratio,'-r', linewidth=1.0, label='Quiescent')
        ax.plot([6,13],[1,1],'--k',linewidth = 1.0)
        #ax.scatter(SF_midbins,SF_spec_completeness_correction,c='b', marker='x',linewidth=0.0)
        #ax.scatter(Q_midbins,Q_spec_completeness_correction,c='b', marker='x',linewidth=0.0)
        ax.errorbar(SF_midbins,SF_spec_completeness_correction,yerr=SF_spec_completeness_correction_err,c='b', marker='x',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.8, mfc='none')
        ax.errorbar(Q_midbins,Q_spec_completeness_correction,yerr=Q_spec_completeness_correction_err,c='r', marker='x',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.8, mfc='none')
        ax.legend(loc='upper right', frameon=False,fontsize=25)
        ax.set_xlabel('$log(M/M_{\odot})$',fontsize=30)
        ax.set_ylim=(-0.5,4.1)
        ax.set_ylabel('Correction factor C$_{s}$',fontsize=30)
        ax.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on', labelleft='on',labelsize=20)
        ax.minorticks_on()
        ax.grid(b=False)#, which='major', axis='both', color = 'k', linestyle = ':')
        ax.set_xlim((range2[0]-0.1),(range2[1]+0.1))
        #
#
# compute how many objects are added to each mass bin as a result of applying the spec_completeness_correction to the *_raw_smf lists. confirm that (# added to SF) + (# added to Q) = (# added to total)
SF_spec_completeness_diff = correction_difference(SF_phot_smf,np.transpose(SF_spec_completeness_correction))  
Q_spec_completeness_diff = correction_difference(Q_phot_smf,np.transpose(Q_spec_completeness_correction))
total_spec_completeness_diff = SF_spec_completeness_diff + Q_spec_completeness_diff
#
#
if (diag_flag_5 == 1 and diag_flag_master == 2) or diag_flag_master == 1:
    # Display correction factors
    print('\nSection 3.2: Spectroscopic completeness correction factors by bin (multiplicative): ')
    print('SF: ',str(np.transpose(SF_spec_completeness_correction)))
    print('Q: ',str(np.transpose(Q_spec_completeness_correction)),'\n')
    # Display some data for total, SF, Q: 
    print('Galaxies added due to SPECTROSCOPIC COMPLETENESS correction')
    print('SF: ',str(np.sum(SF_spec_completeness_diff,axis=0)),'\nTotal: %s'%np.sum(SF_spec_completeness_diff),'   or ',str((np.sum(SF_spec_completeness_diff)/np.sum(SF_raw_smf))*100),'%.\n')
    print('Q: ',str(np.sum(Q_spec_completeness_diff,axis=0)),'\nTotal: %s'%np.sum(Q_spec_completeness_diff),'   or ',str((np.sum(Q_spec_completeness_diff)/np.sum(Q_raw_smf))*100),'%.\n')
    print('Total: ',str(np.sum(total_spec_completeness_diff,axis=0)),'\nTotal: %s'%np.sum(total_spec_completeness_diff),'   or ',str((np.sum(total_spec_completeness_diff)/np.sum(total_raw_smf))*100),'%.\n')
else:
    print('\nSection 3.2: Galaxies added due to SPEC COMPLETENESS correction\nSF: %s'%np.sum(SF_spec_completeness_diff),'\nQ: %s'%np.sum(Q_spec_completeness_diff))
    #
#
#
#
#
## SECTION (3.3): NORMALIZE the SMF lists for each cluster by the TOTAL MASS in that cluster (i.e. integral under the SMF - so cluster*_smf x *_midbins); begin by adding the corrections just computed to the raw totals
#
## Add mass & spec corrections to *_raw_smf lists to form "*corrected_smf" lists
## CLUSTER
SF_corrected_smf = SF_raw_smf + SF_mass_completeness_diff + SF_spec_completeness_diff
Q_corrected_smf = Q_raw_smf + Q_mass_completeness_diff + Q_spec_completeness_diff
#
## Fields
SF_field_corrected_smf = SF_field_raw_smf + SF_field_mass_completeness_diff
Q_field_corrected_smf = Q_field_raw_smf + Q_field_mass_completeness_diff
SF_field_par_corrected_smf = SF_field_par_raw_smf + SF_field_par_mass_completeness_diff 
Q_field_par_corrected_smf = Q_field_par_raw_smf + Q_field_par_mass_completeness_diff
#
#
## The below calculation has been broken out into a separate function "normalize_smf_mass()" in SECTION (0). We maintain the below work for diagnostic information and summary tables, but we forgo this information for the "cluster field sample" and "parallel field sample", since we're only interested in the result obtained: a normalized SMF.
### NOTE: I tested the function versus the hard-coded algorithm explicitly. They are completely equivalent to within 11 decimal places, more than enough accuracy for the purposes of this calculation
#
## fix the shape of these arrays to be [6 clusters ,# of midbins in SMF], i.e. an array with 6 rows, each row containing an array with (#of midbins) data points
SF_corrected_smf = SF_corrected_smf.reshape((6,len(SF_midbins)))
Q_corrected_smf = Q_corrected_smf.reshape((6,len(SF_midbins)))
SF_field_corrected_smf = SF_field_corrected_smf.reshape((6,len(SF_midbins)))
Q_field_corrected_smf = Q_field_corrected_smf.reshape((6,len(SF_midbins)))
SF_field_par_corrected_smf = SF_field_par_corrected_smf.reshape((6,len(SF_midbins)))
Q_field_par_corrected_smf = Q_field_par_corrected_smf.reshape((6,len(SF_midbins)))
#total_corrected_smf = SF_corrected_smf + Q_corrected_smf 
#
N,M = SF_corrected_smf.shape      # store dimensions of raw_smf lists
#
if normalization_flag == 1:
    #
    ## compute the total mass (i.e. sum((# of galaxies in each bin)*(bin mass))
    total_mass = np.array([0]*6,dtype='float32')
    for ii in range(N):          # go through clusters 1 at a time
        total_mass[ii] = np.sum((SF_corrected_smf[ii]+Q_corrected_smf[ii])*np.transpose(SF_midbins))
    #
    ## now compute the raw relative mass fraction in each bin, for SF/Q by cluster; that is, create an array (row_1=SF, row_2=Q)
    SF_rel_mass = np.empty_like(SF_corrected_smf)
    Q_rel_mass = np.empty_like(Q_corrected_smf)
    for ii in range(N):          # go through clusters 1 at a time
        for jj in range(M):      # go through each mass bin one at a time
            SF_rel_mass[ii][jj] = ((SF_corrected_smf[ii][jj]*SF_midbins[jj]) / total_mass[ii])   # mass in the jj'th bin of cluster ii, divided by total mass of cluster ii
            Q_rel_mass[ii][jj] = ((Q_corrected_smf[ii][jj]*SF_midbins[jj]) / total_mass[ii]) # NOT AN ERROR: Q_midbins are offset by 0.05 for plotting purposes. the true value of the point is stored in SF_midbins
    #
    ## compute relative mass fractions of SF/Q by cluster after normalization
    mass_fraction_by_cluster = np.array([[0]*6]*2,dtype='float32')
    #
    mass_fraction_by_cluster[0] = np.sum(SF_rel_mass, axis=1)
    mass_fraction_by_cluster[1] = np.sum(Q_rel_mass, axis=1)
    #
    #
    #
    ### the "*_rel_mass" arrays now hold the relative amount of mass in each bin, where sum(SF_rel_mass) + sum(Q_rel_mass) = 1, i.e. the total mass of the cluster has been normalized to 1. Now sum across all 6 clusters, and divide by 6 s.t. the resulting array has a total area under its curve of 1 (not 6, which would seem rather arbitrary). 
    #
    # Mass in a bin is equal to: (# count in that bin) * (mass value of that bin). The above normalizes by mass, so the arrays "*_normalied_mass" contain the normalized amount of MASS IN EACH MASS BIN. To get the normalized # count, you need to divide (the normalized amount of mass in each bin) by (the mass value of that bin)
    ## Normalize SMFs by cluster. 
    SF_smf_by_cluster = np.empty_like(SF_corrected_smf)     #initialize arrays for FINAL SMF
    Q_smf_by_cluster = np.empty_like(Q_corrected_smf)
    ## normalize each cluster individually, for both SF & Q
    #
    for ii in range(N):          # go through clusters 1 at a time
        for jj in range(M):      # go through each mass bin (within each cluster) 1 at a time
            SF_smf_by_cluster[ii][jj] = (SF_rel_mass[ii][jj] / SF_midbins[jj])  # normalize each cluster by the TOTAL MASS IN EACH CLUSTER (hence the midbins multiplication); 
            Q_smf_by_cluster[ii][jj] = (Q_rel_mass[ii][jj] / SF_midbins[jj]) 
    #
    #
    #
    ## now we check that the mass fractions are still the same
    ## now compute the normalized amount of mass in each bin, for SF/Q by cluster; that is, create an array (row_1=SF, row_2=Q)
    ## compute the total mass (i.e. sum((# of galaxies in each bin)*(bin mass))
    total_norm_mass = np.array([0]*6,dtype='float32')
    for ii in range(N):          # go through clusters 1 at a time
        total_norm_mass[ii] = np.sum((SF_smf_by_cluster[ii]+Q_smf_by_cluster[ii])*np.transpose(SF_midbins))
    #
    SF_norm_mass = np.empty_like(SF_smf_by_cluster)
    Q_norm_mass = np.empty_like(Q_smf_by_cluster)
    for ii in range(N):      # go through each  cluster one at a time
        for jj in range(M):      # go through each mass bin (within each cluster) 1 at a time
            SF_norm_mass[ii][jj] = ((SF_smf_by_cluster[ii][jj]*SF_midbins[jj]) / total_norm_mass[ii])   # mass in the jj'th bin
            Q_norm_mass[ii][jj] = ((Q_smf_by_cluster[ii][jj]*SF_midbins[jj]) / total_norm_mass[ii]) # NOT AN ERROR: Q_midbins are offset by 0.05 for plotting purposes. the true value of the point is store in SF_midbins
    #
    ## compute relative mass fractions of SF/Q by cluster after normalization
    mass_fraction_by_cluster_norm = np.array([[0]*6]*2,dtype='float32')
    #
    mass_fraction_by_cluster_norm[0] = np.sum(SF_norm_mass, axis=1)
    mass_fraction_by_cluster_norm[1] = np.sum(Q_norm_mass, axis=1)
    #
    #
    ## for ease of future calling/plotting, combine all clusters into a sinlge SMF
    ## Clusters
    SF_smf = np.sum(SF_smf_by_cluster,axis=0)
    Q_smf = np.sum(Q_smf_by_cluster,axis=0)
    total_smf = SF_smf + Q_smf
    #
    #
    ###  CALL THE NORMALIZATION FUNCTION HERE
    SF_field_smf,Q_field_smf,total_field_smf = normalize_smf_mass(SF_field_corrected_smf,Q_field_corrected_smf,SF_midbins)
    SF_field_par_smf,Q_field_par_smf,total_field_par_smf = normalize_smf_mass(SF_field_par_corrected_smf,Q_field_par_corrected_smf,SF_midbins)
    ###
    #
    #
elif normalization_flag == 2:       # normalize by VOLUME
    #
    pass        # TO BE WRITTEN
    #
elif normalization_flag == 3:       # normalize by NUMBER COUNT
        #
    #
    ###  CALL THE NORMALIZATION FUNCTION HERE
    # cluster
    SF_smf,Q_smf,total_smf,count_fraction_by_cluster = normalize_smf_count(SF_corrected_smf,Q_corrected_smf,SF_midbins)           # fields               
    SF_field_smf,Q_field_smf,total_field_smf,count_fraction_by_cluster_field = normalize_smf_count(SF_field_corrected_smf,Q_field_corrected_smf,SF_midbins)
    SF_field_par_smf,Q_field_par_smf,total_field_par_smf,count_fraction_by_cluster_par = normalize_smf_count(SF_field_par_corrected_smf,Q_field_par_corrected_smf,SF_midbins)                             
    #
    #
#
## combine into a single SMF
SF_field_smf = (SF_field_smf + SF_field_par_smf)/2
Q_field_smf = (Q_field_smf + Q_field_par_smf)/2
total_field_smf = SF_field_smf + Q_field_smf
#
#
#
## diag_flag_6
# display diagnostics before/after normalization
#
if (diag_flag_6 == 1 and diag_flag_master == 2) or diag_flag_master == 1:
    if normalization_flag == 1:
        # the above should make the area under each of the SF/Q curves equal to 1 (by cluster), so the total area under the curve (sum of all clusters) should be 6. 
        print('\nSection 3.3: Normalization diagnostic:\nMASS fractions below are per cluster. The mass in each bin was calculated as: \n(count in mass bin)*(value of mass bin) = mass in that bin')
        print('\nTotal corrected raw samples\nSF: %s'%np.sum(SF_corrected_smf,axis=1),'\nQ: %s'%np.sum(Q_corrected_smf,axis=1))
        print('\nPre-Normalization:\nTotal SF fraction by mass in each cluster: \n%s'%mass_fraction_by_cluster[0],'\nTotal Q fraction by mass in each cluster: \n%s'%mass_fraction_by_cluster[1])
        print('Total relative mass per cluster (check by summing SF + Q from above): \n%s'%np.sum(mass_fraction_by_cluster,axis=0))
        print('\nPost-Normalization:\nTotal SF fraction by mass in each cluster: \n%s'%mass_fraction_by_cluster_norm[0])
        print('Total Q fraction by mass in each cluster: \n%s'%mass_fraction_by_cluster_norm[1])
        print('Total relative mass per cluster (check by summing SF + Q from above): \n%s'%np.sum(mass_fraction_by_cluster_norm,axis=0))
        #
        ## so we normalized the area under each cluster to 1, then summed all 6 clusters. so the total area under the curve should be 6. CHECK THAT
        SF_area = np.sum(SF_smf*np.transpose(SF_midbins))
        Q_area = np.sum(Q_smf*np.transpose(SF_midbins))
        print('\nCheck on final SMF - total area under curve should be 6\nArea under SF_smf: %s'%SF_area,'\nArea under Q_smf: %s'%Q_area,'\nArea under TOTAL curve: %s'%(SF_area+Q_area))
    #
    elif normalization_flag == 2:
        pass
    elif normalization_flag == 3:
        print('\nSection 3.3: Normalization diagnostic:\nNUMBER COUNT fractions below are per cluster. The final # count in each bin was calculated as: \n(count in mass bin) = (total count of galaxies in cluster) for each bin')
        print('\nTotal corrected raw samples (Pre-Normalization)\nSF: %s'%np.sum(SF_corrected_smf,axis=1),'\nQ: %s'%np.sum(Q_corrected_smf,axis=1))
        print('\nPost-Normalization:\nTotal SF fraction by # count in each cluster: \n%s'%count_fraction_by_cluster[0],'\nTotal Q fraction by # count in each cluster: \n%s'%count_fraction_by_cluster[1])
        print('Total relative count per cluster (check by summing SF + Q from above, s/b 1 in each cluster): \n%s'%np.sum(count_fraction_by_cluster,axis=0))
#
#
#
#
#
## SECTION (3.4) compute QUENCHED FRACTION
## defined as (# of quenched galaxies in a given mass bin) / (total # of galaxies in that mass bin)
#
quenched_fraction = np.array([[0]*len(SF_midbins)]*2,dtype='float32')    # row 1 = cluster;   row 2 = field
#
quenched_fraction[0] = Q_smf / total_smf                             # quenched fraction for cluster
quenched_fraction[1] = Q_field_smf / total_field_smf                 # quenched fraction for cluster
#
#
#
#
#
## SECTION (4) ERROR BARS & relative fractions
#
## METHOD 1: Poissonian error bars will be added to the spec. sample only; the error in the phot sample must propogate the error in the completeness correction; (total rel. error in smf, squared) = (sum of squared rel. errors (spec subsample = phot subsample);
## RECALL key lists: spec/phot subsamples are stored in "SF_spec_list/SF_phot_list" & "Q_spec_list/Q_phot_list", which are lists of lists organized by cluster, corresponding to "SF_spec_smf/SF_phot_smf" & "Q_spec_smf/Q_phot_smf" respectively
#
## Method 2: treat each bin as its own Poisson distribution, and draw samples (i.e. Bootstrap re-sampling with replacement), where the number of galaxies drawn in each realization of the resampling follows a Poisson distribution, with the mean of the distribution equal to the number count in the SMF mass bin. 
#
#
if smf_error_method == 1:                  # method 1: see above
    # initialize arrays
    SF_error_rel = np.empty_like(SF_midbins)
    Q_error_rel = np.empty_like(SF_midbins)
    SF_spec_error_rel = np.empty_like(SF_midbins)
    SF_phot_error_rel = np.empty_like(SF_midbins)
    SF_phot_error_counting_rel = np.empty_like(SF_midbins)
    Q_spec_error_rel = np.empty_like(SF_midbins)
    Q_phot_error_rel = np.empty_like(SF_midbins)
    Q_phot_error_counting_rel = np.empty_like(SF_midbins)
    #
    ## the error on the spec sample is just sqrt(N), i.e. a Poissonian error
    SF_spec_error_rel = np.sqrt(np.sum(SF_spec_smf,axis=0)) / np.sum(SF_spec_smf,axis=0)
    Q_spec_error_rel = np.sqrt(np.sum(Q_spec_smf,axis=0)) / np.sum(Q_spec_smf,axis=0)
    #
    ## reshape arrays to match SMF
    SF_spec_error_rel = SF_spec_error_rel.reshape(len(SF_midbins))
    Q_spec_error_rel = Q_spec_error_rel.reshape(len(SF_midbins))
    #
    ## remove NaNs due to empty bins
    for ii in range(len(SF_spec_error_rel)):
        if np.isnan(SF_spec_error_rel[ii])==1:
            SF_spec_error_rel[ii] = 0
        if np.isnan(Q_spec_error_rel[ii])==1:
            Q_spec_error_rel[ii] = 0
    #
    ## the error on the phot sample is also sqrt(N), plus an additional term for the uncertainty in the spec. completeness correction factor; this accounts for the counting error only still need to add in the uncertainty due to spec. completeness
    SF_phot_error_counting_rel = np.sqrt(np.sum(SF_phot_smf,axis=0))/np.sum(SF_phot_smf,axis=0)     
    Q_phot_error_counting_rel = np.sqrt(np.sum(Q_phot_smf,axis=0))/np.sum(Q_phot_smf,axis=0)       
    #
    ## reshape arrays to match SMF
    SF_phot_error_counting_rel = SF_phot_error_counting_rel.reshape(len(SF_midbins))
    Q_phot_error_counting_rel = Q_phot_error_counting_rel.reshape(len(SF_midbins))
    #
    ## remove NaNs due to empty bins
    for ii in range(len(SF_phot_error_counting_rel)):
        if np.isnan(SF_phot_error_counting_rel[ii])==1:
            SF_phot_error_counting_rel[ii] = 0
        if np.isnan(Q_phot_error_counting_rel[ii])==1:
            Q_phot_error_counting_rel[ii] = 0
    #
    ## compute the relative error on the phot sample
    SF_phot_error_rel = np.sqrt(np.transpose((SF_spec_completeness_correction_err/SF_spec_completeness_correction)**2) + (SF_phot_error_counting_rel)**2)
    Q_phot_error_rel = np.sqrt(np.transpose((Q_spec_completeness_correction_err/Q_spec_completeness_correction)**2) + (Q_phot_error_counting_rel)**2)
    #
    ## reshape arrays to match SMF
    SF_phot_error_rel = SF_phot_error_rel.reshape(len(SF_midbins))
    Q_phot_error_rel = Q_phot_error_rel.reshape(len(SF_midbins))
    #
    ## now put it all together: total relative error (squared) = sum of relative errors (squared)
    SF_error_rel = np.sqrt( SF_spec_error_rel**2 + SF_phot_error_rel**2 )
    Q_error_rel = np.sqrt( Q_spec_error_rel**2 + Q_phot_error_rel**2 )
    #
    ## compute error bars for SMF plot
    SF_error = SF_error_rel * SF_smf
    Q_error = Q_error_rel * Q_smf
    #
    total_error = SF_error + Q_error
    #
    ## reshape arrays to match SMF
    #SF_error = SF_error.reshape(len(SF_midbins))
    #Q_error = Q_error.reshape(len(SF_midbins))
    #total_error = total_error.reshape(len(SF_midbins))
#####
#
#
elif smf_error_method == 2:                # method 2: boostrap resampling
    #
    ## set SCALE of cluster to be sampled (i.e. # of galaxies in the 'typical' cluster for which you're bootstrapping)
    cluster_scale = 10000
    #
    ## set number of bootstrap re-samplings to be done
    num_bootstraps = 100
    #
    #
    #
    ## I need a list of all galaxies in each bin.
    ## RECALL: the list that stores all galaxy masses (before the histogram is created) is called "SF_list/Q_list", which stores galaxies by cluster. Collapse these lists into a single list, sort it by mass, and organize the galaxy masses into bins. The SMF list to follow is "SF_raw_smf/Q_raw_smf" (i.e. if the first bin in "SF_raw_smf" has 8 galaxies, than the first 8 galaxies in the ordered "SF_list" array are the members of that bin.
    #
    ## Step 1: flatten  & sort the SF/Q lists; and flatten SF/Q_raw_smf
    SF_list_flat = [item for sublist in SF_list for item in sublist]
    Q_list_flat = [item for sublist in Q_list for item in sublist]
    ## sort
    SF_list_flat = sorted(SF_list_flat)
    Q_list_flat = sorted(Q_list_flat)
    ## flatten SMFs
    SF_raw_smf = SF_raw_smf.reshape(6,12)
    SF_raw_smf_flat = np.sum(SF_raw_smf,axis=0)
    Q_raw_smf = Q_raw_smf.reshape(6,12)
    Q_raw_smf_flat = np.sum(Q_raw_smf,axis=0)
    #
    ## Step 2: create a list of lists, each sublist is the galaxies contained within that mass bin of the SMF
    SF_raw_smf_bins = [ [] for x in range(len(SF_raw_smf_flat))]      # initialize lists
    Q_raw_smf_bins = [ [] for x in range(len(Q_raw_smf_flat))]
    #
    ## now assign galaxy masses to each bin of the smf
    # SF
    for ii in range(len(SF_raw_smf_bins)):
        if ii == 0:
            index_start = 0
            index_end = SF_raw_smf_flat[ii]
        else:
            index_start = index_end
            index_end = index_end + SF_raw_smf_flat[ii]
        for jj in range(index_start,index_end):
            SF_raw_smf_bins[ii].append(SF_list_flat[jj])
    # Q
    for ii in range(len(Q_raw_smf_bins)):
        if ii == 0:
            index_start = 0
            index_end = Q_raw_smf_flat[ii]
        else:
            index_start = index_end
            index_end = index_end + Q_raw_smf_flat[ii]
        for jj in range(index_start,index_end):
            Q_raw_smf_bins[ii].append(Q_list_flat[jj])
    #
    ## List of galaxies in each bin DONE. Now need to do the bootstrap resampling (100 re-samples for each bin) where the number of galaxies drawn in each re-sampling follow a Poisson distribution with mean equal to the ***CORRECTED & NORMALIZED*** number count in each bin, scaled to a cluster w/ 1000 galaxies. 
    #
    ## first scale the FINALIZED smf from a cluster w/ 1 galaxy to a cluster w/ 1000
    SF_smf_bootstrap = SF_smf * cluster_scale
    Q_smf_bootstrap = Q_smf * cluster_scale
    #
    ## Construct a loop that re-samples the SMF bins one bin at a time. Each bin will be sampled 100 times, with the number of galaxies drawn varying with each re-sample; the distribution of galaxies drawn will follow a Poisson distribution with mean equal to the number count of galaxies for that SMF mass bin
    # SF
    SF_error_bootstrap = np.empty_like(SF_midbins)
    SF_other_bootstrap = np.array([[0]*2]*len(SF_midbins),dtype='float32')
    #
    for ii in range(len(SF_raw_smf_bins)):
        bootstrap_array = np.array(SF_raw_smf_bins[ii])
        mean = SF_smf_bootstrap[ii]           # this is the mean of the Poisson distribution for how many galaxies to draw
        if mean == 0:
            SF_error_bootstrap[ii] = 0
        else:
            num_to_draw = np.random.poisson(mean, num_bootstraps)
            #
            ## initialize an array to store the result of bootstrapping, & a bootstrapping function for the statistic you want returned
            bootstrap_result = np.array([[0.0]*3]*num_bootstraps)    # col1= # of gal drawn; col2 = mean; col3= std dev
            bootstrap_statistic = lambda x: (np.mean(x), np.std(x))
            #
            ## now loop through "num_to_draw" one at a time, each time doing a single bootstrap, and storing the result
            ## add error-handling for poisson samples which are zero - pass; store zeros as result
            for jj in range(len(num_to_draw)):
                if num_to_draw[jj] == 0:
                    bootstrap_result[jj] = [0,0,0]
                else:
                    boot = bootstrap(bootstrap_array, 1, samples = num_to_draw[jj], bootfunc=bootstrap_statistic)
                    boot = boot.reshape(2)
                    bootstrap_result[jj][0] = num_to_draw[jj]
                    bootstrap_result[jj][1] = boot[0]
                    bootstrap_result[jj][2] = boot[-1]
            # we now have the results of 100 bootstrap re-samples for a single SMF mass bin. save the mean std. dev. as the error on that SMF mass bin point (i.e. the errorbars on the y-values of the SMF). 
            bootstrap_means = np.mean(bootstrap_result,axis=0)
            SF_error_bootstrap[ii] = bootstrap_means[-1]         # this is the ERROR
            SF_other_bootstrap[ii][0] = bootstrap_means[0]       # this is the # of galaxies drawn in the bootstrap re-sample
            SF_other_bootstrap[ii][1] = bootstrap_means[1]       # this is the mean drawn
        #
    SF_error_bootstrap = SF_error_bootstrap.reshape(12)
    #
    # Q
    Q_error_bootstrap = np.empty_like(SF_midbins)
    Q_other_bootstrap = np.array([[0]*2]*len(SF_midbins),dtype='float32')    #
    for ii in range(len(Q_raw_smf_bins)):
        bootstrap_array = np.array(Q_raw_smf_bins[ii])
        mean = Q_smf_bootstrap[ii]           # this is the mean of the Poisson distribution for how many galaxies to draw
        if mean == 0:
            Q_error_bootstrap[ii] = 0
        else:
            num_to_draw = np.random.poisson(mean, num_bootstraps)
            #
            ## initialize an array to store the result of bootstrapping, & a bootstrapping function for the statistic you want returned
            bootstrap_result = np.array([[0.0]*3]*num_bootstraps)    # col1= # of gal drawn; col2 = mean; col3= std dev
            bootstrap_statistic = lambda x: (np.mean(x), np.std(x))
            #
            ## now loop through "num_to_draw" one at a time, each time doing a single bootstrap, and storing the result
            ## add error-handling for poisson samples which are zero - pass; store zeros as result
            for jj in range(len(num_to_draw)):
                if num_to_draw[jj] == 0:
                    bootstrap_result[jj] = [0,0,0]
                else:
                    boot = bootstrap(bootstrap_array, 1, samples = num_to_draw[jj], bootfunc=bootstrap_statistic)
                    boot = boot.reshape(2)
                    bootstrap_result[jj][0] = num_to_draw[jj]
                    bootstrap_result[jj][1] = boot[0]
                    bootstrap_result[jj][2] = boot[-1]
            # we now have the results of 100 bootstrap re-samples for a single SMF mass bin. save the mean std. dev. as the error on that SMF mass bin point (i.e. the errorbars on the y-values of the SMF). 
            bootstrap_means = np.mean(bootstrap_result,axis=0)
            Q_error_bootstrap[ii] = bootstrap_means[-1]         # this is the ERROR
            Q_other_bootstrap[ii][0] = bootstrap_means[0]       # this is the # of galaxies drawn in the bootstrap re-sample
            Q_other_bootstrap[ii][1] = bootstrap_means[1]       # this is the mean drawn
        #
    Q_error_bootstrap = Q_error_bootstrap.reshape(12)
    #
    ## I realized after the fact that I don't actually want the below output. What I'm really interested in is the relative error found from bootstrapping. This may still be useful diagnostic info tho, so change the criteria for the diag_flag if-statement from ==1 to ==99.
    if (diag_flag_4 == 99 and diag_flag_master == 2) or diag_flag_master == 1:
        ## We now have errors estimated on an SMF containing 10,000 galaxies. display this information
        print('\n(scaled) SF SMF for bootstrapping: \n%s'%SF_smf_bootstrap)
        print('\nSF error from bootstrapping: \n%s'%SF_error_bootstrap)
        print('\nSF: [# count drawn in bin, mean mass]: \n%s'%SF_other_bootstrap)
        print('\n(scaled) Q SMF for bootstrapping: \n%s'%Q_smf_bootstrap)
        print('\nQ error from bootstrapping: \n%s'%Q_error_bootstrap)
        print('\nQ: [# count drawn in bin, mean mass]: \n%s'%Q_other_bootstrap)
        
    #
    #
    ## Now compute the relative error, which will be applied to the SMF of arbitrary normalization
    SF_rel_error = SF_error_bootstrap / SF_smf_bootstrap
    Q_rel_error = Q_error_bootstrap / Q_smf_bootstrap
    #
    if (diag_flag_4 == 1 and diag_flag_master == 2) or diag_flag_master == 1:
        print('REL SF error from bootstrapping: \n%s'%SF_rel_error)
        print('REL Q error from bootstrapping: \n%s'%Q_rel_error)
    #
    ## compute the error bars on the SMF mass bins
    SF_error = SF_rel_error * SF_smf
    Q_error = Q_rel_error * Q_smf
    total_error = SF_error + Q_error
    #
elif smf_error_method == 3:                  # method 3: treat entire SMF as one sample for counting errors
    # initialize arrays
    SF_error_rel = np.empty_like(SF_midbins)
    Q_error_rel = np.empty_like(SF_midbins)
    SF_error = np.empty_like(SF_midbins)
    Q__error = np.empty_like(SF_midbins)
    SF_error_counting_rel = np.empty_like(SF_midbins)
    Q_error_counting_rel = np.empty_like(SF_midbins)
    #
    ## the COUNTING error on the RAW SMF is just sqrt(N), i.e. a Poissonian error
    SF_error_counting_rel = np.sqrt(np.sum(SF_raw_smf,axis=0)) / np.sum(SF_raw_smf,axis=0)
    Q_error_counting_rel = np.sqrt(np.sum(Q_raw_smf,axis=0)) / np.sum(Q_raw_smf,axis=0)
    #
    ## reshape arrays to match SMF
    SF_error_counting_rel = SF_error_counting_rel.reshape(len(SF_midbins))
    Q_error_counting_rel = Q_error_counting_rel.reshape(len(SF_midbins))
    #
    ## remove NaNs due to empty bins
    for ii in range(len(SF_error_counting_rel)):
        if np.isnan(SF_error_counting_rel[ii])==1:
            SF_error_counting_rel[ii] = 0
        if np.isnan(Q_error_counting_rel[ii])==1:
            Q_error_counting_rel[ii] = 0
    #
    ## compute the relative error on the SMF
    SF_error_rel = np.sqrt(np.transpose((SF_spec_completeness_correction_err/SF_spec_completeness_correction)**2) + (SF_error_counting_rel)**2)
    Q_error_rel = np.sqrt(np.transpose((Q_spec_completeness_correction_err/Q_spec_completeness_correction)**2) + (Q_error_counting_rel)**2)
    #
    ## reshape arrays to match SMF
    SF_error_rel = SF_error_rel.reshape(len(SF_midbins))
    Q_error_rel = Q_error_rel.reshape(len(SF_midbins))
    #
    ## compute error bars for SMF plot
    SF_error = SF_error_rel * SF_smf
    Q_error = Q_error_rel * Q_smf
    #
    total_error = SF_error + Q_error
    #
    ## reshape arrays to match SMF
    #SF_error = SF_error.reshape(len(SF_midbins))
    #Q_error = Q_error.reshape(len(SF_midbins))
    #total_error = total_error.reshape(len(SF_midbins))
#####
#
#
#
#
#
#
#
#
if (plot_flag_2 == 1 and project_plot_flag ==2) or project_plot_flag == 1: # plot interpolated/extrapolated points on top of computed correction fractions
    if project_plot_flag == 0:
        pass
    else:
        SF_field_error = np.zeros_like(SF_smf)
        Q_field_error = np.zeros_like(Q_smf)
        total_field_error = np.zeros_like(total_smf)
        quenched_err = np.zeros_like(quenched_fraction[0])
    ## upper: SMF for cluster, field;       lower: fractions of SF/Q in cluster, field
    #
        #plt.close()
        SMF = plt.figure()
        string = 'Spec: %s'%z_cutoff[0]+'  Phot: %s'%z_cutoff[1]+'  Method: %i'%smf_error_method
        SMF.suptitle(string, fontsize=30)
        gs = gridspec.GridSpec(2,2, wspace=0, hspace=0, width_ratios=[1,1], height_ratios=[2,1])   #make a tiled-plot like vdB2013 w/ fractions below, this line sets the proporitons of plots in the figure
        #gs = gridspec.GridSpec(2,3, width_ratios=[1,1,1], height_ratios=[2,1])   #make a tiled-plot like vdB2013 w/ fractions below, this line sets the proporitons of plots in the figure
        #
        ## CLUSTER
        ax0 = plt.subplot(gs[0])      
        ax0.errorbar(SF_midbins,SF_smf,yerr=SF_error, fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5, label='Star-forming')#yerr=SF_error,
        ax0.errorbar(Q_midbins,Q_smf,yerr=Q_error,fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5,label='Quiescent')
        ax0.errorbar(SF_midbins,total_smf,yerr=total_error,fmt='.k',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5,label='Total')
        ## Plot Schechter fits:  (uncomment 5 hashtags when fits complete)
        ######plt.plot(x_plot_Q,Q_model_ml_plot, ':r')
        #####plt.plot(x_plot_Q,Q_model_mcmc_plot, '--r')
        ######plt.plot(x_plot_SF,SF_model_ml_plot, ':c', label = 'Max. Likelihood', linewidth = 0.5)
        #####plt.plot(x_plot_SF,SF_model_mcmc_plot, '--b', label = 'MCMC', linewidth = 0.5)
        #####plt.plot(x_plot_T,T_model_mcmc_plot, 'k')
        ax0.set_xlabel('$log(M/M_{\odot})$')
        ax0.set_xscale('linear')
        ax0.minorticks_on()
        ax0.set_xlim(7,12.5)
        ax0.set_yscale('log')
        ax0.set_ylim(6e-5,0.5)
        ax0.minorticks_on()
        ax0.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=False,labelbottom=False,grid_alpha=0.4,grid_linestyle=':')
        ax0.yaxis.set_label_position("left")
        ax0.set_ylabel('???')
        ax0.set_title('Cluster')
        ax0.legend(scatterpoints=1,loc='lower left', frameon=False, fontsize = 'x-small')
        ax0.grid(b=False)#, which='major', axis='both', color = 'k', linestyle = '--')
        #
        ## cluster fraction
        ax2 = plt.subplot(gs[2])    
        #plt.plot(SF_midbins,frac_smf[0],'.b',linewidth=0.5)
        #plt.plot(SF_midbins,frac_smf[1],'.r',linewidth=0.5)
        ax2.errorbar(SF_midbins,quenched_fraction[0],yerr=quenched_err, fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
        #ax2.errorbar(SF_midbins,frac_smf[1],yerr=frac_error[1], fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
        ax2.set_xscale('linear')
        ax2.set_xlabel('$log(M/M_{\odot})$')
        ax2.set_xlim(7,12.5)
        ax2.set_yscale('linear')
        ax2.set_ylim(-0.1,1.1)
        ax2.minorticks_on()
        ax2.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=False)
        ax2.yaxis.set_label_position("left")
        ax2.set_ylabel('Quenched fraction')
        #
        ## FIELD
        ax1 = plt.subplot(gs[1])      
        ax1.errorbar(SF_midbins,SF_field_smf,yerr=SF_field_error, fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5, label='Star-forming')#yerr=SF_error,
        ax1.errorbar(Q_midbins,Q_field_smf,yerr=Q_field_error,fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5,label='Quiescent')
        ax1.errorbar(SF_midbins,total_field_smf,yerr=total_field_error,fmt='.k',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5,label='Total')
        ## Plot Schechter fits:  (uncomment 5 hashtags when fits complete)
        ######plt.plot(x_plot_Q,Q_model_ml_plot, ':r')
        #####plt.plot(x_plot_Q,Q_model_mcmc_plot, '--r')
        ######plt.plot(x_plot_SF,SF_model_ml_plot, ':c', label = 'Max. Likelihood', linewidth = 0.5)
        #####plt.plot(x_plot_SF,SF_model_mcmc_plot, '--b', label = 'MCMC', linewidth = 0.5)
        #####plt.plot(x_plot_T,T_model_mcmc_plot, 'k')
        ax1.set_xlabel('$log(M/M_{\odot})$')
        ax1.set_xscale('linear')
        ax1.minorticks_on()
        ax1.set_xlim(7,12.5)
        ax1.set_yscale('log')
        ax1.set_ylim(6e-5,0.5)
        ax1.minorticks_on()
        ax1.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=True,labelleft=False,labelbottom=False,grid_alpha=0.4,grid_linestyle=':')
        ax1.yaxis.set_label_position("right")
        ax1.set_ylabel('???')
        ax1.set_title('Field')
        #ax3.legend(scatterpoints=1,loc='lower left', frameon=False, fontsize = 'x-small')
        ax1.grid(b=False)#, which='major', axis='both', color = 'k', linestyle = '--')
        #
        ## field fraction
        ax3 = plt.subplot(gs[3])    
        ax3.errorbar(SF_midbins,quenched_fraction[1],yerr=quenched_err, fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
        #ax2.errorbar(SF_midbins,frac_smf[1],yerr=frac_error[1], fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
        ax3.set_xscale('linear')
        ax3.set_xlabel('$log(M/M_{\odot})$')
        ax3.set_xlim(7,12.5)
        ax3.set_yscale('linear')
        ax3.set_ylim(-0.1,1.1)
        ax3.minorticks_on()
        ax3.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=True)
        ax3.yaxis.set_label_position("right")
        ax3.set_ylabel('Quenched fraction')
        #
        plt.show()
#
#
#
#
#
#
#
#
################
#                   THIS IS WHERE I'M AT
################
#
#
#
#
if where_im_at_flag == 1:
    #
    ########################
    ########################
    #
    #
    # BREAK    
    #
    #
    ########################
    ########################
    #
    #
    #
    ## SECTION (5)    EMCEE simulation; see emcee_chi2_final.py;
    #
    ######## This section has been broken out into its own program, called "emcee_chi2", which uses chi-squared as the cost function and is the code to be implemented in the final run for this project 
    #
    #
    #
    ## The following summarizes the result of the MCMC simulation and sets up the appropriate arrays for plotting
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
    #
    #    
    #
    #
    #
    ## SECTION (6): build SCHECHTER MODELS of MCMC fits
    #
    ## define x array to generate points to plot Schechter fit
    x = np.linspace(SF_midbins[0],SF_midbins[len(SF_midbins)-1],num=1000)#
    #
    #
    ## build a model for plotting
    SF_model_mcmc = np.log(10)*SFphi_mcmc*(10**((x-SFM_star_mcmc)*(1+SFalpha_mcmc)))*np.exp(-10**(x-SFM_star_mcmc))
    ## single-schechter Q pop
    Q_model_mcmc = np.log(10)*Qphi_mcmc*(10**((x-QM_star_mcmc)*(1+Qalpha_mcmc)))*np.exp(-10**(x-QM_star_mcmc))
    ## double-schechter Q pop
    #Q_model_mcmc = np.log(10)*np.exp(-10**(x-QM_star_mcmc))*((Qphi1_mcmc*(10**((x-QM_star_mcmc)*(1+Qalpha1_mcmc))))+(Qphi2_mcmc*(10**((x-QM_star_mcmc)*(1+Qalpha2_mcmc)))))
    T_model_mcmc = np.log(10)*Tphi_mcmc*(10**((x-TM_star_mcmc)*(1+Talpha_mcmc)))*np.exp(-10**(x-TM_star_mcmc))
    #
    ## create arrays for plotting purposes to display values down to a y_min = 1
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
    #
    #
    #
    ## SECTION (7): create PLOTS for SMF
    #
    #
    ### MAY NEED TO EDIT: plot_flag
    ## upper: SMF for cluster, field;       lower: fractions of SF/Q in cluster, field
    #
    plt.close()
    SMF = plt.figure(num=1)
    smf = gridspec.GridSpec(2,2, wspace=0, hspace=0, width_ratios=[1,1], height_ratios=[2,1])   #make a tiled-plot like vdB2013 w/ fractions below, this line sets the proporitons of plots in the figure
    #gs = gridspec.GridSpec(2,3, width_ratios=[1,1,1], height_ratios=[2,1])   #make a tiled-plot like vdB2013 w/ fractions below, this line sets the proporitons of plots in the figure
    # Cluster
    cluster = plt.subplot(gs[0])      
    xa = plt.errorbar(SF_midbins,SF_smf,yerr=SF_error, fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)#yerr=SF_error,
    xb = plt.errorbar(Q_midbins,Q_smf,yerr=Q_error,fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    xc = plt.errorbar(SF_midbins,total_smf,yerr=total_error,fmt='.k',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    ## Plot Schechter fits:  (uncomment 5 hashtags when fits complete)
    ######plt.plot(x_plot_Q,Q_model_ml_plot, ':r')
    #####plt.plot(x_plot_Q,Q_model_mcmc_plot, '--r')
    ######plt.plot(x_plot_SF,SF_model_ml_plot, ':c', label = 'Max. Likelihood', linewidth = 0.5)
    #####plt.plot(x_plot_SF,SF_model_mcmc_plot, '--b', label = 'MCMC', linewidth = 0.5)
    #####plt.plot(x_plot_T,T_model_mcmc_plot, 'k')
    plt.xscale('linear')
    plt.xlim=(8,12.25)
    plt.yscale('log')
    plt.ylim=(1,1e3)
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
    ## SECTION (????): compare CLUSTERvFIELD BY POPULATION
    #
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

