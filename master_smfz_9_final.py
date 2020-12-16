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
### (4)    add ERROR BARS for scatter plot;
### (4.1)  compute QUENCHED FRACTION (i.e. # Q galaxies / (total # SF+Q galaxies) in each mass bin
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
import math
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
            if range2_flag == 0 or range2_flag == 1:
                if limiting_masses[jj] <= mass_bin_edges[ii]:    # count # of clusters complete at each mass bin
                    mass_completeness_correction_factors[ii]+=1
            elif range2_flag == 2:
                if limiting_masses[jj] <= (mass_bin_edges[ii] + 0.06):    # count # of clusters complete at each mass bin
                    mass_completeness_correction_factors[ii]+=1
                # if limiting_masses[jj] <= (mass_bin_edges[ii]):# + 0.06):    # count # of clusters complete at each mass bin
                #     mass_completeness_correction_factors[ii]+=1
                # elif limiting_masses[jj] > (mass_bin_edges[ii]) and limiting_masses[jj] <= (mass_bin_edges[ii+1]):# + 0.06):    # count # of clusters complete at each mass bin
                #     relative_completeness = ( bin_width - (limiting_masses[jj] - mass_bin_edges[ii]) ) / bin_width
                #     mass_completeness_correction_factors[ii]+=relative_completeness
                #
    # print(mass_completeness_correction_factors)
    mass_completeness_correction_factors = np.transpose(len(limiting_masses)/mass_completeness_correction_factors) # take recipricol for multiplicative factors
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
### TEMPORARY WRITING FLAG - comments out all code past the point where the flag is invoked; should be set to ==0
where_im_at_flag = 0
#
#
## MAY NEED TO EDIT: NORMALIZATION FLAG
## Choose your normalization:   1=by mass;  2=by volume;   3=by number count (i.e. # of galaxies in each cluster); FINAL DECISION: flag == 3 (by # count)
normalization_flag = 3
#
#
## MAY NEED TO EDIT: choose the method to CALCULATE ERRORS
smf_error_method = 1             #   1 = propagate errors;   2 = bootstrap resampling; FINAL DECISION: flag == 1
spec_completeness_error_flag = 2     # 1 = sum of rel. errors (DEPRECATED - incorrect); 2 = (rel. error of numerator)^2 + (rel. error of denom.)^2 (THIS FLAG SHOULD ALWAYS BE SET TO =2 )
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
plot_flag_1 = 0            # spec completeness correction factors
plot_flag_2 = 0           # SMF
cluster_only_plot_flag = 0  # create a plot of just the cluster data alone
field_only_plot_flag = 0  # create a plot of just the field data alone
field_construction_flag = 0 # create plots of HFF field and UVC field separately, and then the combined plot
UVC_fit_flag = 0            # create plot fitting just the UVC sample to a schechter function
schechter_plot_flag = 1     # plot the best-fitting schechter function: 0==off;  1==on
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
    limiting_mass = [7.62,7.63,8.2,7.5,7.15,7.34] # clusters 1,2,3,4,5,6, see IDs below;
#
SF_pos_lost = np.array([0]*6)        # to track SF/Q false pos/neg objects lost due to their being below the mass limit, by cluster
SF_neg_lost = np.array([0]*6)
Q_pos_lost = np.array([0]*6)
Q_neg_lost = np.array([0]*6)
other_lost = np.array([0]*6)    #objects below limiting mass other than false pos/neg
below_range2 = np.array([[0]*6]*2)
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
                    if master_cat['lmass'][counter] >= range2[0]:
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
                    elif master_cat['lmass'][counter] < range2[0]:
                        if master_cat['type'][counter] == 1:
                            below_range2[0][cluster]+=1                    # SF neg
                        elif master_cat['type'][counter] ==2:
                            below_range2[1][cluster]+=1
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
                #
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
    member_smf_names = Column(['TOTAL Members (master_data*.py)','Total SF >lim. mass','SF - Phot.','SF - Spec.','SF <lim. mass','SF < range2','Total Q >lim. mass','Q - Phot.','Q - Spec.','Q <lim. mass','Q < range2','SUM'],name='Property')
    col_names = cluster_names
    # SF table
    member_smf0 = Column([np.sum([mem_phot,mem_spec]),num_SF,num_SF_phot,num_SF_spec,num_SF_lost,np.sum(below_range2[0]),num_Q,num_Q_phot,num_Q_spec,num_Q_lost,np.sum(below_range2[1]),np.sum([num_SF_phot,num_SF_spec,num_Q_phot,num_Q_spec,num_SF_lost,num_Q_lost,np.sum(below_range2[0]),np.sum(below_range2[1])])],name='Total')  # total column
    member_smf_stats = Table([member_smf_names,member_smf0])
    for ii in range(len(mem_spec[0])):
        col = Column([np.sum([mem_phot[0][ii],mem_phot[1][ii],mem_spec[0][ii],mem_spec[1][ii]]),len(SF_list[ii]),len(SF_phot_list[ii]),len(SF_spec_list[ii]),len(SF_lost[ii]),below_range2[0][ii],len(Q_list[ii]),len(Q_phot_list[ii]),len(Q_spec_list[ii]),len(Q_lost[ii]),below_range2[1][ii],np.sum([len(SF_phot_list[ii]),len(SF_spec_list[ii]),len(SF_lost[ii]),len(Q_phot_list[ii]),len(Q_spec_list[ii]),len(Q_lost[ii]),below_range2[0][ii],below_range2[1][ii]])],name=col_names[ii])
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
##### use "SF_field_uvc_list/Q_field_uvc_list" from 'UVC_master_data.py'
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
if range2_flag == 0:
    range2 = [range2[0],max(max_SF,max_Q)]
elif range2_flag != 0:
    pass
# range2 = [7.3,12.2]             # bin edges hard-coded. lowest/hightest mass objects are 7.36/12.18
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
print('SF: ',str(np.sum(SF_raw_smf,axis=0)))
print('Q: ',str(np.sum(Q_raw_smf,axis=0)))
print('Total: ',str(np.sum(total_raw_smf)),'\n')
print('\nRAW totals - Field\nSF field (clu): %s'%np.sum(SF_field_raw_smf,axis=0),'\nSF field (par): %s'%np.sum(SF_field_par_raw_smf,axis=0))
print('Q field (clu): %s'%np.sum(Q_field_raw_smf,axis=0),'\nQ field (par): %s'%np.sum(Q_field_par_raw_smf,axis=0))
print('Total (clu): %s'%np.sum(total_field_raw_smf,axis=0),'\nTotal (par): %s'%np.sum(total_field_par_raw_smf,axis=0))
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
print('\nSF_midbins: %s'%SF_midbins+'\n')
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
    print('\nSection 2.1: Differences between raw cluster count and (spec + phot) subsamples, by cluster')
    for ii in range(len(SF_raw_smf)):
        SF_diff[ii] = SF_raw_smf[ii] - (SF_phot_smf[ii] + SF_spec_smf[ii])
        Q_diff[ii] = Q_raw_smf[ii] - (Q_phot_smf[ii] + Q_spec_smf[ii])
        print('SF',str(ii+1),' difference: ',str(np.sum(SF_diff[ii])))
        print('Q',str(ii+1),' difference: ',str(np.sum(Q_diff[ii])))
        print('Total difference: %s'%np.sum([SF_diff,Q_diff]),'\n')
    print('\nSF spec completeness: %s'%(np.sum(SF_spec_smf,axis=0)/np.sum(SF_raw_smf,axis=0)))
    print('Q spec completeness: %s'%(np.sum(Q_spec_smf,axis=0)/np.sum(Q_raw_smf,axis=0))+'\n')
    print('NOTE: discrepancy in totals for Table 1B due to galaxies above limiting mass of cluster, but below limiting mass of study.')
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
##
## cluster  limiting mass corrections
mass_completeness_correction = mass_completeness_correction_function(mass_bins,limiting_mass)   # function is the same as steps taken above, but now in function form
#
#
## field limiting mass corrections: UPDATE
## the field sample is taken across all 12 images (6 cluster fields, 6 parallel fields). The completeness correction then needs to be applied to the full HFF field SMF (i.e. the cluster_field_smf + parallel_field_smf), and computed on the basis of how many images of the 12 are complete at each mass bin
#
## first, compute the mass_completeness_correction_field factors. you need to combine limiting_mass array with limiting_mass_par array into a single array of length 12
if cluster_field_inclusion_flag == 0:
    limiting_mass_field = limiting_mass_par
elif cluster_field_inclusion_flag == 1:
    limiting_mass_field = np.concatenate([limiting_mass_cluster_field,limiting_mass_par])
#
mass_completeness_correction_field = mass_completeness_correction_function(mass_bins,limiting_mass_field)
## first, combine the cluster_field_smf + parallel_field_smf to obtain the full HFF RAW smf
SF_field_raw_smf_HFF = SF_field_raw_smf + SF_field_par_raw_smf
Q_field_raw_smf_HFF = Q_field_raw_smf + Q_field_par_raw_smf
total_field_raw_smf_HFF = SF_field_raw_smf_HFF + Q_field_raw_smf_HFF
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
SF_field_mass_completeness_diff = correction_difference(SF_field_raw_smf_HFF,mass_completeness_correction_field)
Q_field_mass_completeness_diff = correction_difference(Q_field_raw_smf_HFF,mass_completeness_correction_field)
total_field_mass_completeness_diff = SF_field_mass_completeness_diff + Q_field_mass_completeness_diff
#
### MAY NEED TO EDIT: diag_flag_3
#
if (diag_flag_3 == 1 and diag_flag_master == 2) or diag_flag_master == 1:
# Display correction factors
    print('/nSection 3.1: Mass completeness correction factors by bin\nCLUSTER SMF: ',str(mass_completeness_correction),'\nFIELD SMF: %s'%mass_completeness_correction_field)
    #
    # Display some data for total, SF, Q:
    print('Section 3.1: Galaxies added due to MASS COMPLETENESS correction - CLUSTER')
    print('CLUSTER SF: ',str(np.sum(SF_mass_completeness_diff,axis=0)),'\nTotal: %s'%np.sum(SF_mass_completeness_diff),', or ',str((np.sum(SF_mass_completeness_diff)/np.sum(SF_raw_smf))*100),'%\n')
    print('CLUSTER Q: ',str(np.sum(Q_mass_completeness_diff,axis=0)),'\nTotal: %s'%np.sum(Q_mass_completeness_diff),', or ',str((np.sum(Q_mass_completeness_diff)/np.sum(Q_raw_smf))*100),'%\n')
    print('CLUSTER Total: ',str(np.sum(total_mass_completeness_diff,axis=0)),'\nTotal: %s'%np.sum(total_mass_completeness_diff),', or ',str((np.sum(total_mass_completeness_diff)/np.sum(total_raw_smf))*100),'%\n')
    #
    # Display some data for total, SF, Q:
    print('Section 3.1: Galaxies added due to MASS COMPLETENESS correction - FIELD')
    print('FIELD SF: ',str(np.sum(SF_field_mass_completeness_diff,axis=0)),'\nTotal: %s'%np.sum(SF_field_mass_completeness_diff),', or ',str((np.sum(SF_field_mass_completeness_diff)/np.sum(SF_field_raw_smf_HFF))*100),'%\n')
    print('FIELD Q: ',str(np.sum(Q_field_mass_completeness_diff,axis=0)),'\nTotal: %s'%np.sum(Q_field_mass_completeness_diff),', or ',str((np.sum(Q_field_mass_completeness_diff)/np.sum(Q_field_raw_smf_HFF))*100),'%\n')
    print('FIELD Total: ',str(np.sum(total_field_mass_completeness_diff,axis=0)),'\nTotal: %s'%np.sum(total_field_mass_completeness_diff),', or ',str((np.sum(total_field_mass_completeness_diff)/np.sum(total_field_raw_smf_HFF))*100),'%\n')
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
#
## Now compute errors on the interpolated correction factors
SF_spec_completeness_correction_err = interpolate_errors(SF_ratio,SF_ratio_err,SF_spec_completeness_correction,SF_ratio_midbins,SF_midbins)
Q_spec_completeness_correction_err = interpolate_errors(Q_ratio,Q_ratio_err,Q_spec_completeness_correction,Q_ratio_midbins,SF_midbins)
#
#
## do NOT apply membership correction to 1st (up to two) bins due to spectroscopic incompleteness
if bins_exempt_from_membership_correction == 4:
    for ii in range(len(SF_spec_completeness_correction)):
        SF_spec_completeness_correction[ii] = 1
        Q_spec_completeness_correction[ii] = 1
        SF_spec_completeness_correction_err[ii] = 0
        Q_spec_completeness_correction_err[ii] = 0
elif bins_exempt_from_membership_correction == 0:
    pass
else:
    for ii in range(bins_exempt_from_membership_correction):
        SF_spec_completeness_correction[ii] = 1
        Q_spec_completeness_correction[ii] = 1
        SF_spec_completeness_correction_err[ii] = 0
        Q_spec_completeness_correction_err[ii] = 0
#
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
        string = 'Spec: %s'%z_cutoff[0]+'  Phot: %s'%z_cutoff[1]+'  Method: %i'%membership_correction_binning_flag
        fig.suptitle(string, fontsize=30)
        ax = fig.add_subplot(1, 1, 1)
        ax.errorbar(SF_ratio_midbins,SF_ratio,yerr=SF_ratio_err, fmt='sb',lolims=False, uplims=False, linewidth=0.0, elinewidth=2, mfc='b', ms=20)
        ax.errorbar(Q_ratio_midbins,Q_ratio,yerr=Q_ratio_err, fmt='sr',lolims=False, uplims=False, linewidth=0.0, elinewidth=2, mfc='r', ms=20)
        ax.plot(SF_ratio_midbins,SF_ratio,'-b', linewidth=2.0, label='Star-forming')
        ax.plot(Q_ratio_midbins,Q_ratio,'-r', linewidth=2.0, label='Quiescent')
        ax.plot([6,13],[1,1],'--k',linewidth = 2.0)
        #ax.scatter(SF_midbins,SF_spec_completeness_correction,c='b', marker='x',linewidth=0.0)
        #ax.scatter(Q_midbins,Q_spec_completeness_correction,c='b', marker='x',linewidth=0.0)
        SF_spec_completeness_correction_err = SF_spec_completeness_correction_err.reshape(len(SF_midbins),)
        Q_spec_completeness_correction_err = Q_spec_completeness_correction_err.reshape(len(SF_midbins),)
        ax.errorbar(SF_midbins,SF_spec_completeness_correction,yerr=SF_spec_completeness_correction_err,c='b', marker='x',lolims=False, uplims=False, linewidth=0.0, elinewidth=2, mfc='none',ms=20)
        ax.errorbar(Q_midbins,Q_spec_completeness_correction,yerr=Q_spec_completeness_correction_err,c='r', marker='x',lolims=False, uplims=False, linewidth=0.0, elinewidth=2, mfc='none',ms=20)
        ax.legend(loc='upper right', frameon=False,fontsize=25)
        ax.set_xlabel('$log(M/M_{\odot})$',fontsize=30)
        ax.set_ylim(-0.5,4.1)
        ax.set_ylabel('Correction factor C$_{s}$',fontsize=30)
        ax.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on', labelleft='on',labelsize=25)
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
# ## Fields
SF_field_corrected_smf_HFF = SF_field_raw_smf_HFF + SF_field_mass_completeness_diff
Q_field_corrected_smf_HFF = Q_field_raw_smf_HFF + Q_field_mass_completeness_diff
#
#
## The below calculation has been broken out into a separate function "normalize_smf_mass()" in SECTION (0). We maintain the below work for diagnostic information and summary tables, but we forgo this information for the "cluster field sample" and "parallel field sample", since we're only interested in the result obtained: a normalized SMF.
### NOTE: I tested the function versus the hard-coded algorithm explicitly. They are completely equivalent to within 11 decimal places, more than enough accuracy for the purposes of this calculation
#
## fix the shape of these arrays to be [6 clusters ,# of midbins in SMF], i.e. an array with 6 rows, each row containing an array with (#of midbins) data points
SF_corrected_smf = SF_corrected_smf.reshape((6,len(SF_midbins)))
Q_corrected_smf = Q_corrected_smf.reshape((6,len(SF_midbins)))
SF_field_corrected_smf_HFF = SF_field_corrected_smf_HFF.reshape((6,len(SF_midbins)))
Q_field_corrected_smf_HFF = Q_field_corrected_smf_HFF.reshape((6,len(SF_midbins)))
# SF_field_par_corrected_smf = SF_field_par_corrected_smf.reshape((6,len(SF_midbins)))
# Q_field_par_corrected_smf = Q_field_par_corrected_smf.reshape((6,len(SF_midbins)))
total_corrected_smf = SF_corrected_smf + Q_corrected_smf
total_field_corrected_smf_HFF = SF_field_corrected_smf_HFF + Q_field_corrected_smf_HFF
#
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
    #SF_field_smf,Q_field_smf,total_field_smf = normalize_smf_mass(SF_field_corrected_smf,Q_field_corrected_smf,SF_midbins)
    #SF_field_par_smf,Q_field_par_smf,total_field_par_smf = normalize_smf_mass(SF_field_par_corrected_smf,Q_field_par_corrected_smf,SF_midbins)
    ###
    #
    #
elif normalization_flag == 2:       # normalize by VOLUME
    #
    pass        # USED ONLY FOR FIELD SAMPLE - see program "field_normalization.py"
    #
elif normalization_flag == 3:       # normalize by NUMBER COUNT: UPDATE: THIS IS HOW WE DO IT! flag should be set always to ==3
        #
    #
    ###  CALL THE NORMALIZATION FUNCTION HERE
    # cluster
    SF_smf,Q_smf,total_smf,count_fraction_by_cluster = normalize_smf_count(SF_corrected_smf,Q_corrected_smf,SF_midbins)           # fields
    #SF_field_smf,Q_field_smf,total_field_smf,count_fraction_by_cluster_field = normalize_smf_count(SF_field_corrected_smf,Q_field_corrected_smf,SF_midbins)
    #SF_field_par_smf,Q_field_par_smf,total_field_par_smf,count_fraction_by_cluster_par = normalize_smf_count(SF_field_par_corrected_smf,Q_field_par_corrected_smf,SF_midbins)
    #
    #
#
#
#
#
#
## SECTION (3.4): NORMALIZE the FIELD SMF lists for each cluster by the TOTAL VOLUME subtended by each survey (for field galaxies only, i.e. not the volume occupied by the cluster; see "field_normalization.py" for calculation; RECALL: the UVC field sample is contained in lists "SF_field_uvc_list" & "Q_field_uvc_list", taken from "UVC_master_data.py";
#
#
## combine into a single SMF
SF_field_smf_HFF = ( np.sum(SF_field_corrected_smf_HFF,axis=0) )
Q_field_smf_HFF = ( np.sum(Q_field_corrected_smf_HFF,axis=0) )
#
## make UVC lists into histograms
SF_field_smf_UVC, mass_bins = np.histogram(SF_field_uvc_list, bins=num_bins,range=range2)
Q_field_smf_UVC, mass_bins = np.histogram(Q_field_uvc_list, bins=num_bins,range=range2)
#
## compute relative errors (Poissonian)
#
SF_field_smf_relerr_HFF = (np.sqrt(np.sum(SF_field_raw_smf_HFF,axis=0)) / np.sum(SF_field_raw_smf_HFF,axis=0))
SF_field_smf_relerr_HFF = SF_field_smf_relerr_HFF.reshape(len(SF_midbins))
Q_field_smf_relerr_HFF = (np.sqrt(np.sum(Q_field_raw_smf_HFF,axis=0)) / np.sum(Q_field_raw_smf_HFF,axis=0))
Q_field_smf_relerr_HFF = Q_field_smf_relerr_HFF.reshape(len(SF_midbins))
SF_field_smf_relerr_UVC = (np.sqrt(SF_field_smf_UVC) / SF_field_smf_UVC)
Q_field_smf_relerr_UVC = (np.sqrt(Q_field_smf_UVC) / Q_field_smf_UVC)
#
## run "field_normalization.py" to calculate volume subtended by each survey. volumes are stored in variables "vol_HFF" & "vol_UVC"
#
exec(open('field_normalization.py').read())
#
## Normalize HFF volume to 1 Mpc^3
#
SF_field_smf_HFF = SF_field_smf_HFF / vol_HFF
Q_field_smf_HFF = Q_field_smf_HFF / vol_HFF
total_field_smf_HFF = SF_field_smf_HFF + Q_field_smf_HFF
#
#
## Normalize UVC volume to 1 Mpc^3
#
SF_field_smf_UVC = SF_field_smf_UVC / vol_UVC
Q_field_smf_UVC = Q_field_smf_UVC / vol_UVC
total_field_smf_UVC = SF_field_smf_UVC + Q_field_smf_UVC
#
### NECESSARY DIVERSION: plot both the HFF and UVC field SMFs to check their shape and normalization
#
empty_error = np.zeros_like(SF_field_smf_HFF)
SF_field_err_temp_HFF = SF_field_smf_relerr_HFF * SF_field_smf_HFF
SF_field_err_temp_HFF = SF_field_err_temp_HFF.reshape(len(SF_midbins))
SF_field_err_temp_UVC = SF_field_smf_relerr_UVC * SF_field_smf_UVC
SF_field_err_temp_UVC = SF_field_err_temp_UVC.reshape(len(SF_midbins))
Q_field_err_temp_HFF = Q_field_smf_relerr_HFF * Q_field_smf_HFF
Q_field_err_temp_HFF = Q_field_err_temp_HFF.reshape(len(SF_midbins))
Q_field_err_temp_UVC = Q_field_smf_relerr_UVC * Q_field_smf_UVC
Q_field_err_temp_UVC = Q_field_err_temp_UVC.reshape(len(SF_midbins))
#
#
if field_construction_flag == 1:
    fig = plt.figure()
    string = 'HFF field SMF: %s'%np.round(z_field_bounds[0],decimals=2),' < z < %s'%np.round(z_field_bounds[1],decimals=2)
    fig.suptitle(string, fontsize=30)
    ax = fig.add_subplot(1, 1, 1)
    ax.errorbar(SF_midbins,SF_field_smf_HFF,yerr=SF_field_err_temp_HFF, fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0, label='Star-forming', ms=15)#yerr=SF_error,
    ax.errorbar(Q_midbins,Q_field_smf_HFF,yerr=Q_field_err_temp_HFF,fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0,label='Quiescent', ms=15)
    ax.errorbar(SF_midbins,total_field_smf_HFF,yerr=empty_error,fmt='.k',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0,label='Total', ms=15)
    ax.set_xlabel('$log(M/M_{\odot})$',fontsize=25)
    ax.set_xscale('linear')
    ax.minorticks_on()
    ax.set_xlim(7,12.5)
    ax.set_yscale('log')
    ax.set_ylim(1e-6,0.8)
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both',direction='in',color='k',top=True,left=True, right=True,labelleft=True,labelright=False,labelbottom=True,labelsize=18)
    ax.yaxis.set_label_position("left")
    ax.set_ylabel('???',fontsize=20)
    ax.legend(scatterpoints=1,loc='lower left', frameon=False, fontsize = 'x-small')
    ax.grid(b=False)#, which='major', axis='both', color = 'k', linestyle = '--')
    #
    #
    fig = plt.figure()
    string = 'UVC field SMF: %.2f'%z_field_bounds[0],' < z < %.2f'%z_field_bounds[1],' at UVC cutoff: %.2f'%limiting_mass_uvc
    fig.suptitle(string, fontsize=30)
    ax = fig.add_subplot(1, 1, 1)
    ax.errorbar(SF_midbins,SF_field_smf_UVC,yerr=SF_field_err_temp_UVC, fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0, label='Star-forming', ms=15)#yerr=SF_error,
    ax.errorbar(Q_midbins,Q_field_smf_UVC,yerr=Q_field_err_temp_UVC,fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0,label='Quiescent', ms=15)
    ax.errorbar(SF_midbins,total_field_smf_UVC,yerr=empty_error,fmt='.k',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0,label='Total', ms=15)
    ax.set_xlabel('$log(M/M_{\odot})$',fontsize=25)
    ax.set_xscale('linear')
    ax.minorticks_on()
    ax.set_xlim(7,12.5)
    ax.set_yscale('log')
    ax.set_ylim(1e-6,0.8)
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both',direction='in',color='k',top=True,left=True, right=True,labelleft=True,labelright=False,labelbottom=True,labelsize=18)
    ax.yaxis.set_label_position("left")
    ax.set_ylabel('???',fontsize=20)
    ax.legend(scatterpoints=1,loc='lower left', frameon=False, fontsize = 'x-small')
    ax.grid(b=False)#, which='major', axis='both', color = 'k', linestyle = '--')
#
#
#
#
some_flag = 1           # 1== avg over 3 points near UVC cutoff; 2 == hard transition from HFF curve to UVC curve
if some_flag == 1:
    ## COMBINE HFF & UVC into a single SMF: start by finding the index corresponding to the first midbin with mass >10^9. Then cut the two lists above/below that bin index, then add the two together
    #
    ## search for the first mass bin above 10^9 (i.e. "limiting_mass_uvc"), and return its index
    index = 0
    for ii in range(len(SF_midbins)):
        if SF_midbins[ii] < (limiting_mass_uvc + (bin_width/2)):  # nominally, pts <9 are HFF, >10.2 of UVC, and the two in between are 50/50 (for range2 = [8,12.4])
            index = ii
    print('Index for HFF/UVC cutoff is: %s'%(index+1)+'out of: %s'%len(SF_midbins)+' bins for UVC mass cutoff: %s'%limiting_mass_uvc)
    #
    ## set all entries in HFF >index to zero, and all entries <index to zero for UVC
    for ii in range(len(SF_field_smf_HFF)):
        if ii > (index+2):
            SF_field_smf_HFF[ii] = 0
            Q_field_smf_HFF[ii] = 0
            SF_field_smf_relerr_HFF[ii] = 0
            Q_field_smf_relerr_HFF[ii] = 0
        elif ii == (index+2) or ii == (index+1):
            SF_field_smf_HFF[ii] = 0.5*SF_field_smf_HFF[ii]
            Q_field_smf_HFF[ii] = 0.5*Q_field_smf_HFF[ii]
            SF_field_smf_relerr_HFF[ii] = 0.5*SF_field_smf_relerr_HFF[ii]
            Q_field_smf_relerr_HFF[ii] = 0.5*Q_field_smf_relerr_HFF[ii]
            SF_field_smf_UVC[ii] = 0.5*SF_field_smf_UVC[ii]
            Q_field_smf_UVC[ii] = 0.5*Q_field_smf_UVC[ii]
            SF_field_smf_relerr_UVC[ii] = 0.5*SF_field_smf_relerr_UVC[ii]
            Q_field_smf_relerr_UVC[ii] = 0.5*Q_field_smf_relerr_UVC[ii]
        # if ii > (index+1):
        #     SF_field_smf_HFF[ii] = 0
        #     Q_field_smf_HFF[ii] = 0
        #     SF_field_smf_relerr_HFF[ii] = 0
        #     Q_field_smf_relerr_HFF[ii] = 0
        # elif ii == (index+1):
        #     SF_field_smf_HFF[ii] = 0.5*SF_field_smf_HFF[ii]
        #     Q_field_smf_HFF[ii] = 0.5*Q_field_smf_HFF[ii]
        #     SF_field_smf_relerr_HFF[ii] = 0.5*SF_field_smf_relerr_HFF[ii]
        #     Q_field_smf_relerr_HFF[ii] = 0.5*Q_field_smf_relerr_HFF[ii]
        #     SF_field_smf_UVC[ii] = 0.5*SF_field_smf_UVC[ii]
        #     Q_field_smf_UVC[ii] = 0.5*Q_field_smf_UVC[ii]
        #     SF_field_smf_relerr_UVC[ii] = 0.5*SF_field_smf_relerr_UVC[ii]
        #     Q_field_smf_relerr_UVC[ii] = 0.5*Q_field_smf_relerr_UVC[ii]
        elif ii <= (index+0):
                SF_field_smf_UVC[ii] = 0
                Q_field_smf_UVC[ii] = 0
                SF_field_smf_relerr_UVC[ii] = 0
                Q_field_smf_relerr_UVC[ii] = 0
    #
elif some_flag == 2:
    ## COMBINE HFF & UVC into a single SMF: start by finding the index corresponding to the first midbin with mass >10^9. Then cut the two lists above/below that bin index, then add the two together
    #
    ## search for the first mass bin above 10^9 (i.e. "limiting_mass_uvc"), and return its index
    index = 0
    for ii in range(len(SF_midbins)):
        if SF_midbins[ii] < limiting_mass_uvc:
            index = ii
    print('Index for HFF/UVC cutoff is: %s'%(index+1)+'out of: %s'%len(SF_midbins)+' bins for UVC mass cutoff: %s'%limiting_mass_uvc)
    #
    ## set all entries in HFF >index to zero, and all entries <index to zero for UVC
    for ii in range(len(SF_field_smf_HFF)):
        if ii > index:
            SF_field_smf_HFF[ii] = 0
            Q_field_smf_HFF[ii] = 0
            SF_field_smf_relerr_HFF[ii] = 0
            Q_field_smf_relerr_HFF[ii] = 0
        elif ii < index:
            SF_field_smf_UVC[ii] = 0
            Q_field_smf_UVC[ii] = 0
            SF_field_smf_relerr_UVC[ii] = 0
            Q_field_smf_relerr_UVC[ii] = 0
    #
    ## Now add the two together
#
## check: all errors should be finite. replace NaN entries with zeros
for ii in range(len(SF_field_smf_relerr_UVC)):
    if np.isnan(SF_field_smf_relerr_HFF[ii]) == 1:
        SF_field_smf_relerr_HFF[ii] = 0
    if np.isnan(SF_field_smf_relerr_UVC[ii]) == 1:
        SF_field_smf_relerr_UVC[ii] = 0
    if np.isnan(Q_field_smf_relerr_HFF[ii]) == 1:
        Q_field_smf_relerr_HFF[ii] = 0
    if np.isnan(Q_field_smf_relerr_UVC[ii]) == 1:
        Q_field_smf_relerr_UVC[ii] = 0
#
SF_field_smf = (SF_field_smf_HFF + SF_field_smf_UVC)
Q_field_smf = (Q_field_smf_HFF + Q_field_smf_UVC)
total_field_smf = SF_field_smf + Q_field_smf
#
## compute the actual error bars
SF_field_relerr = SF_field_smf_relerr_HFF + SF_field_smf_relerr_UVC
Q_field_relerr = Q_field_smf_relerr_HFF + Q_field_smf_relerr_UVC
#
SF_field_error = SF_field_relerr * SF_field_smf
Q_field_error = Q_field_relerr * Q_field_smf
total_field_error = SF_field_error + Q_field_error
#
#
#
#
if field_construction_flag == 1:
    fig = plt.figure()
    string = 'COMBINED field SMF: %.3f'%z_field_bounds[0],' < z < %.3f'%z_field_bounds[1],' at UVC cutoff: %.2f'%limiting_mass_uvc
    fig.suptitle(string, fontsize=30)
    ax = fig.add_subplot(1, 1, 1)
    ax.errorbar(SF_midbins,SF_field_smf,yerr=SF_field_error, fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0, label='Star-forming', ms=15)#yerr=SF_error,
    ax.errorbar(Q_midbins,Q_field_smf,yerr=Q_field_error,fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0,label='Quiescent', ms=15)
    ax.errorbar(SF_midbins,total_field_smf,yerr=total_field_error,fmt='.k',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0,label='Total', ms=15)
    ax.set_xlabel('$log(M/M_{\odot})$',fontsize=25)
    ax.set_xscale('linear')
    ax.minorticks_on()
    ax.set_xlim(7,12.5)
    ax.set_yscale('log')
    ax.set_ylim(1e-6,0.8)
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both',direction='in',color='k',top=True,left=True, right=True,labelleft=True,labelright=False,labelbottom=True,labelsize=18)
    ax.yaxis.set_label_position("left")
    ax.set_ylabel('???',fontsize=20)
    ax.legend(scatterpoints=1,loc='lower left', frameon=False, fontsize = 'x-small')
    ax.grid(b=False)#, which='major', axis='both', color = 'k', linestyle = '--')
#
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
## Method 2: treat each bin as its own Poisson distribution, and draw samples (i.e. Bootstrap re-sampling with replacement), where the number of galaxies drawn in each realization of the resampling follows a Poisson distribution, with the mean of the distribution equal to the number count in the SMF mass bin;     DEPRECATED;
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
    SF_spec_completeness_correction = SF_spec_completeness_correction.reshape(len(SF_midbins))
    Q_spec_completeness_correction = Q_spec_completeness_correction.reshape(len(SF_midbins))
    SF_phot_error_rel = np.sqrt(((np.transpose(SF_spec_completeness_correction_err)/SF_spec_completeness_correction)**2) + (SF_phot_error_counting_rel)**2)
    Q_phot_error_rel = np.sqrt(((np.transpose(Q_spec_completeness_correction_err)/Q_spec_completeness_correction)**2) + (Q_phot_error_counting_rel)**2)
    #
    ## compute the weight factor for spec/phot; weight_spec = (# of spec members) / (total # of raw members); weight_phot = 1 - weight_spec
    error_weight_SF = np.array([[0.0]*len(SF_midbins)]*2)  # row1=spec;  row2=phot
    error_weight_Q = np.array([[0.0]*len(SF_midbins)]*2)
    #
    error_weight_SF[0] = np.sum(SF_spec_smf,axis=0) / np.sum(SF_raw_smf,axis=0)
    error_weight_Q[0] = np.sum(SF_spec_smf,axis=0) / np.sum(Q_raw_smf,axis=0)
    #
    ## remove NaNs due to empty bins
    for ii in range(len(error_weight_SF[0])):
        if np.isnan(error_weight_SF[0][ii])==1:
            error_weight_SF[0][ii] = 0
        if np.isnan(error_weight_Q[0][ii])==1:
            error_weight_Q[0][ii] = 0
    #
    error_weight_SF[1] = 1 - error_weight_SF[0]
    error_weight_Q[1] = 1 - error_weight_Q[0]
    #
    ## reshape arrays to match SMF
    SF_phot_error_rel = SF_phot_error_rel.reshape(len(SF_midbins))
    Q_phot_error_rel = Q_phot_error_rel.reshape(len(SF_midbins))
    #
    ## now put it all together: total relative error (squared) = sum of relative errors (squared)
    SF_error_rel = np.sqrt( error_weight_SF[0]*SF_spec_error_rel**2 + error_weight_SF[1]*SF_phot_error_rel**2 )
    Q_error_rel = np.sqrt( error_weight_Q[0]*Q_spec_error_rel**2 + error_weight_Q[1]*Q_phot_error_rel**2 )
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
elif smf_error_method == 3:                  # method 3: treat entire SMF as one sample for counting errors; DEPRECATED - smf_error_method should be set to ==1
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
## SECTION (4.1) compute QUENCHED FRACTION
## defined as (# of quenched galaxies in a given mass bin) / (total # of galaxies in that mass bin)
#
quenched_fraction = np.array([[0]*len(SF_midbins)]*2,dtype='float32')    # row 1 = cluster;   row 2 = field
quenched_fraction_err = np.array([[0]*len(SF_midbins)]*2,dtype='float32')
#
quenched_fraction[0] = Q_smf / total_smf                             # quenched fraction for cluster
quenched_fraction[1] = Q_field_smf / total_field_smf                 # quenched fraction for cluster
#
#
## compute errors on the quenched fractions
#
quenched_fraction_err[0] = np.sqrt( (Q_error / Q_smf)**2 + (total_error / total_smf)**2 ) * quenched_fraction[0]
quenched_fraction_err[1] = np.sqrt( Q_field_relerr**2 + (total_field_error / total_field_smf)**2 ) * quenched_fraction[1]
#
#
#
#### ADDITIONAL STEP ADDED AT LAST MINUTE: need to normalize by bin width, so that axes of SMF are per dex
# # cluster
SF_smf = SF_smf / bin_width
Q_smf = Q_smf / bin_width
total_smf = total_smf / bin_width
#
## compute error bars for SMF plot
SF_error = SF_error_rel * SF_smf
Q_error = Q_error_rel * Q_smf
#
total_error = SF_error + Q_error
#
# #field
SF_field_smf = SF_field_smf / bin_width
Q_field_smf = Q_field_smf / bin_width
total_field_smf = total_field_smf / bin_width
#
## compute error bars
SF_field_error = SF_field_relerr * SF_field_smf
Q_field_error = Q_field_relerr * Q_field_smf
total_field_error = SF_field_error + Q_field_error
#
#
#
### ADDED LAST MINUTE: remove the lowest-mass bin from the field samples
#
if mcmc_flag == 1:
    if range2_flag == 0 and lim_mass_offset_flag == 0:           # if not using offset, delete 1st data point
        SF_field_midbins = np.delete(SF_midbins,0)
        Q_field_midbins = np.delete(Q_midbins,0)
        SF_midbins = np.delete(SF_midbins,0)
        Q_midbins = np.delete(Q_midbins,0)
        SF_smf = np.delete(SF_smf,0)
        SF_error = np.delete(SF_error,0)
        Q_smf = np.delete(Q_smf,0)
        Q_error = np.delete(Q_error,0)
        total_smf = np.delete(total_smf,0)
        total_error = np.delete(total_error,0)
        quenched_fraction = np.delete(quenched_fraction,0,axis=1)
        quenched_fraction_err = np.delete(quenched_fraction_err,0,axis=1)
        SF_field_smf = np.delete(SF_field_smf,0)
        Q_field_smf = np.delete(Q_field_smf,0)
        total_field_smf = np.delete(total_field_smf,0)
        SF_field_error = np.delete(SF_field_error,0)
        Q_field_error = np.delete(Q_field_error,0)
        total_field_error = np.delete(total_field_error,0)
    elif range2_flag == 0 and lim_mass_offset_flag == 1:           # if not using offset, delete 1st data point
        SF_field_midbins = np.delete(SF_midbins,[0,1])
        Q_field_midbins = np.delete(Q_midbins,[0,1])
        SF_midbins = np.delete(SF_midbins,[0,1])
        Q_midbins = np.delete(Q_midbins,[0,1])
        SF_smf = np.delete(SF_smf,[0,1])
        SF_error = np.delete(SF_error,[0,1])
        Q_smf = np.delete(Q_smf,[0,1])
        Q_error = np.delete(Q_error,[0,1])
        total_smf = np.delete(total_smf,[0,1])
        total_error = np.delete(total_error,[0,1])
        quenched_fraction = np.delete(quenched_fraction,[0,1],axis=1)
        quenched_fraction_err = np.delete(quenched_fraction_err,[0,1],axis=1)
        SF_field_smf = np.delete(SF_field_smf,[0,1])
        Q_field_smf = np.delete(Q_field_smf,[0,1])
        total_field_smf = np.delete(total_field_smf,[0,1])
        SF_field_error = np.delete(SF_field_error,[0,1])
        Q_field_error = np.delete(Q_field_error,[0,1])
        total_field_error = np.delete(total_field_error,[0,1])
    elif range2_flag == 1 or range2_flag == 2:
        SF_field_midbins = SF_midbins
        Q_field_midbins = Q_midbins
#
elif mcmc_flag == 0:
    SF_field_midbins = SF_midbins
    Q_field_midbins = Q_midbins
#
## SECTION (5)    EMCEE simulation; see emcee_chi2_final.py;
#
######## This section has been broken out into its own program, called "emcee_chi2", which uses chi-squared as the cost function and is the code to be implemented in the final run for this project
#
#
if mcmc_flag == 1:
    #
    #
    #
    ## Call "spec_membership_selection.py" to determine the preliminary spectroscopic sample
    #
    #
    if mcmc_field_flag == 0:
        exec(open('emcee_chi2_final.py').read())        # fits a single schechter function, for CLUSTER
        #
        #
        print('\n"master_smfz*.py" Section 5 MCMC complete for binning method %i .\n\nPROGRAM SHOULD EXIT AFTER PRINTING THIS STATEMENT'%membership_correction_binning_flag)
        sys.exit()
        print('PROGRAM SHOULD HAVE EXITED BEFORE PRINTING THIS STATEMENT')
        #
    elif mcmc_field_flag == 1:
        if field_smf_schechter_flag == 1:
            exec(open('emcee_chi2_field_single.py').read())        # fits a single schechter function, for CLUSTER
        exec(open('emcee_chi2_double.py').read())        # fits a double schechter function, for FIELD
        #
        #
        print('\n"master_smfz*.py" Section 5 MCMC complete for binning method %i .\n\nPROGRAM SHOULD EXIT AFTER PRINTING THIS STATEMENT'%membership_correction_binning_flag)
        sys.exit()
        print('PROGRAM SHOULD HAVE EXITED BEFORE PRINTING THIS STATEMENT')
    #
    #
    #
    #
#
#
else:
## The following summarizes the result of the MCMC simulation and sets up the appropriate arrays for plotting
    #
    #
    if range2_flag == 0 and z_cutoff_field[1] == 0.10:
        #
        ## pass
        if lim_mass_offset_flag == 0 and cluster_field_inclusion_flag == 0:            # par only
            #
            # ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.07932761e+01,2.67133425e-03,-1.52211218e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [11.02475088,0.0548576,-1.11912586]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.08698214,0.04611485,-1.20737296]
            #
            ## Field
            # #double-schechter Q pop
            SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.06250408e+01,1.34868286e-03 ,-1.15654789e+00,1.33074635e-05,-2.09371535e+00]     # 100 walkers, 25,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08911525e+01,1.57509772e-03,-5.33607622e-01,3.53911493e-06,-2.03959241e+00]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09842577e+01,1.68980938e-03,-1.08541784e+00,1.05363224e-05,-2.04646381e+00]
            #
        elif lim_mass_offset_flag == 0 and cluster_field_inclusion_flag == 1:      # clu + par
            #
            # ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.07932761e+01,2.67133425e-03,-1.52211218e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [11.02475088,0.0548576,-1.11912586]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.08698214,0.04611485,-1.20737296]
            #
            ## Field
            # #double-schechter Q pop
            SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.06267271e+01,1.30963180e-03,-1.11303459e+00,5.34621184e-05,-2.02012805e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08423249e+01,1.67748950e-03,-3.85220256e-01,2.18232373e-05,-1.94929595e+00]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09492543e+01,1.81649458e-03,-9.70100384e-01,6.33980357e-05,-1.94398905e+00]
            #
        elif lim_mass_offset_flag == 1 and cluster_field_inclusion_flag == 0:      # par
            #
            # ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.08210401e+01,2.51004581e-03,-1.53053501e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [10.97345436,0.06593025,-1.07729047]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.06057385,0.05202444,-1.18526766]
            #
            ## Field
            #double-schechter Q pop
            SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.06165938e+01,1.34523402e-03 ,-1.11780109e+00,2.31932949e-05,-2.02703239e+00]     # 100 walkers, 25,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08645556e+01,1.60700496e-03,-4.94498499e-01,4.20169577e-06,-2.03800999e+00]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09523912e+01,1.78624096e-03,-1.03944747e+00,1.77838455e-05,-1.99695055e+00]
            #
        elif lim_mass_offset_flag == 1 and cluster_field_inclusion_flag == 1:      # clu + par
            #
            # pass
            # ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.08210401e+01,2.51004581e-03,-1.53053501e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [10.97345436,0.06593025,-1.07729047]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.06057385,0.05202444,-1.18526766]
            ##
            # ## Field
            # # #double-schechter Q pop
            SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.05951424e+01,1.26243947e-03,-9.49451317e-01,1.92649294e-04,-1.80047401e+00]     # 100 walkers, 15,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08201445e+01,1.69449383e-03,-3.56489254e-01,2.20711044e-05,-1.96054438e+00]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09123466e+01,1.89791151e-03,-8.89690837e-01,1.11715509e-04,-1.85419976e+00]
            #
    elif range2_flag == 0 and z_cutoff_field[1] == 0.125:
        #
        ## pass
        if lim_mass_offset_flag == 0 and cluster_field_inclusion_flag == 0:            # par only
            #
            # pass
            # # ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.07932761e+01,2.67133425e-03,-1.52211218e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [11.02475088,0.0548576,-1.11912586]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.08698214,0.04611485,-1.20737296]
            # #
            # ## Field
            # # #double-schechter Q pop
            SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.06036239e+01,1.37317063e-03, -1.08517157e+00,6.28338502e-05,-1.74087920e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08905300e+01,1.57370132e-03, -5.33208764e-01,3.39966486e-06,-1.96669094e+00]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09749600e+01,1.72709079e-03,-1.07322214e+00,1.02367573e-05,-1.96273061e+00]
            # #
        elif lim_mass_offset_flag == 0 and cluster_field_inclusion_flag == 1:      # clu + par
            #
            # ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.07932761e+01,2.67133425e-03,-1.52211218e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [11.02475088,0.0548576,-1.11912586]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.08698214,0.04611485,-1.20737296]
            #
            ## Field
            # #double-schechter Q pop
            SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.06043947e+01,  1.34618096e-03,-1.05285636e+00,1.03917778e-04,-1.80501764e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08683497e+01,1.62721439e-03,-4.67210415e-01,1.18449523e-05,-1.92368574e+00 ]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09604065e+01,1.77169438e-03,-1.01291735e+00,5.16342526e-05,-1.86079496e+00 ]
            #
        elif lim_mass_offset_flag == 1 and cluster_field_inclusion_flag == 0:      # par only
            #
            # pass
            # ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.08210401e+01,2.51004581e-03,-1.53053501e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [10.97345436,0.06593025,-1.07729047]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.06057385,0.05202444,-1.18526766]
            # #
            # ## Field
            # #double-schechter Q pop
            SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.05994151e+01,1.33094741e-03,-1.06482905e+00,9.56666563e-05,-1.67073543e+00]     # 100 walkers, 15,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08642570e+01,1.61001292e-03,-4.96012787e-01,3.54679674e-06,-1.97554077e+00]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09443440e+01,1.80217089e-03,-1.02401208e+00,2.42873749e-05,-1.85036362e+00]
            # #
        elif lim_mass_offset_flag == 1 and cluster_field_inclusion_flag == 1:      # clu + par
            #
            # pass
            ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.08210401e+01,2.51004581e-03,-1.53053501e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [10.97345436,0.06593025,-1.07729047]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.06057385,0.05202444,-1.18526766]
            ##
            # ## Field
            # # #double-schechter Q pop
            SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.05906712e+01,1.28674740e-03,-9.59819997e-01,1.96579727e-04,-1.69923770e+00]     # 100 walkers, 15,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08482652e+01,1.64494637e-03,-4.45895352e-01,9.74155829e-06,-1.98149026e+00]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09281553e+01,1.85946130e-03,-9.58377257e-01,7.32820586e-05,-1.82061029e+00]
            #
    elif range2_flag == 0 and z_cutoff_field[1] == 0.14:
        #
        ## pass
        if lim_mass_offset_flag == 0 and cluster_field_inclusion_flag == 0:            # par only
            #
            pass                              #BAD VALUES
            # # # ## Cluster
            # SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.07932761e+01,2.67133425e-03,-1.52211218e+00]     # 500 walkers, 10,000 steps
            # QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [11.02475088,0.0548576,-1.11912586]
            # TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.08698214,0.04611485,-1.20737296]
            # # #
            # # ## Field
            # # # #double-schechter Q pop
            # SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.06036239e+01,1.37317063e-03, -1.08517157e+00,6.28338502e-05,-1.74087920e+00]     # 500 walkers, 10,000 steps
            # QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08905300e+01,1.57370132e-03, -5.33208764e-01,3.39966486e-06,-1.96669094e+00]
            # TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09749600e+01,1.72709079e-03,-1.07322214e+00,1.02367573e-05,-1.96273061e+00]
            # #
        elif lim_mass_offset_flag == 0 and cluster_field_inclusion_flag == 1:      # clu + par
            #
            pass
            # ## Cluster                        BAD VALUES
            # #
            # SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.07932761e+01,2.67133425e-03,-1.52211218e+00]     # 500 walkers, 10,000 steps
            # QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [11.02475088,0.0548576,-1.11912586]
            # TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.08698214,0.04611485,-1.20737296]
            # #
            # ## Field
            # # #double-schechter Q pop
            # SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.06043947e+01,  1.34618096e-03,-1.05285636e+00,1.03917778e-04,-1.80501764e+00]     # 500 walkers, 10,000 steps
            # QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08683497e+01,1.62721439e-03,-4.67210415e-01,1.18449523e-05,-1.92368574e+00 ]
            # TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09604065e+01,1.77169438e-03,-1.01291735e+00,5.16342526e-05,-1.86079496e+00 ]
            #
        elif lim_mass_offset_flag == 1 and cluster_field_inclusion_flag == 0:      # par only
            #
            # pass
            # ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.08210401e+01,2.51004581e-03,-1.53053501e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [10.97345436,0.06593025,-1.07729047]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.06057385,0.05202444,-1.18526766]
            # #
            # ## Field
            # #double-schechter Q pop
            SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.06094033e+01,6.28581886e-04,-1.10823989e+00, 1.60904659e-04,-1.46677539e+00]     # 100 walkers, 15,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08655999e+01,8.97294580e-04,-5.02153708e-01,1.29752547e-06,-2.06251664e+00]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09468384e+01,9.88342022e-04,-1.02731089e+00,1.52601179e-05,-1.83983304e+00]
            # #
        elif lim_mass_offset_flag == 1 and cluster_field_inclusion_flag == 1:      # clu + par
            #
            # pass
            ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.08210401e+01,2.51004581e-03,-1.53053501e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [10.97345436,0.06593025,-1.07729047]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.06057385,0.05202444,-1.18526766]
            ##
            # ## Field
            # # #double-schechter Q pop
            SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.05922684e+01,6.84527683e-04,-9.34621112e-01,1.37744160e-04,-1.67381947e+00]     # 100 walkers, 15,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08501226e+01,9.15300471e-04,-4.50894716e-01,5.11694635e-06,-2.00802381e+00]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09302091e+01,1.02501846e-03,-9.55305573e-01, 4.48724043e-05,-1.81946069e+00]
            #
    elif range2_flag == 0 and z_cutoff_field[1] == 0.15:
        #
        ## pass
        if range2_flag == 0 and lim_mass_offset_flag == 0 and cluster_field_inclusion_flag == 0:            # par only
            #
            pass                              #BAD VALUES
            # # # ## Cluster
            # SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.07932761e+01,2.67133425e-03,-1.52211218e+00]     # 500 walkers, 10,000 steps
            # QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [11.02475088,0.0548576,-1.11912586]
            # TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.08698214,0.04611485,-1.20737296]
            # # #
            # # ## Field
            # # # #double-schechter Q pop
            # SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.06036239e+01,1.37317063e-03, -1.08517157e+00,6.28338502e-05,-1.74087920e+00]     # 500 walkers, 10,000 steps
            # QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08905300e+01,1.57370132e-03, -5.33208764e-01,3.39966486e-06,-1.96669094e+00]
            # TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09749600e+01,1.72709079e-03,-1.07322214e+00,1.02367573e-05,-1.96273061e+00]
            # #
        elif range2_flag == 0 and lim_mass_offset_flag == 0 and cluster_field_inclusion_flag == 1:      # clu + par
            #
            pass
            # ## Cluster                        BAD VALUES
            # #
            # SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.07932761e+01,2.67133425e-03,-1.52211218e+00]     # 500 walkers, 10,000 steps
            # QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [11.02475088,0.0548576,-1.11912586]
            # TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.08698214,0.04611485,-1.20737296]
            # #
            # ## Field
            # # #double-schechter Q pop
            # SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.06043947e+01,  1.34618096e-03,-1.05285636e+00,1.03917778e-04,-1.80501764e+00]     # 500 walkers, 10,000 steps
            # QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08683497e+01,1.62721439e-03,-4.67210415e-01,1.18449523e-05,-1.92368574e+00 ]
            # TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09604065e+01,1.77169438e-03,-1.01291735e+00,5.16342526e-05,-1.86079496e+00 ]
            #
        elif range2_flag == 0 and lim_mass_offset_flag == 1 and cluster_field_inclusion_flag == 0:      # par only
            #
            # pass
            # ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.08210401e+01,2.51004581e-03,-1.53053501e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [10.97345436,0.06593025,-1.07729047]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.06057385,0.05202444,-1.18526766]
            # #
            # ## Field
            # #double-schechter Q pop
            SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.06094033e+01,6.28581886e-04,-1.10823989e+00, 1.60904659e-04,-1.46677539e+00]     # 100 walkers, 15,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08655999e+01,8.97294580e-04,-5.02153708e-01,1.29752547e-06,-2.06251664e+00]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09468384e+01,9.88342022e-04,-1.02731089e+00,1.52601179e-05,-1.83983304e+00]
            # #
        elif range2_flag == 0 and lim_mass_offset_flag == 1 and cluster_field_inclusion_flag == 1:      # clu + par
            #
            # pass
            ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.08210401e+01,2.51004581e-03,-1.53053501e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [10.97345436,0.06593025,-1.07729047]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.06057385,0.05202444,-1.18526766]
            ##
            # ## Field
            # # #double-schechter Q pop
            SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.05922684e+01,6.84527683e-04,-9.34621112e-01,1.37744160e-04,-1.67381947e+00]     # 100 walkers, 15,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08501226e+01,9.15300471e-04,-4.50894716e-01,5.11694635e-06,-2.00802381e+00]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09302091e+01,1.02501846e-03,-9.55305573e-01, 4.48724043e-05,-1.81946069e+00]
            #
    elif range2_flag == 1 and z_cutoff_field[1] == 0.125:
        #
        ## pass
        if lim_mass_offset_flag == 1 and cluster_field_inclusion_flag == 0:      # par only
            #
            # pass
            # ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.08210401e+01,2.51004581e-03,-1.53053501e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [10.97345436,0.06593025,-1.07729047]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.06057385,0.05202444,-1.18526766]
            # #
            # ## Field
            # #double-schechter Q pop
            SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.06094033e+01,6.28581886e-04,-1.10823989e+00, 1.60904659e-04,-1.46677539e+00]     # 100 walkers, 15,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08655999e+01,8.97294580e-04,-5.02153708e-01,1.29752547e-06,-2.06251664e+00]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09468384e+01,9.88342022e-04,-1.02731089e+00,1.52601179e-05,-1.83983304e+00]
            # #
        elif lim_mass_offset_flag == 1 and cluster_field_inclusion_flag == 1:      # clu + par
            #
            # pass
            ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.08210401e+01,2.51004581e-03,-1.53053501e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [10.97345436,0.06593025,-1.07729047]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.06057385,0.05202444,-1.18526766]
            ##
            # ## Field
            # # #double-schechter Q pop
            SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.05922684e+01,6.84527683e-04,-9.34621112e-01,1.37744160e-04,-1.67381947e+00]     # 100 walkers, 15,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08501226e+01,9.15300471e-04,-4.50894716e-01,5.11694635e-06,-2.00802381e+00]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09302091e+01,1.02501846e-03,-9.55305573e-01, 4.48724043e-05,-1.81946069e+00]
            #
    elif range2_flag == 1 and z_cutoff_field[1] == 0.14:
    #
    ## pass
        if lim_mass_offset_flag == 1 and cluster_field_inclusion_flag == 0:      # par only
            #
            # pass
            # ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.08210401e+01,2.51004581e-03,-1.53053501e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [10.97345436,0.06593025,-1.07729047]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.06057385,0.05202444,-1.18526766]
            # #
            # ## Field
            # #double-schechter Q pop
            SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.06094033e+01,6.28581886e-04,-1.10823989e+00, 1.60904659e-04,-1.46677539e+00]     # 100 walkers, 15,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08655999e+01,8.97294580e-04,-5.02153708e-01,1.29752547e-06,-2.06251664e+00]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09468384e+01,9.88342022e-04,-1.02731089e+00,1.52601179e-05,-1.83983304e+00]
            # #
        elif lim_mass_offset_flag == 1 and cluster_field_inclusion_flag == 1:      # clu + par
            #
            # pass
            ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.08210401e+01,2.51004581e-03,-1.53053501e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [10.97345436,0.06593025,-1.07729047]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.06057385,0.05202444,-1.18526766]
            ##
            # ## Field
            # # #double-schechter Q pop
            SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.05922684e+01,6.84527683e-04,-9.34621112e-01,1.37744160e-04,-1.67381947e+00]     # 100 walkers, 15,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08501226e+01,9.15300471e-04,-4.50894716e-01,5.11694635e-06,-2.00802381e+00]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09302091e+01,1.02501846e-03,-9.55305573e-01, 4.48724043e-05,-1.81946069e+00]
            #
    elif range2_flag == 1 and z_cutoff_field[1] == 0.15:
    #
    ## pass
        if lim_mass_offset_flag == 1 and cluster_field_inclusion_flag == 1:      # clu + par
            #
            # pass
            ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.10828184e+01,9.65982352e-04,-1.70834634e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [11.07427988,0.05569934,-1.13798742]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.20677233,0.03700197,-1.27802838]
            ##
            # ## Field
            # # #double-schechter Q pop
            SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.06103073e+01,4.22485481e-04,-1.64737127e+00,1.11162273e-03,-8.25053805e-01]     # 100 walkers, 15,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08493729e+01,6.18911676e-06,-2.05586677e+00,1.77212356e-03,-4.49397658e-01]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09203010e+01,1.49657661e-04,-1.77108574e+00,1.99200276e-03,-8.89938215e-01]
            #
    elif range2_flag == 2 and z_cutoff_field[1] == 0.15:    ####  THIS IS THE FINAL SETTING  ####
    #
    ## pass
        if lim_mass_offset_flag == 1 and cluster_field_inclusion_flag == 1:      # clu + par
            #
            # pass
            # ## Cluster
            SFM_star_mcmc, SFphi_mcmc, SFalpha_mcmc = [1.07672040e+01,4.01123411e-03,-1.50639438e+00]     # 500 walkers, 10,000 steps
            QM_star_mcmc, Qphi_mcmc, Qalpha_mcmc = [11.08304502,0.07109424,-1.13096696]
            TM_star_mcmc, Tphi_mcmc, Talpha_mcmc = [11.1388822,0.06171816,-1.20479295]
            # ## Uncertainties
            SFM_star_mcmc_upper_err, SFphi_mcmc_upper_err, SFalpha_mcmc_upper_err = [0,0,0]
            SFM_star_mcmc_lower_err, SFphi_mcmc_lower_err, SFalpha_mcmc_lower_err = [0,0,0]
            QM_star_mcmc_upper_err, Qphi_mcmc_upper_err, Qalpha_mcmc_upper_err = [0,0,0]
            QM_star_mcmc_lower_err, Qphi_mcmc_lower_err, Qalpha_mcmc_lower_err = [0,0,0]
            FM_star_mcmc_upper_err, Tphi_mcmc_upper_err, Talpha_mcmc_upper_err = [0,0,0]
            TM_star_mcmc_lower_err, Tphi_mcmc_lower_err, Talpha_mcmc_lower_err = [0,0,0]
            # #
            # ## Field
            # #double-schechter
            if field_smf_schechter_flag == 1:
                SFM_star_mcmc_field, SFphi_mcmc_field, SFalpha_mcmc_field = [1.08133812e+01,7.64683293e-04,-1.48788260e+00]     # 500 walkers, 10,000 steps
                # UNCERTAINTIES
                SFM_star_field_mcmc_upper_err, SFphi_field_mcmc_upper_err, SFalpha_field_mcmc_upper_err = [1.68426542e-02,3.93628286e-05,1.18640899e-02]
                SFM_star_field_mcmc_lower_err, SFphi_field_mcmc_lower_err, SFalpha_field_mcmc_lower_err = [1.72854603e-02,3.78249826e-05,1.17287409e-02]
            elif field_smf_schechter_flag == 2:
                SFM_star_mcmc_field, SFphi1_mcmc, SFalpha1_mcmc, SFphi2_mcmc, SFalpha2_mcmc = [1.06202533e+01,3.13208121e-04,-1.70734791e+00,1.17692393e-03,-9.34529707e-01]     # 100 walkers, 15,000 steps
            QM_star_mcmc_field, Qphi1_mcmc, Qalpha1_mcmc, Qphi2_mcmc, Qalpha2_mcmc = [1.08933871e+01,2.88605325e-06,-2.19388020e+00,1.64264220e-03,-5.46546062e-01]
            TM_star_mcmc_field, Tphi1_mcmc, Talpha1_mcmc, Tphi2_mcmc, Talpha2_mcmc = [1.09727877e+01,7.43104284e-05,-1.87471712e+00,1.77631629e-03,-1.01697399e+00]
            # UNCERTAINTIES
            QM_star_field_mcmc_upper_err, Qphi1_field_mcmc_upper_err, Qalpha1_field_mcmc_upper_err, Qphi2_field_mcmc_upper_err, Qalpha2_field_mcmc_upper_err = [2.04041047e-02,2.93094916e-06,1.17326104e-01,6.58025860e-05,4.26811050e-02]
            QM_star_field_mcmc_lower_err, Qphi1_field_mcmc_lower_err, Qalpha1_field_mcmc_lower_err, Qphi2_field_mcmc_lower_err, Qalpha2_field_mcmc_lower_err = [2.03420057e-02,1.58991070e-06,1.30903570e-01,6.64712572e-05,4.08266666e-02]
            TM_star_field_mcmc_upper_err, Tphi1_field_mcmc_upper_err, Talpha1_field_mcmc_upper_err, Tphi2_field_mcmc_upper_err, Talpha2_field_mcmc_upper_err = [3.18149949e-02,6.35458958e-05,9.07088767e-02,1.44743656e-04,8.11498153e-02]
            TM_star_field_mcmc_lower_err, Tphi1_field_mcmc_lower_err, Talpha1_field_mcmc_lower_err, Tphi2_field_mcmc_lower_err, Talpha2_field_mcmc_lower_err = [3.23796480e-02,3.88214605e-05,1.06004426e-01,1.38635909e-04,6.91648155e-02]
            # #
    ## define x array to generate points to plot Schechter fit
    x_plot = np.linspace(7.5,12.5,num=1000)#
    x_plot = x_plot.reshape(len(x_plot),)
    #
    ## Cluster models
    SF_model_mcmc_plot = np.log(10) * SFphi_mcmc * (10**((x_plot-SFM_star_mcmc)*(1+SFalpha_mcmc))) * np.exp(-10**(x_plot-SFM_star_mcmc))#
    Q_model_mcmc_plot = np.log(10) * Qphi_mcmc * (10**((x_plot-QM_star_mcmc)*(1+Qalpha_mcmc))) * np.exp(-10**(x_plot-QM_star_mcmc))#
    total_model_mcmc_plot = np.log(10) * Tphi_mcmc * (10**((x_plot-TM_star_mcmc)*(1+Talpha_mcmc))) * np.exp(-10**(x_plot-TM_star_mcmc))#
    ## Field models
    if field_smf_schechter_flag == 1:
        SF_model_field_mcmc_plot = np.log(10) * SFphi_mcmc_field * (10**((x_plot-SFM_star_mcmc_field)*(1+SFalpha_mcmc_field))) * np.exp(-10**(x_plot-SFM_star_mcmc_field))#
    elif field_smf_schechter_flag == 2:
        SF_model_field_mcmc_plot = np.log(10) * math.e**(-10**(x_plot-SFM_star_mcmc_field)) * ( (SFphi1_mcmc*(10**(x_plot-SFM_star_mcmc_field))**(1+SFalpha1_mcmc))  + (SFphi2_mcmc*(10**(x_plot-SFM_star_mcmc_field))**(1+SFalpha2_mcmc)) )
    # Q_model_field_mcmc_plot_single = np.log(10) * Qphi_field_mcmc * (10**((x_plot-QM_field_star_mcmc)*(1+Qalpha_field_mcmc))) * np.exp(-10**(x_plot-QM_field_star_mcmc))#
    Q_model_field_mcmc_plot_double =  np.log(10) * math.e**(-10**(x_plot-QM_star_mcmc_field)) * ( (Qphi1_mcmc*(10**(x_plot-QM_star_mcmc_field))**(1+Qalpha1_mcmc))  + (Qphi2_mcmc*(10**(x_plot-QM_star_mcmc_field))**(1+Qalpha2_mcmc)) )
    # total_model_field_mcmc_plot = np.log(10) * Tphi_field_mcmc * (10**((x_plot-TM_field_star_mcmc)*(1+Talpha_field_mcmc))) * np.exp(-10**(x_plot-TM_field_star_mcmc))#
    total_model_field_mcmc_plot_double = np.log(10) * math.e**(-10**(x_plot-TM_star_mcmc_field)) * ( (Tphi1_mcmc*(10**(x_plot-TM_star_mcmc_field))**(1+Talpha1_mcmc))  + (Tphi2_mcmc*(10**(x_plot-TM_star_mcmc_field))**(1+Talpha2_mcmc)) )
#
#
if UVC_fit_flag == 1:
    #
    SF_midbins_UVC = np.delete(SF_midbins,[0,1,2,3])
    SF_field_smf_UVC = np.delete(SF_field_smf_UVC,[0,1,2,3])
    Q_field_smf_UVC = np.delete(Q_field_smf_UVC,[0,1,2,3])
    #
    ## the following fits the UVC curve alone
    SFM_star_mcmc_UVC, SFphi_mcmc_UVC, SFalpha_mcmc_UVC = [ 1.04188073e+01 , 9.76448133e-04, -8.52526314e-01]
    QM_star_mcmc_UVC, Qphi_mcmc_UVC, Qalpha_mcmc_UVC = [ 1.08724617e+01 , 6.57997600e-04 ,-5.09139371e-01]
    ## define x array to generate points to plot Schechter fit
    x_plot_UVC = np.linspace(9.5,12,num=1000)#
    x_plot_UVC = x_plot.reshape(len(x_plot_UVC),)
    #
    ## UVC models
    SF_model_mcmc_plot_UVC = np.log(10) * SFphi_mcmc_UVC * (10**((x_plot_UVC-SFM_star_mcmc_UVC)*(1+SFalpha_mcmc_UVC))) * np.exp(-10**(x_plot_UVC-SFM_star_mcmc_UVC))#
    Q_model_mcmc_plot_UVC = np.log(10) * Qphi_mcmc_UVC * (10**((x_plot_UVC-QM_star_mcmc_UVC)*(1+Qalpha_mcmc_UVC))) * np.exp(-10**(x_plot_UVC-QM_star_mcmc_UVC))#
    fig = plt.figure()
    string = 'UVC field SMF: %.2f'%z_field_bounds[0],' < z < %.2f'%z_field_bounds[1],' at UVC cutoff: %.2f'%limiting_mass_uvc
    fig.suptitle(string, fontsize=30)
    ax = fig.add_subplot(1, 1, 1)
    ax.errorbar(SF_midbins_UVC,SF_field_smf_UVC,yerr=np.sqrt(SF_field_smf_UVC), fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0, label='Star-forming', ms=15)#yerr=SF_error,
    ax.errorbar(SF_midbins_UVC,Q_field_smf_UVC,yerr=np.sqrt(Q_field_smf_UVC),fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0,label='Quiescent', ms=15)
    ax.errorbar(SF_midbins_UVC,total_field_smf_UVC,yerr=empty_error,fmt='.k',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0,label='Total', ms=15)
    ax.plot(x_plot_UVC,SF_model_mcmc_plot_UVC,'-b')
    ax.plot(x_plot_UVC,Q_model_mcmc_plot_UVC,'-r')
    ax.set_xlabel('$log(M/M_{\odot})$',fontsize=25)
    ax.set_xscale('linear')
    ax.minorticks_on()
    ax.set_xlim(7,12.5)
    ax.set_yscale('log')
    #ax.set_ylim(1e2,2.5e8)
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both',direction='in',color='k',top=True,left=True, right=True,labelleft=True,labelright=False,labelbottom=True,labelsize=18)
    ax.yaxis.set_label_position("left")
    ax.set_ylabel('???',fontsize=20)
    ax.legend(scatterpoints=1,loc='lower left', frameon=False, fontsize = 'x-small')
    ax.grid(b=False)#, which='major', axis='both', color = 'k', linestyle = '--')
#



#
if cluster_only_plot_flag == 1:
    fig = plt.figure()
    ax0 = fig.add_subplot(111)
    ax0.errorbar(SF_midbins,SF_smf,yerr=SF_error, fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0, ms=15)#label='Star-forming', ms=15)#yerr=SF_error,
    ax0.errorbar(Q_midbins,Q_smf,yerr=Q_error,fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0, ms=15)#,label='Quiescent', ms=15)
    ax0.errorbar(SF_midbins,total_smf,yerr=total_error,fmt='.k',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0, ms=15)#,label='Total', ms=15)
    ## Plot Schechter fits:  (uncomment 5 hashtags when fits complete)
    # ax0.plot(x_plot,SF_model_mcmc_plot, '-b', label = 'MCMC - SF', linewidth = 2.0)
    # ax0.plot(x_plot,Q_model_mcmc_plot, '-r', label = 'MCMC - Q', linewidth = 2.0)
    # ax0.plot(x_plot,total_model_mcmc_plot, '-k', label = 'MCMC - Total', linewidth = 2.0)
    ax0.set_xlabel('$log(M/M_{\odot})$',fontsize=25)
    ax0.set_xscale('linear')
    ax0.minorticks_on()
    ax0.set_xlim(7,12.5)
    ax0.set_yscale('log')
    ax0.set_ylim(5e-4,0.7)
    ax0.minorticks_on()
    ax0.tick_params(axis='both', which='both',direction='in',color='k',top=True,left=True, right=True,labelleft=True,labelright=False,labelbottom=True,labelsize=18)
    ax0.yaxis.set_label_position("left")
    ax0.set_ylabel('# density',fontsize=20)
    if lim_mass_offset_flag == 0:
        ax0.set_title('HFF Cluster SMF - w/o offset',fontsize=30)
    elif lim_mass_offset_flag == 1:
        ax0.set_title('HFF Cluster SMF - w/ offset',fontsize=30)
    # SF_string = 'SF: $M^*$: %.3f'%SFM_star_mcmc+'; alpha: %.3f'%SFalpha_mcmc
    # Q_string = 'Q: $M^*$: %.3f'%QM_star_mcmc+'; alpha: %.3f'%Qalpha_mcmc
    # total_string = 'Total: $M^*$: %.3f'%TM_star_mcmc+'; alpha: %.3f'%Talpha_mcmc
    range2_string = 'lim. masses: %s'%limiting_mass
    ax0.text(7.05,4e-2,SF_string,color='b',fontsize='medium')
    ax0.text(7.05,3e-2,Q_string,color='r',fontsize='medium')
    ax0.text(7.05,2e-2,total_string,color='k',fontsize='medium')
    ax0.text(7.05,5e-3,range2_string,color='k',fontsize='medium')
    ax0.legend(scatterpoints=1,loc='lower left', frameon=False, fontsize = 'medium')
    ax0.grid(b=False)#, which='major', axis='both', color = 'k', linestyle = '--')
    #
#
#
if field_only_plot_flag == 1:
    fig = plt.figure()
    ax0 = fig.add_subplot(111)
    ax0.errorbar(SF_field_midbins,SF_field_smf,yerr=SF_field_error, fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0, ms=15)#label='Star-forming', ms=15)#yerr=SF_error,
    ax0.errorbar(Q_field_midbins,Q_field_smf,yerr=Q_field_error,fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0, ms=15)#,label='Quiescent', ms=15)
    ax0.errorbar(SF_field_midbins,total_field_smf,yerr=total_field_error,fmt='.k',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0, ms=15)#,label='Total', ms=15)
    ## Plot Schechter fits:  (uncomment 5 hashtags when fits complete)
    # ax0.plot(x_plot,SF_model_field_mcmc_plot, '-b', label = 'MCMC - SF', linewidth = 2.0)
    # ax0.plot(x_plot,Q_model_field_mcmc_plot_single, '-r', label = 'MCMC - Q', linewidth = 2.0)
    # ax0.plot(x_plot,total_model_field_mcmc_plot, '-k', label = 'MCMC - Total', linewidth = 2.0)
    # ax0.plot(x_plot,SF_model_field_mcmc_plot_double, '-b', label = 'MCMC - Q', linewidth = 2.0)
    # ax0.plot(x_plot,Q_model_field_mcmc_plot_double, '-r', label = 'MCMC - Q', linewidth = 2.0)
    # ax0.plot(x_plot,total_model_field_mcmc_plot_double, '-k', label = 'MCMC - Total', linewidth = 2.0)
    ax0.set_xlabel('$log(M/M_{\odot})$',fontsize=25)
    ax0.set_xscale('linear')
    ax0.minorticks_on()
    ax0.set_xlim(7,12.5)
    ax0.set_yscale('log')
    ax0.set_ylim(1e-6,0.8)
    ax0.minorticks_on()
    ax0.tick_params(axis='both', which='both',direction='in',color='k',top=True,left=True, right=True,labelleft=True,labelright=False,labelbottom=True,labelsize=18)
    ax0.yaxis.set_label_position("left")
    ax0.set_ylabel('# density',fontsize=20)
    if lim_mass_offset_flag == 0:
        ax0.set_title('HFF+UVC FIELD - w/o offset',fontsize=30)
    elif lim_mass_offset_flag == 1:
        ax0.set_title('HFF+UVC FIELD SMF - w/ offset',fontsize=30)
    if cluster_field_inclusion_flag == 0:
        ax0.text(7.05,1e-4,'Field: Just Parallel',fontsize=15)
    elif cluster_field_inclusion_flag == 1:
        ax0.text(7.05,1e-4,'Field: Cluster + Parallel',fontsize=15)
    # SF_string = 'SF: $M^*$: %.3f'%SFM_field_star_mcmc+'; alpha: %.3f'%SFalpha_field_mcmc
    # Q_string = 'Q: $M^*$: %.3f'%QM_field_star_mcmc+'; alpha: %.3f'%Qalpha_field_mcmc
    # total_string = 'Total: $M^*$: %.3f'%TM_field_star_mcmc+'; alpha: %.3f'%Talpha_field_mcmc
    # ax0.text(7.05,4e-4,SF_string,color='b',fontsize='medium')
    # ax0.text(7.05,3e-4,Q_string,color='r',fontsize='medium')
    # ax0.text(7.05,2e-4,total_string,color='k',fontsize='medium')
    ax0.legend(scatterpoints=1,loc='lower left', frameon=False, fontsize = 'medium')
    ax0.grid(b=False)#, which='major', axis='both', color = 'k', linestyle = '--')
    #
#
#
#
if (plot_flag_2 == 1 and project_plot_flag ==2) or project_plot_flag == 1: # plot interpolated/extrapolated points on top of computed correction fractions
    if project_plot_flag == 0:
        pass
    else:
    ## upper: SMF for cluster, field;       lower: fractions of SF/Q in cluster, field
    #
        #plt.close()
        SMF = plt.figure()
        # string = 'Spec: %s'%z_cutoff[0]+'  Phot: %s'%z_cutoff[1]+'  Method: %i'%membership_correction_binning_flag
        if lim_mass_offset_flag == 0:
            string = 'W/O offset'
        elif lim_mass_offset_flag == 1:
            string = 'W/ offset'
        SMF.suptitle(string, fontsize=30)
        gs = gridspec.GridSpec(2,2, wspace=0, hspace=0, width_ratios=[1,1], height_ratios=[2,1])   #make a tiled-plot like vdB2013 w/ fractions below, this line sets the proporitons of plots in the figure
        #gs = gridspec.GridSpec(2,3, width_ratios=[1,1,1], height_ratios=[2,1])   #make a tiled-plot like vdB2013 w/ fractions below, this line sets the proporitons of plots in the figure
        #
        ## CLUSTER
        ax0 = plt.subplot(gs[0])
        ax0.errorbar(SF_midbins,SF_smf,yerr=SF_error, fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0, label='Star-forming', ms=15)#yerr=SF_error,
        ax0.errorbar(Q_midbins,Q_smf,yerr=Q_error,fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0,label='Quiescent', ms=15)
        ax0.errorbar(SF_midbins,total_smf,yerr=total_error,fmt='.k',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0,label='Total', ms=15)
        ## Plot Schechter fits:  (uncomment 5 hashtags when fits complete)
        if schechter_plot_flag == 1:
            ax0.plot(x_plot,SF_model_mcmc_plot, '-b',  linewidth = 2.0)#label = 'MCMC - SF',
            ax0.plot(x_plot,Q_model_mcmc_plot, '-r', linewidth = 2.0)#, label = 'MCMC - Q'
            ax0.plot(x_plot,total_model_mcmc_plot, '-k', linewidth = 2.0)#, label = 'MCMC - Total'
        ax0.set_xlabel('$log($M_{*}$/M_{\odot})$',fontsize=20)
        ax0.set_xscale('linear')
        ax0.minorticks_on()
        ax0.set_xlim(7.8,12.5)
        ax0.set_yscale('log')
        ax0.set_ylim(5e-4,0.8)
        ax0.minorticks_on()
        ax0.tick_params(axis='both', which='both',direction='in',color='k',top=True,left=True, right=True,labelleft=True,labelright=False,labelbottom=False,labelsize=18)
        ax0.yaxis.set_label_position("left")
        ax0.set_ylabel('N density',fontsize=20)
        ax0.set_title('Cluster',fontsize=30)
        if schechter_plot_flag == 1:
            SF_string = 'SF: $M_{*}$: %.3f'%SFM_star_mcmc+'; alpha: %.3f'%SFalpha_mcmc
            Q_string = 'Q: $M_{*}$: %.3f'%QM_star_mcmc+'; alpha: %.3f'%Qalpha_mcmc
            total_string = 'Total: $M_{*}$: %.3f'%TM_star_mcmc+'; alpha: %.3f'%Talpha_mcmc
            ax0.text(7.85,1.8e-2,SF_string,color='b',fontsize='medium')
            ax0.text(7.85,1.4e-2,Q_string,color='r',fontsize='medium')
            ax0.text(7.85,1.1e-2,total_string,color='k',fontsize='medium')
        limmass_string = 'Cluster:\nNames: [M0416,M1149,M0717,A370,A1063,A2744] \nlim. masses: %s'%limiting_mass+'\nz_cluster: %s'%z_cluster+'\ndelta_z_cutoff_clu [spec,phot]: %s'%z_cutoff
        ax0.text(7.85,3e-3,limmass_string,color='k',fontsize='medium')
        range_string = 'SMF range: %s'%np.round(range2,decimals=2)
        ax0.text(7.85,2e-3,range_string,color='k',fontsize='medium')
        ax0.legend(scatterpoints=1,loc='lower left', frameon=False, fontsize = 'medium')
        ax0.grid(b=False)#, which='major', axis='both', color = 'k', linestyle = '--')
        #
        ## cluster fraction
        ax2 = plt.subplot(gs[2])
        #plt.plot(SF_midbins,frac_smf[0],'.b',linewidth=0.5)
        #plt.plot(SF_midbins,frac_smf[1],'.r',linewidth=0.5)
        ax2.errorbar(SF_midbins,quenched_fraction[0],yerr=quenched_fraction_err[0], fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0, ms=15)
        #ax2.errorbar(SF_midbins,frac_smf[1],yerr=frac_error[1], fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
        ax2.set_xscale('linear')
        ax2.set_xlabel('$log($M_{*}$/M_{\odot})$',fontsize=20)
        ax2.set_xlim(7.8,12.5)
        ax2.set_yscale('linear')
        ax2.set_ylim(0,1.1)
        ax2.minorticks_on()
        ax2.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=False,labelsize=18)
        ax2.yaxis.set_label_position("left")
        ax2.set_ylabel('Quenched fraction',fontsize=20)
        #
        ## FIELD
        ax1 = plt.subplot(gs[1])
        ax1.errorbar(SF_field_midbins,SF_field_smf,yerr=SF_field_error, fmt='.b',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0, label='Star-forming', ms=15)#yerr=SF_error,
        ax1.errorbar(Q_field_midbins,Q_field_smf,yerr=Q_field_error,fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0,label='Quiescent', ms=15)
        ax1.errorbar(SF_field_midbins,total_field_smf,yerr=total_field_error,fmt='.k',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0,label='Total', ms=15)
        ## Plot Schechter fits:  (uncomment 5 hashtags when fits complete)
        if schechter_plot_flag == 1:
            ax1.plot(x_plot,SF_model_field_mcmc_plot, '-b')#, label = 'MCMC - Q', linewidth = 2.0)
            ax1.plot(x_plot,Q_model_field_mcmc_plot_double, '-r')#, label = 'MCMC - Q', linewidth = 2.0)
            ax1.plot(x_plot,total_model_field_mcmc_plot_double, '-k')#, label = 'MCMC - Q', linewidth = 2.0)
        range2_string_par = 'lim. masses (par): %s'%limiting_mass_par
        bounds_string_par = 'z_field_bounds: %s'%np.round(z_field_bounds,decimals=3)
        # ax1.text(7.85,5e-5,range2_string_par,color='k',fontsize='medium')
        ax1.text(7.85,3.25e-6,bounds_string_par,color='k',fontsize='medium')
        if cluster_field_inclusion_flag == 0:
            ax1.text(7.85,1.5e-4,'Field: Just Parallel',fontsize='medium')
        elif cluster_field_inclusion_flag == 1:
            range2_string_cluster_field = 'lim. masses (cluster field): %s'%limiting_mass_cluster_field
            ax1.text(7.85,0.8e-4,'Field: \nHFF: Cluster + Parallel ($M_{*}$ <$10^{9.2}$)\n 50/50 - HFF/UltraVISTA ($10^{9.2}$ < $M_{*}$ < $10^{10}$)\nCOSMOS/UltraVISTA ($M_{*}$ >$10^{10}$)',fontsize='medium')
            # ax1.text(7.85,8e-5,range2_string_cluster_field,color='k',fontsize='medium')
        #
        if schechter_plot_flag == 1:
            if field_smf_schechter_flag == 1:
                SF_string = 'SF: $M_{*}$: %.3f'%SFM_star_mcmc_field+'; alpha: %.3f'%SFalpha_mcmc_field
            elif field_smf_schechter_flag == 2:
                SF_string = 'SF: $M_{*}$: %.3f'%SFM_star_mcmc_field+'; alpha1: %.3f'%SFalpha2_mcmc+'; alpha2: %.3f'%SFalpha1_mcmc
            Q_string = 'Q: $M^*$: %.3f'%QM_star_mcmc_field+'; alpha1: %.3f'%Qalpha2_mcmc+'; alpha2: %.3f'%Qalpha1_mcmc
            total_string = 'Total: $M_{*}$: %.3f'%TM_star_mcmc_field+'; alpha1: %.3f'%Talpha2_mcmc+'; alpha2: %.3f'%Talpha1_mcmc
            ax1.text(7.85,1.5e-5,SF_string,color='b',fontsize='medium')
            ax1.text(7.85,9.5e-6,Q_string,color='r',fontsize='medium')
            ax1.text(7.85,6e-6,total_string,color='k',fontsize='medium')
        #
        ax1.text(7.85,2e-6,'delta_z_cutoff_field [phot]: %s'%z_cutoff_field[1],color='k',fontsize='medium')
        ax1.set_xlabel('$log($M_{*}$/M_{\odot})$',fontsize=20)
        ax1.set_xscale('linear')
        ax1.minorticks_on()
        ax1.set_xlim(7.8,12.5)
        ax1.set_yscale('log')
        ax1.set_ylim(1e-6,0.8)
        ax1.minorticks_on()
        ax1.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=True,labelleft=False,labelbottom=False,labelsize=18)
        ax1.yaxis.set_label_position("right")
        ax1.set_ylabel('${\phi}$ [$Mpc^{-3}$ $dex^{-1}$]',fontsize=20)
        ax1.set_title('Field',fontsize=30)
        #ax3.legend(scatterpoints=1,loc='lower left', frameon=False, fontsize = 'x-small')
        ax1.grid(b=False)#, which='major', axis='both', color = 'k', linestyle = '--')
        #
        ## field fraction
        ax3 = plt.subplot(gs[3])
        ax3.errorbar(SF_midbins,quenched_fraction[1],yerr=quenched_fraction_err[1], fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=2.0, ms=15)
        #ax2.errorbar(SF_midbins,frac_smf[1],yerr=frac_error[1], fmt='.r',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
        ax3.set_xscale('linear')
        ax3.set_xlabel('$log($M_{*}$/M_{\odot})$',fontsize=25)
        ax3.set_xlim(7.8,12.5)
        ax3.set_yscale('linear')
        ax3.set_ylim(0,1.1)
        ax3.minorticks_on()
        ax3.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelleft=False,labelright=True,labelsize=18)
        ax3.yaxis.set_label_position("right")
        ax3.set_ylabel('Quenched fraction',fontsize=20)
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
#
#
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
if where_im_at_flag == 1:
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
print('"master_smfz_9_final.py" terminated successfully.')
#
