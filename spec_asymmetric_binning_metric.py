# Created on Fri Jul 03 03:43:46 2020
#
################## spec_asymmetric_binning_metric.py ##################
#
#
## This program used to be Section (5) of "spec_completeness_binning.py", but gets called multiple times so now it's its own script. The file COMPUTES MIDBINS, false pos/neg RATIOS, and the "METRIC OF MERIT" for binning methods 2&3 (asymmetric binning)
#
#
#
#
#
###################     PROGRAM START     ###################
#
#
#
## SF
#
# make histograms
SF_pos_hist, bins_SF = np.histogram(SF_pos, bins=num_bins_to_try[number], range=range2)
SF_neg_hist, bins_SF = np.histogram(SF_neg, bins=num_bins_to_try[number], range=range2)
#
bins_SF = np.round(bins_SF,decimals=3)
#
SF_ratio = np.round((SF_pos_hist/SF_neg_hist),decimals=3)
for jj in range(len(SF_pos_hist)):
    if SF_pos_hist[jj] == 0 and SF_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
        SF_ratio[jj] = float('NaN')
#
## compute midbins for spec. mass completeness plot (i.e. plot of false pos/false neg ratios)
SF_ratio_midbins = midbins(bins_SF)
#
#
if np.sum(np.isnan(SF_ratio)) > 0:     # if the ratios are bad, don't both interpolating the correction factors
    SF_var = float('NaN')
else:                                  # call correction factors file and compute metric of merit
    #
    #
    #    
    ## call a new file: "correction_factors.py" to interpolate/extrapolate the correction factors to the SMF, and return the metric of merit
    #
    exec(open('spec_correction_factors.py').read())      #opens and executes the script
    #
    #
    ## compute variance of SF/Q ratios from 1
    SF_var = SF_metric        
#
#       
## Q
Q_pos_hist, bins_Q = np.histogram(Q_pos, bins=num_bins_to_try[number], range=range2)
Q_neg_hist, bins_Q = np.histogram(Q_neg, bins=num_bins_to_try[number], range=range2)
#
bins_Q = np.round(bins_Q,decimals=3)
#
Q_ratio = np.round((Q_pos_hist/Q_neg_hist),decimals=3)
for jj in range(len(Q_pos_hist)):
    if Q_pos_hist[jj] == 0 and Q_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
        Q_ratio[jj] = 1
    elif Q_pos_hist[jj] == 0 or Q_neg_hist[jj] == 0:
        Q_ratio[jj] = float('NaN')
#
# compute midbins for spec. mass completeness plot (i.e. plot of false pos/false neg ratios)
Q_ratio_midbins = midbins(bins_Q)
#
if np.sum(np.isnan(Q_ratio)) > 0:     # if the ratios are bad, don't both interpolating the correction factors
    Q_var = float('NaN')
else:                                  # call correction factors file and compute metric of merit
    #
    #
    #    
    ## call a new file: "correction_factors.py" to interpolate/extrapolate the correction factors to the SMF, and return the metric of merit
    #
    exec(open('spec_correction_factors.py').read())      #opens and executes the script
    #
    #
    ## compute variance of SF/Q ratios from 1
    Q_var = Q_metric      
#
#
#
#
#
#
#
###################     PROGRAM END     ###################    