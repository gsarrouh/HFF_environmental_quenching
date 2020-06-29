#Created on Sun Jun 28 17:48:55 2020
#
################## correction_factors.py ##################
#
#
#
### This program is based on code from the latter half of "master_smfz*.py" Section (3.2). The purpose of this file is to be called iteratively by "spec_completeness_binning.py" during the diagnostic loop varying bin numbers/techniques for various redshift cuts to cluster membership. This file will look at the spectroscopic completeness CORRECTION FACTORS calculated in "spec_completeness_binning.py", INTERPOLATE/EXTRAPOLATE the correction factors to the SMF mass bins, calculate and return the "metric of merit", which is the average squared deviation from 1 for all correction factors (i.e. {sum_of squares_across_all_bins_of (1 - correction_factor_for_a_given_bin)} divided by the number of bins.
#
#
#
### Section summary:
#
### PROGRAM START
#
### EVALUATE DIFFERENT BIN #s & BINNING TECHNIQUES:
### (0)    ;
### (1)    ; 
### (2)    ;
### (3)    ;
### (4)    compute METRIC OF MERIT;
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
#
#
def metric(correction_factors):
    sum_of_sq_deviations = 0
    if correction_factors[ii] < 0:
        return float('NaN')
    else:
        for ii in range(len(correction_factors)):
            sum_of_sq_deviations = sum_of_sq_deviations + (1 - correction_factors[ii])**2
        metric_of_merit = (np.round(sum_of_sq_deviations[ii] / len(correction_factors),decimals=5))
        return metric_of_merit
#
#
## SECTION (1)
#
## Now interpolate/extrapolate between the SF_ratio/Q_ratio data points from "spec_completeness_binning.py"
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
    if SF_spec_completeness_correction[ii] == 0:     # only write in cell once, don't overwrite cell once correction factor has been computed
        if SF_midbins[ii] < SF_frac_midbins[0]:      # extrapolate below lowest mass bin
            SF_spec_completeness_correction[ii] = m_SF[0]*SF_midbins[ii] + b_SF[0]    
        elif SF_midbins[ii] > SF_frac_midbins[-1]:    # extrapolate above highest mass bin
            SF_spec_completeness_correction[ii] = m_SF[-1]*SF_midbins[ii] + b_SF[-1]    
        elif SF_midbins[ii] > SF_frac_midbins[0] and SF_midbins[ii] < SF_frac_midbins[-1]:    # interpolate in between all other points
            for jj in range(len(SF_frac_midbins)-1):
                if SF_midbins[ii] > SF_frac_midbins[jj] and SF_midbins[ii] < SF_frac_midbins[jj+1]:
                    SF_spec_completeness_correction[ii] = m_SF[jj]*SF_midbins[ii] + b_SF[jj]
        else:
            SF_spec_completeness_correction[ii] = float('NaN')
            print('Error in SF spec completeness correction computation. ABORT')
            break
#
## now compute the metric of merit
SF_metric = metric(SF_spec_completeness_correction)
#
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
## now compute the metric of merit
Q_metric = metric(Q_spec_completeness_correction)
#




