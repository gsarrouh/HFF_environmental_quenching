# Created on Sun Jun 28 17:48:55 2020
#
################## spec_correction_factors.py ##################
#
#
### This program is based on code from the latter half of "master_smfz*.py" Section (3.2). The purpose of this file is to be CALLED ITERATIVELY BY "spec_completeness_binning.py" during the diagnostic loop varying bin numbers/binning methods for various redshift cuts to cluster membership. 
#
### This program will take as INPUTS the spectroscopic FALSE POS/NEG RATIO calculated in "spec_completeness_binning.py", INTERPOLATE/EXTRAPOLATE the CORRECTION FACTORS to the SMF mass bins, calculate and return the "METRIC OF MERIT", which is the average squared deviation from 1 for all correction factors (i.e. {sum_of_squares_across_all_bins_of (1 - correction_factor_for_a_given_bin)} divided by the_number_of_bins.
#
#
#
### Section summary:
#
### PROGRAM START
#
### EVALUATE DIFFERENT BIN #s & BINNING TECHNIQUES:
### (0)    import modules, definitions, flags
### (1)    INTERPOLATE/EXTRAPOLATE SMF midbins onto false pos/neg ratios to determine spec. completeness correction factor by mass bin for the SMF;  
### (2):   compute METRIC of MERIT;
#
### PROGRAM END
#
#
#
#
###################     PROGRAM START     ###################
#
#
## SECTION (0): modules, definitions, flags
#
## MODULES
#
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
#
## DEFINITIONS
#
def metric(correction_factors):
    for ii in range(len(correction_factors)):
        sum_of_sq_deviations = 0
        if correction_factors[ii] < 0 or np.sum(np.isnan(correction_factors)) > 0:
            return float('NaN')
        else:
            sum_of_sq_deviations = sum_of_sq_deviations + (1 - correction_factors[ii])**2
            metric_of_merit = (np.round(sum_of_sq_deviations[ii] / len(correction_factors),decimals=5))
            return metric_of_merit
#
## FLAGS
#
## Diagnostic flag: 0=off, suppress all diagnostic output;  1=on, enable diagnostic output;
diag_flag = 0     
#
## SECTION (1)
#
## Now interpolate/extrapolate between the SF_ratio/Q_ratio data points from "spec_completeness_binning.py"
#
# initialize arrays to store slopes/intercepts for extrapolation/interpolation of spec mass completeness correction factors
m_SF = np.zeros((len(SF_ratio_midbins)-1))     
b_SF = np.zeros((len(SF_ratio_midbins)-1))
m_Q = np.zeros((len(Q_ratio_midbins)-1))     
b_Q = np.zeros((len(Q_ratio_midbins)-1))
SF_spec_completeness_correction = np.zeros_like(SF_midbins,dtype='float32')   # "SF_midbins are the bin midpoints of SF STELLAR MASS FUNCTION
Q_spec_completeness_correction = np.zeros_like(SF_midbins,dtype='float32')     # NOT A MISTAKE: the "Q_midbins" array has been offset from the SF_midbins"" array by 0.05 for plotting/visual clarity, but the bin points themselves are the same, so the correction factors should be computed for "SF_midbins", not "Q_midbins"
#
## SF
for ii in range(len(SF_ratio_midbins)-1):
    m_SF[ii] = (SF_ratio[ii+1] - SF_ratio[ii]) / (SF_ratio_midbins[ii+1] - SF_ratio_midbins[ii]) # calc slope
    b_SF[ii] = SF_ratio[ii] - (SF_ratio_midbins[ii]*m_SF[ii])   # calc intercept
#
for ii in range(len(SF_midbins)):
    if SF_spec_completeness_correction[ii] == 0:     # only write in cell once, don't overwrite cell once correction factor has been computed
        if SF_midbins[ii] < SF_ratio_midbins[0]:      # extrapolate below lowest mass bin
            SF_spec_completeness_correction[ii] = m_SF[0]*SF_midbins[ii] + b_SF[0]    
        elif SF_midbins[ii] > SF_ratio_midbins[-1]:    # extrapolate above highest mass bin
            SF_spec_completeness_correction[ii] = m_SF[-1]*SF_midbins[ii] + b_SF[-1]    
        elif SF_midbins[ii] > SF_ratio_midbins[0] and SF_midbins[ii] < SF_ratio_midbins[-1]:    # interpolate in between all other points
            for jj in range(len(SF_ratio_midbins)-1):
                if SF_midbins[ii] > SF_ratio_midbins[jj] and SF_midbins[ii] < SF_ratio_midbins[jj+1]:
                    SF_spec_completeness_correction[ii] = m_SF[jj]*SF_midbins[ii] + b_SF[jj]
        else:
            SF_spec_completeness_correction[ii] = float('NaN')
            if diag_flag == 1 or project_diagnostic_flag == 1:
                if project_diagnostic_flag == 0:
                    pass
                else:                
                    print('Error in SF spec completeness correction computation. ABORT')
            pass#break
#
#
## Q: NOTE: it is NOT A MISTAKE that the below code uses "SF_midbins" instead of "Q_midbins". see comment above next to "Q_spec_completeness_correction" initialization
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
            Q_spec_completeness_correction[ii] = float('NaN')
            if diag_flag == 1 or project_diagnostic_flag == 1:
                if project_diagnostic_flag == 0:
                    pass
                else:                
                    print('Error in Q spec completeness correction computation. ABORT')
            pass#break   
#
#
## SECTION (2): compute METRIC of MERIT
# 
SF_metric = metric(SF_spec_completeness_correction)
Q_metric = metric(Q_spec_completeness_correction)
#
#
#
#
#
#
#                        
###################     PROGRAM END     ###################



