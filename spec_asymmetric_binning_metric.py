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

### SECTION (5): compute BIN EDGES, CORRECTION FACTORS, & METRIC OF MERIT
## The goal now is to produce a list of bin edges that will be applied to both the false pos/neg lists. for the ii'th bin edge from each of the false pos/neg lists, test all possible values for the ii'th bin between false_pos_bin_edge[ii] & false_neg_bin_edge[ii] in increments of 0.01. choose the bin edge which yields the lowest average deviation from 1. then call (i.e. execute the file) "correction_factors.py", which is based on code from "master_smfz*" S3.2. It interpolates/extrapolates SMF bin correction factors from the false pos/neg ratios computed below, and returns the "metric of merit", i.e. (the sum of squared deviations from 1 for each pos/neg ratio bin) divided by (the number of bins). 
#
## ignore 'divide by zero' errors, since they are handled with exceptions below
np.seterr(divide='ignore')
#
#
## SF
#
## add a loop to ignore the case where the number of bins exceeds the population of false pos/neg.
pop_SF_pos = len(SF_pos)
pop_SF_neg = len(SF_neg)
 #
if num_bins_to_try[number] > pop_SF_pos or num_bins_to_try[number] > pop_SF_neg:
    if diag_flag == 1 or project_diagnostic_flag == 1:
        if project_diagnostic_flag == 0:
            pass
        else:
            print('ERROR: # of bins exceeds population of SF false pos/neg for %s'%num_bins_to_try[number],' bins and cutoffs:\n%s'%z_cutoff[0],' spec and %s'%z_cutoff[1],' phot\nUSE FEWER BINS')
            SF_var = float('NaN')
#
else:
    #
    for ii in range(len(bin_edge_means_SF)):
        #
        if ii != (len(bin_edge_means_SF)-1):                    # if not the last bin edge (i.e. upper limit of mass range for study)...
            if bin_edge_means_SF[ii] > bin_edge_means_SF[(ii+1)]:   # if the next bin edge is lower than the current bin edge (in mass)...
              #  FIX THIS HERE: s/b looking at bin edges in num_bins_*, not bin_edge_means*
                #
                ## enable diagnostic output
                if (diag_flag == 1 and project_diagnostic_flag ==2) or project_diagnostic_flag == 1:
                    if project_diagnostic_flag == 0:
                        pass
                    else:
                        print("ERROR 2: the next bin value is lower than the current bin value for cutoffs - spec: %s"%z_cutoff[0],";  phot: %s"%z_cutoff[1],";   for %s"%num_bins_to_try[number]," bins, method = %s"%method,'\nThe current bin is being set to the value of the next bin.\nUSE FEWER BINS')
                #
                bin_edge_means_SF[ii] = bin_edge_means_SF[(ii+1)] 
            else:
                #
                #
                SF_var_list = []
                if num_bins_SF[0][ii] < num_bins_SF[1][ii]:    # this if statement won't apply to the first/last bin, since they have been set to equal "range2", the mass range for the study
                    if ii == 0:
                        bin_edge_SF_start = num_bins_SF[0][ii]     # store the lower bound of the ii'th bin edge b/w the false pos/neg bin edges
                    else:
                        bin_edge_SF_start = max(num_bins_SF[0][ii],bin_edge_means_SF[(ii-1)])
                    bin_edge_SF_end = num_bins_SF[1][ii]
                elif num_bins_SF[0][ii] > num_bins_SF[1][ii]:
                    if ii == 0:
                        bin_edge_SF_start = num_bins_SF[1][ii]
                    else:
                        bin_edge_SF_start = max(num_bins_SF[1][ii],bin_edge_means_SF[(ii-1)])
                    bin_edge_SF_end = num_bins_SF[0][ii]
                else:                                          # for when the i'th bin edge is the same for both false pos/neg lists
                    equal_flag = 1
                    bin_edge_SF_start = num_bins_SF[0][ii]
                    bin_edge_SF_end = num_bins_SF[0][ii]
                    mid_pt_to_use = num_bins_SF[0][ii]
                #
                if bin_edge_SF_start == bin_edge_SF_end:
                    equal_flag = 1
                else:
                    equal_flag = 0
                #
                ## enable diagnostic output
                if (diag_flag == 1 and project_diagnostic_flag ==2) or project_diagnostic_flag == 1:
                    if project_diagnostic_flag == 0:
                        pass
                    else:
                        print('Start of : %.2f'%bin_edge_SF_start,' for SF method: %s'%method,' for z_cutoff = %s'%z_cutoff)
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
                        SF_var = np.sum((1 - SF_ratio)**2) / len(SF_ratio)      # i.e. divide by num_bins_to_try[number]
                        #
                        ## store variance in array for comparison of different values of bin edges
                        SF_var_list.append([SF_var,bin_edge_means_SF[ii],num_bins_to_try[number]])
                        #print('No error. appended to list for # of bins: %s'%num_bins_to_try[number])
                        #break
                    except:
                         #
                         ## enable diagnostic output
                        if (diag_flag == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
                            if project_diagnostic_flag == 0:
                                pass
                            else:
                                print("ERROR 1: SF false pos/neg bins overlap s.t. bin edges don't increase\nmonotonically for cutoffs - spec: %s"%z_cutoff[0],";  phot: %s"%z_cutoff[1],";   for %s"%num_bins_to_try[number]," bins, method = %s"%method,'\nUSE FEWER BINS\n"NaN" appended to list.')
                            #print('error_flag=1 appended list for # of bins: %s'%num_bins_to_try[number])
                        bins_SF = np.array([float('NaN')]*len(num_bins_SF))
                        SF_var_bin_edges = float('NaN')
                        SF_var_list.append([SF_var_bin_edges,bin_edge_means_SF[ii],num_bins_to_try[number]])
                    #print('error_flag=1 appended list for # of bins: %s'%num_bins_to_try[number])
                    #
                    #
                    if (diag_flag == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
                        if project_diagnostic_flag == 0:
                            pass
                        else:
                            print('SF bin edge start: %s'%bin_edge_SF_start)
                            print('SF bin edges: %s'%bin_edge_means_SF)
                    #             
                    bin_edge_SF_start = np.round(bin_edge_SF_start+0.01,decimals=2)
                    #
                if equal_flag == 0:     # i.e. if the corresponding bin edges b/w pos/neg lists aren't equal...
                        #
                    SF_var_list.sort(key=lambda x: x[0])
                    mid_pt_to_use = SF_var_list[0][1]        # assign best midpoint to bin_edge* array 
                    if ii == 0 or ii == (len(bin_edge_means_SF)-1):   # skip edges of array
                        pass
                    else:
                        bin_edge_means_SF[ii] = mid_pt_to_use 
                    print('Bin edge to use: bin_edge_means_SF[%s]'%ii,' = %s'%bin_edge_means_SF[ii])
                    #        
                    #
                elif equal_flag == 1:
                    bin_edge_means_SF[ii] = mid_pt_to_use
                    print('Bin edge to use: bin_edge_means_SF[%s]'%ii,' = %s'%bin_edge_means_SF[ii])
                #
        #########        
#########                
#
## SF cont'd, once bin edges set
#
# make histograms
SF_pos_hist, bins_SF = np.histogram(SF_pos, bins=bin_edge_means_SF, range=range2)
SF_neg_hist, bins_SF = np.histogram(SF_neg, bins=bin_edge_means_SF, range=range2)
#
bins_SF = np.round(bins_SF,decimals=3)
#
SF_ratio = np.round((SF_pos_hist/SF_neg_hist),decimals=3)
for jj in range(len(SF_pos_hist)):
    if SF_pos_hist[jj] == 0 and SF_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
        SF_ratio[jj] = float('NaN')
#
## compute midbins for spec. mass completeness plot (i.e. plot of false pos/false neg ratios)
SF_ratio_midbins = np.round(midbins(bins_SF),decimals=3)
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
#
## add a loop to ignore the case where the number of bins exceeds the population of false pos/neg.
pop_Q_pos = len(Q_pos)
pop_Q_neg = len(Q_neg)
 #
if num_bins_to_try[number] > pop_Q_pos or num_bins_to_try[number] > pop_Q_neg:
    if diag_flag == 1 or project_diagnostic_flag == 1:
        if project_diagnostic_flag == 0:
            pass
        else:
            print('ERROR: # of bins exceeds population of Q false pos/neg for %s'%num_bins_to_try[number],' bins and cutoffs:\n%s'%z_cutoff[0],' spec and %s'%z_cutoff[1],' phot\nUSE FEWER BINS')
            Q_var = float('NaN')
#
else:
    #
    for ii in range(len(bin_edge_means_Q)):
        #
        if ii != (len(bin_edge_means_Q)-1):                    # if not the last bin edge (i.e. upper limit of mass range for study)...
            if bin_edge_means_Q[ii] > bin_edge_means_Q[(ii+1)]:   # if the next bin edge is lower than the current bin edge (in mass)...
              #  FIX THIS HERE: s/b looking at bin edges in num_bins_*, not bin_edge_means*
                #
                ## enable diagnostic output
                if (diag_flag == 1 and project_diagnostic_flag ==2) or project_diagnostic_flag == 1:
                    if project_diagnostic_flag == 0:
                        pass
                    else:
                        print("ERROR 2: the next bin value is lower than the current bin value for cutoffs - spec: %s"%z_cutoff[0],";  phot: %s"%z_cutoff[1],";   for %s"%num_bins_to_try[number]," bins, method = %s"%method,'\nThe current bin is being set to the value of the next bin.\nUSE FEWER BINS')
                #
                bin_edge_means_Q[ii] = bin_edge_means_Q[(ii+1)] 
            else:
                #
                #
                Q_var_list = []
                if num_bins_Q[0][ii] < num_bins_Q[1][ii]:    # this if statement won't apply to the first/last bin, since they have been set to equal "range2", the mass range for the study
                    if ii == 0:
                        bin_edge_Q_start = num_bins_Q[0][ii]     # store the lower bound of the ii'th bin edge b/w the false pos/neg bin edges
                    else:
                        bin_edge_Q_start = max(num_bins_Q[0][ii],bin_edge_means_Q[(ii-1)])
                    bin_edge_Q_end = num_bins_Q[1][ii]
                elif num_bins_Q[0][ii] > num_bins_Q[1][ii]:
                    if ii == 0:
                        bin_edge_Q_start = num_bins_Q[1][ii]
                    else:
                        bin_edge_Q_start = max(num_bins_Q[1][ii],bin_edge_means_Q[(ii-1)])
                    bin_edge_Q_end = num_bins_Q[0][ii]
                else:                                          # for when the i'th bin edge is the same for both false pos/neg lists
                    equal_flag = 1
                    bin_edge_Q_start = num_bins_Q[0][ii]
                    bin_edge_Q_end = num_bins_Q[0][ii]
                    mid_pt_to_use = num_bins_Q[0][ii]
                #
                if bin_edge_Q_start == bin_edge_Q_end:
                    equal_flag = 1
                else:
                    equal_flag = 0
                #
                ## enable diagnostic output
                if (diag_flag == 1 and project_diagnostic_flag ==2) or project_diagnostic_flag == 1:
                    if project_diagnostic_flag == 0:
                        pass
                    else:
                        print('Start of : %.2f'%bin_edge_Q_start,' for Q method: %s'%method,' for z_cutoff = %s'%z_cutoff)
                        #print(bin_edge_Q_start)
                        #print(bin_edge_Q_end)
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
                        Q_var = np.sum((1 - Q_ratio)**2) / len(Q_ratio)      # i.e. divide by num_bins_to_try[number]
                        #
                        ## store variance in array for comparison of different values of bin edges
                        Q_var_list.append([Q_var,bin_edge_means_Q[ii],num_bins_to_try[number]])
                        #print('No error. appended to list for # of bins: %s'%num_bins_to_try[number])
                        #break
                    except:
                         #
                         ## enable diagnostic output
                        if (diag_flag == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
                            if project_diagnostic_flag == 0:
                                pass
                            else:
                                print("ERROR 1: Q false pos/neg bins overlap s.t. bin edges don't increase\nmonotonically for cutoffs - spec: %s"%z_cutoff[0],";  phot: %s"%z_cutoff[1],";   for %s"%num_bins_to_try[number]," bins, method = %s"%method,'\nUSE FEWER BINS\n"NaN" appended to list.')
                            #print('error_flag=1 appended list for # of bins: %s'%num_bins_to_try[number])
                        bins_Q = np.array([float('NaN')]*len(num_bins_Q))
                        Q_var_bin_edges = float('NaN')
                        Q_var_list.append([Q_var_bin_edges,bin_edge_means_Q[ii],num_bins_to_try[number]])
                    #print('error_flag=1 appended list for # of bins: %s'%num_bins_to_try[number])
                    #
                    #
                    if (diag_flag == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
                        if project_diagnostic_flag == 0:
                            pass
                        else:
                            print('Q bin edge start: %s'%bin_edge_Q_start)
                            print('Q bin edges: %s'%bin_edge_means_Q)
                    #             
                    bin_edge_Q_start = np.round(bin_edge_Q_start+0.01,decimals=2)
                    #
                if equal_flag == 0:     # i.e. if the corresponding bin edges b/w pos/neg lists aren't equal...
                        #
                    Q_var_list.sort(key=lambda x: x[0])
                    mid_pt_to_use = Q_var_list[0][1]        # assign best midpoint to bin_edge* array 
                    if ii == 0 or ii == (len(bin_edge_means_Q)-1):   # skip edges of array
                        pass
                    else:
                        bin_edge_means_Q[ii] = mid_pt_to_use 
                    print('Bin edge to use: bin_edge_means_Q[%s]'%ii,' = %s'%bin_edge_means_Q[ii])
                    #        
                    #
                elif equal_flag == 1:
                    bin_edge_means_Q[ii] = mid_pt_to_use
                    print('Bin edge to use: bin_edge_means_Q[%s]'%ii,' = %s'%bin_edge_means_Q[ii])
                #
        #########        
#########                
#
## Q cont'd, once bin edges set
#
# make histograms
Q_pos_hist, bins_Q = np.histogram(Q_pos, bins=bin_edge_means_Q, range=range2)
Q_neg_hist, bins_Q = np.histogram(Q_neg, bins=bin_edge_means_Q, range=range2)
#
bins_Q = np.round(bins_Q,decimals=3)
#
Q_ratio = np.round((Q_pos_hist/Q_neg_hist),decimals=3)
for jj in range(len(Q_pos_hist)):
    if Q_pos_hist[jj] == 0 and Q_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
        Q_ratio[jj] = float('NaN')
#
## compute midbins for spec. mass completeness plot (i.e. plot of false pos/false neg ratios)
Q_ratio_midbins = np.round(midbins(bins_Q),decimals=3)
#
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
print('\n"spec_asymmetric_binning_metric.py" for z_cutoff = %s'%z_cutoff,' terminated successfully.\n')
#
#
#
#
###################     PROGRAM END     ###################    