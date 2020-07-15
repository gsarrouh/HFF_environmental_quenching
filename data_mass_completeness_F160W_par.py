#Created on Fri Jul 10 00:00:01 2020
#
#####################  data_mass_completeness_F160W_par.py  #####################
#
## THIS IS A CARBON COPY OF "data_mass_completeness_F160W*.py", but adapted to handle the parallel field arrays. Next time, turn the program into a function to be called!!!
#
#
###  This program takes as input the limiting magnitude in each cluster (ref. Shipley et al. 2018, Section 3.12 & Table 7 (p.18)). It then converts fluxes to magnitudes, and plots mass vs magnitude (in HST F160W filter) to determine the range of masses visible at the limiting magnitude for each cluster. The limiting mass is then (conservatively) the upper end of this mass range. Limiting masses are needed to determine sample sizes and membership, as well as false pos/neg populations for the spectroscopic completeness correction. 
#
### This file is a complete overhaul of data_completeness_4.py. Version 4 made plots that have since been deprecated (e.g. an initial, "hand-wavey" attempt at determining limiting magnitude) and are not reproduced here. The GOAL OF THIS FILE is to produce a mass vs magnitude plot with as little extraneous information as possible
#
### FIGURES (this file has no sections)
#
##  Fig. 1: mass-to-luminosity plot; first, convert 'f_F160W' fluxes to magnitudes. then determine limiting mass by cluster based on unique limiting magnitude of each cluster in the F160W filter; if you can't make it look pretty don't worry, at best it will be an appedix in the paper. 
##  NOTE: it may be useful to create two figures: one for a single cluster (i.e. an example, visually appealing, to demonstrate the procedure of limiting mag-to-limiting mass determination), and another a tiled plot of all 6 clusters. 
#
#
#
### NOTE: To find a section, search "SECTION (*)" in all caps, where * is the section number you want. To find fields which require user input, search "MAY NEED TO EDIT" in all caps. Some of these fields may direct you to manually enter input into a sub-program. 
#
### Section summary:
#
### PROGRAM START
#
### (0)    import modules, define functions; define LIMITING MAGNITUDES (ref. Shipley et al. 2018)
### (1)    compute MAGNITUDES, collect lists, check accounting of all selected cluster members
### (1.1)   prepare SUMMARY TABLES
### (2)    VISUALIZATION; individually (PLOT_FLAG_1) and tiled subplot (PLOT_FLAG_2) 
#
#
#
### NOTE: if DIAG_FLAG_1 IS NOT TURNED ON, THIS PROGRAM WILL NOT COMPUTE LIMITING MASSES. This is a pretty big design flaw on my part. 
#
#
#
###################     PROGRAM START     ###################
#
#
time_flag = 0     # 0= all timers off;   1=on, time entire program;    2=off, allow individual sections to be timed
#
if time_flag == 1:
    if project_time_flag == 1:
        pass
    else:
        start_time = time.time()
#
#
## SECTION (0): Modules, Definitions, FLAGS
#
## import modules & definitions
#
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
import time
#
## DEFINITIONS 
#
# the idea: feed in list of magnitudes, and find_nearest tells you the index of the one closest to the limiting magnitude. that way if the clostest objects mag. is above the limiting mag., you can just choose the next object in the list, since the list is sorted. 
def find_nearest(list, value):
    array = np.asarray(list)
    index = (np.abs(array - value)).argmin()
    return index
#
## for finding the max/min value in a nested list where each row has X values to unpack
def max_nested_list(nested_list,axis):
    max_value = 0
    index_returned = 0
    for index in range(len(nested_list)):
        item = nested_list[index][axis]
        if item > max_value:
            max_value = item
            index_returned = index
    return nested_list[index_returned], index
#
def min_nested_list(nested_list,axis):
    min_value = 1e10
    index_returned = 0
    for index in range(len(nested_list)):
        item = nested_list[index][axis]
        if item < min_value:
            min_value = item
            index_returned = index
    return nested_list[index_returned], index
#
#def std_nested_list(nested_list,axis):
#    std_array = np.array([0]*len(nested_list))
#    for ii in range(len(nested_list)):
#        std_array[ii] = nested_list[ii][axis]
#    std = np.std(std_array)
#
#
limiting_mag_par = np.array([27.3,27.6,26.4,27.3,27.8,27.1])    # NOTE: order is [M0416,M1149,M0717,A370,A1063,A274]
cluster_names_par = ['M0416','M1149','M0717','A370','A1063','A2744']   # in the order corresponding to limiting mags above
#
## FLAGS - search flag name to find it, or search "MAY NEED TO EDIT"
#
summary_flag = 1          # print summary table comparing total & cluster members here to those in 'master_data*.py'
#
diag_flag_1 = 0           # display cluster stats and comparison w/ "master_data*.py" analysis;
diag_flag_2 = 1           # display result: limiting mag, limiting mass, magnitude of object chosen as limiting mass
diag_flag_3 = 0           # prints the indices of objects w/ "NaN" mass estimates;  
diag_flag_4 = 1           # Count galaxies sorted through at each step of initial 'for' loop
diag_flag_5 = 0           # AVAILABLE
#
plot_flag_1 = 0           # produce Fig.1: mass vs mag for 1 cluster at a time (plots 6 figures total);
plot_flag_2 = 0           # produce Fig.2: tiled subplot of all 6 clusters zoomed in (single figure);   
#
#
#
## SECTION (1): compute MAGNITUDES & perform accounting of all objects in catalogue, specifically selected members
#
## start by adding a column to master_cat store magnitudes, if it doesn't already exist (add try/except for repeatablility)
#
try: 
    F = Column([-99]*len(master_cat_par),name='F160W_mag',dtype='float32')   # add a column in which to store magnitudes
    master_cat_par.add_column(F)
    print("This is the first time you ran this program today, isn't it ghassan?")
except:
    pass
#
#                **** ****                **** ****      LOOK          **** ****                **** ****
#
## I orginally did the following in two loops, but it made more sense to just combine it in one. unfortuantely, the documentation got a bit unruly, and some of the below comments are deprecated and i don't want to go through it to update. IN A NUTSHELL: [mag,mass] pairs are stored as lists in "mag_by_cluster_member", a list of lists (i.e. one list for each cluster, the entries of which are the paired lists [mag,mass] just described). these are what will be scatter-plotted by cluster for visualization. "masses_at_lim_mag" is a dedicated list to storing [mag,mass] pairs which are within +/- "TOL" of the limiting magnitude for that cluster (again, stores each [mag,mass] pair in one of 6 cluster lists). the limiting mass is the max. mass in this last list, for each individual cluster. 

## UPDATE: to keep track of the total catalogue population and cluster members separately, I went back to two separate loops.
#
#                **** ****                **** ****      HERE          **** ****                **** ****
#
#
## to store RESULT of program
limiting_mass_par = np.array([0]*6,dtype='float32')                    
mag_of_limiting_mass_par = np.array([0]*6,dtype='float32')             # to store magnitude of the object which sets the limiting mass
z_of_limiting_mass_par =  np.array([0]*6,dtype='float32')             # to store del_z photometric redshift difference ('z_clusterphot') of the object which sets the limiting mass (for comparison w/ redshift cuts)
#
masses_at_lim_mag_par =  [ [], [], [], [], [], [] ]    # array to store masses at (or near) the limiting magnitude, to determine a range of masses corresponding to the cluster's limiting mag (along w/ mag for a diagnostic check)
#
## MAY NEED TO EDIT: tolerance ("TOL")
TOL = 0.05                                         # this is how near the limiting mag. i'm willing to accept an object
#
mag_by_cluster_par = [ [], [], [], [], [], [] ]        # two lists: 1-all objects in the ii'th image; 1-cluster members only
mag_by_cluster_member_par = [ [], [], [], [], [], [] ] # each list (col) contains lists with 3 entries: [mag,mass,z_phot]
# for storing info on stats (diagnostic)
avg_mag_by_cluster_par = np.array([0.0]*6)               # diagnostic
max_mag_by_cluster_par = np.array([0.0]*6)               # i.e. store the dimmest (least bright) source
# for storing the max/min masses at the limiting magnitude of each cluster (for plotting)
max_min_mass_par = np.array([[0]*2]*6,dtype='float32')       # [min,max]
#
mag_spec1_par = np.array([0]*6) 
mag_phot1_par = np.array([0]*6) 
bad_flag1_par = np.array([0]*6)                         # track how many objects are flagged by cluster
bad_flux1_par = np.array([0]*6)                         # track # of objects w/ bad fluxes not picked up by 'flag_F160W'
count_nans1_par = np.array([0]*6)                       # keep track of objects w/ NaN mass estimates by cluster (entire image)
count_non_member_par = np.array([0]*6) 
#
counting_array1_par = np.array([0]*9)
#
#
for counter in range(len(master_cat_par)):             # loop through catalogue
    if master_cat_par['flag_F160W'][counter] == 0:         # avoid objects w/ bad phot; good galaxies marked by flag_band == 0
        counting_array1_par[0]+=1  # good flag galaxies
        if master_cat_par['f_F160W'][counter] > 0:         # double check, since bad phots are set to == -99
            counting_array1_par[1]+=1  # good flux galaxies
            master_cat_par['F160W_mag'][counter] = ((-2.5*np.log10(master_cat_par['f_F160W'][counter]))+25)     # AB_mag system zero point = 25, confirmed in Shipley et al. 
            if (diag_flag_3 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:   # diagnostic check
                if project_diagnostic_flag == 0:
                    pass
                else:
                    if np.isnan(master_cat_par['F160W_mag'][counter]) == 1:      # to identify bad mass estimates (i.e NaNs)
                        print('Diag. 3: NaNs in mag estimates\nNan index in master_cat (mag) = %s'%counter)
                    elif np.isnan(master_cat_par['lmass'][counter]) == 1:
                        print('Diag. 3: NaNs in mass estimates\nNan index in master_cat (mass) = %s'%counter)
                        #
            #
            for ii in range(len(limiting_mag_par)):
                if master_cat_par['cluster'][counter] == (ii+1):                 # track objects by cluster
                    if np.isnan(master_cat_par['lmass'][counter]) == 1:          # track all nans in IMAGE
                        count_nans1_par[ii]+=1
                        counting_array1_par[2]+=1 # galaxies w/ bad mass estimates
                    else:
                        #masses[ii].append(master_cat['lmass'][counter])      # store all the masses (regardless of type)
                        mag_by_cluster_par[ii].append([master_cat_par['F160W_mag'][counter],master_cat_par['lmass'][counter],master_cat_par['z_peak'][counter]])     # look at all objects in the ii'th image
                        if np.isnan(master_cat_par['lmass'][counter]) == 1:      # track all nans in CLUSTER
                            counting_array1_par[3]+=1
                        else:
                            if master_cat_par['member'][counter] == 1:    # now look only at cluster MEMBERS only  
                                counting_array1_par[4]+=1
                                if master_cat_par['sub'][counter] == 1:     # (spec+phot)
                                    mag_spec1_par[ii]+=1
                                elif master_cat_par['sub'][counter] == 2:   # phot only
                                    mag_phot1_par[ii]+=1
                                #
                                #
                                ###  if you want to track SF/Q TYPE galaxies, add a list here
                                #
                                #
                                ## now track objects at the limiting magnitude of the cluster
                                #
                            else:
                                counting_array1_par[8]+=1   # all non-member galaxies            
                                for ii in range(len(limiting_mag_par)):
                                    if master_cat_par['cluster'][counter] == (ii+1):
                                        count_non_member_par[ii]+=1
                            #  
        # 
        else:
            for cluster in range(len(limiting_mag_par)):
                if master_cat_par['cluster'][counter] == (cluster+1):
                    bad_flux1_par[cluster]+=1
                    counting_array1_par[5]+=1    # bad flux (<0)
    else:
        counting_array1_par[6]+=1     # bad flag ALL
        for cluster in range(len(limiting_mag_par)):
            if master_cat_par['cluster'][counter] == (cluster+1):
                bad_flag1_par[cluster]+=1
        for cluster in range(len(limiting_mag_par)):
            if master_cat_par['cluster'][counter] == (cluster+1):
                if master_cat_par['member'][counter] == 0:     # we're only interesting in comparing for cluster members
                    #
                    counting_array1_par[7]+=1     # bad flag member
            #
#
if (diag_flag_4 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    print('\n[good_flag,good_flux,NaN-mass,NaN-mass(member),member,bad_flux,bad_flag,bad_flag_member,non_member]:\n%s'%counting_array1_par)
#
#
mag_spec2_par = np.array([0]*6) 
mag_phot2_par = np.array([0]*6) 
bad_flag2_par = np.array([0]*6) 
bad_flux2_par = np.array([0]*6) 
outside_parent_sample_par = np.array([0]*6) 
count_nans2_par = np.array([0]*6)                       # keep track of objects w/ NaN mass estimates by cluster (entire image)
counting_array2_par = np.array([0]*6)
#
for counter in range(len(master_cat_par)):
    if master_cat_par['member'][counter] == 1:                        # number of galaxies in desired redshift range
        counting_array2_par[0]+=1                                     # number of members
        if master_cat_par['flag_F160W'][counter] == 0:                # good flag
            counting_array2_par[1]+=1                       
            if master_cat_par['f_F160W'][counter] > 0:                # good flux
                counting_array2_par[2]+=1
                if np.isnan(master_cat_par['lmass'][counter]) == 1:
                    for ii in range(len(limiting_mag_par)):
                        if master_cat_par['cluster'][counter] == (ii+1):
                            count_nans2_par[ii]+=1
                elif master_cat_par['sub'][counter] == 1:     # (spec+phot)
                    for ii in range(len(limiting_mag_par)):
                        if master_cat_par['cluster'][counter] == (ii+1):
                            mag_by_cluster_member_par[ii].append([master_cat_par['F160W_mag'][counter],master_cat_par['lmass'][counter],master_cat_par['z_peak'][counter]])
                            mag_spec2_par[ii]+=1
                            if master_cat_par['F160W_mag'][counter] < (limiting_mag_par[ii]+TOL) and master_cat_par['F160W_mag'][counter] > (limiting_mag_par[ii]-TOL):           # if flux is equal to limiting magnitude, +/- 'TOL', store mass value...
                                masses_at_lim_mag_par[ii].append([master_cat_par['F160W_mag'][counter],master_cat_par['lmass'][counter],master_cat_par['z_peak'][counter]])
                                    #
                elif master_cat_par['sub'][counter] == 2:   # phot only
                    for ii in range(len(limiting_mag_par)):
                        if master_cat_par['cluster'][counter] == (ii+1):
                            mag_by_cluster_member_par[ii].append([master_cat_par['F160W_mag'][counter],master_cat_par['lmass'][counter],master_cat_par['z_peak'][counter]])
                            mag_phot2_par[ii]+=1
                            if master_cat_par['F160W_mag'][counter] < (limiting_mag_par[ii]+TOL) and master_cat_par['F160W_mag'][counter] > (limiting_mag_par[ii]-TOL):           # if flux is equal to limiting magnitude, +/- 'TOL', store mass value...
                                masses_at_lim_mag_par[ii].append([master_cat_par['F160W_mag'][counter],master_cat_par['lmass'][counter],master_cat_par['z_peak'][counter]])
                                #
                        #
                else:                                                     # outside parent sample
                    for ii in range(len(limiting_mag_par)):
                        if master_cat_par['cluster'][counter] == (ii+1):  # outside parent sample
                            outside_parent_sample_par[ii]+=1
            else:
                for ii in range(len(limiting_mag_par)):
                    if master_cat_par['cluster'][counter] == (ii+1):  # bad flux
                        bad_flux2_par[ii]+=1
                counting_array2_par[4]+=1                            
        else:
            for ii in range(len(limiting_mag_par)):
                if master_cat_par['cluster'][counter] == (ii+1):      # bad flag
                    bad_flag2_par[ii]+=1
            counting_array2_par[5]+=1                           
            
#
if (diag_flag_4 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    print('\n[members,good_flag,good_flux,not_in_parent,bad_flux,bad_flag]:\n%s'%counting_array2_par)
#
## DIAG_FLAG_4: ACCOUNTING for each galaxy at each condition imposed within the above for loop, which sets up the lists of galaxies to be plotted
if (diag_flag_4 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    print('\n"flag_F160W==0" OK: %s'%counting_array1_par[0],'\n"flag_F160W!=0" BAD: %s'%counting_array1_par[6],'\nGood flag + Bad flag= %s'%counting_array1_par[0],' + %s'%counting_array1_par[6],' = %s'%np.sum([counting_array1_par[0],counting_array1_par[6]]))
    print('\nOf the "Good flag" galaxies\nflux_F160W>0: %s'%counting_array1_par[1],'\nflux_F160W<0: %s'%counting_array1_par[5],'\n(flux>0) + (flux<0) = %s'%counting_array1_par[1],' + %s'%counting_array1_par[5],' = %s'%(np.sum([counting_array1_par[1],counting_array1_par[5]])))
    print('\nlmass = NaN (all): %s'%counting_array1_par[2],'\nCluster members selected: %s'%counting_array1_par[4],'\nlmass = NaN (members): %s'%counting_array1_par[3],'\nBad flag: %s'%counting_array1_par[7],'\nBad flux (<0): %s'%counting_array1_par[5])

#
for ii in range(len(limiting_mag_par)):
    mag_by_cluster_par[ii].sort(key=lambda x: x)                # sort data by magnitude within each clust
    mag_by_cluster_member_par[ii].sort(key=lambda x: x)
    # temp_*_array(s) store the contents of "mag_by_cluster_*" for a single cluster
    temp_mag_array_par = np.array([[0]*3]*len(mag_by_cluster_par[ii]),dtype='float32')        # find max magnitude of image & cluser
    temp_mag_member_array_par = np.array([[0]*3]*len(mag_by_cluster_member_par[ii]),dtype='float32')
    # similar to temp*arrays, but dedicated for magnitude
    store_mags_par = np.array([0]*len(mag_by_cluster_par[ii]),dtype='float32') 
    store_mags_member_par = np.array([0]*len(mag_by_cluster_member_par[ii]),dtype='float32') 
    #
    for jj in range(len(mag_by_cluster_par[ii])):
        temp_mag_array_par[jj][0] = mag_by_cluster_par[ii][jj][0]            # assign fluxes from a given cluster to a dedicated temporary list
        temp_mag_array_par[jj][1] = mag_by_cluster_par[ii][jj][1]
        temp_mag_array_par[jj][2] = mag_by_cluster_par[ii][jj][2]
        store_mags_par[jj] = mag_by_cluster_par[ii][jj][0]                   # magnitude
    for jj in range(len(mag_by_cluster_member_par[ii])):
        temp_mag_member_array_par[jj][0] = mag_by_cluster_member_par[ii][jj][0]  # do the same for cluster members - MAG
        temp_mag_member_array_par[jj][1] = mag_by_cluster_member_par[ii][jj][1]  # do the same for cluster members - MASS
        temp_mag_member_array_par[jj][2] = mag_by_cluster_member_par[ii][jj][2]  # do the same for cluster members - z_peak (i.e. z_phot)
        store_mags_member_par[jj] = mag_by_cluster_member_par[ii][jj][0]                   # magnitude
    #
    if len(masses_at_lim_mag_par[ii]) == 0:                     # if no objects were found within the tolerance range...
        temp_list_par = []
        for jj in range(len(temp_mag_member_array_par)):
            if temp_mag_member_array_par[jj][0] <= limiting_mag_par[ii]:
                temp_list_par.append([temp_mag_member_array_par[jj][0],temp_mag_member_array_par[jj][1],temp_mag_member_array_par[jj][2]])   #[mag,mass,redshift]
        temp_limit_par = min(temp_list_par, key=lambda x: abs(x[0]-limiting_mag_par[ii])) 
        while temp_limit_par[0] > limiting_mag_par[ii]:
            print('Closest limiting mass found was above limiting mag. Looking for a new one...')
            index = np.where(temp_mag_member_array_par == temp_limit_par[0])
            temp_mag_member_array_par = np.delete(temp_mag_member_array_par, (index), axis=0)
            temp_limit_par = min(temp_mag_member_array_par, key=lambda x: abs(x[0]-limiting_mag_par[ii]))
        print('WARNING: NO GALAXIES FOUND AT CLUSTER %s'%(ii+1),' LIMITING MAG. of %s'%limiting_mag_par[ii],' +/- %.2f'%TOL,'\nClosest thing Old Faithful could find was %.2f'%temp_limit_par[1],' w/ mag %.2f'%temp_limit_par[0])
        limiting_mass_par[ii] = temp_limit_par[1]  
        mag_of_limiting_mass_par[ii] = temp_limit_par[0]
        z_of_limiting_mass_par[ii] = temp_limit_par[2]
        max_min_mass_par[ii][0] = temp_limit_par[1]               # for PLOTTING
        max_min_mass_par[ii][1] = temp_limit_par[1]
    elif len(masses_at_lim_mag_par[ii]) == 1:
        limiting_mass_par[ii] = masses_at_lim_mag_par[ii][0][1]
        mag_of_limiting_mass_par[ii] = masses_at_lim_mag_par[ii][0][0]
        z_of_limiting_mass_par[ii] = masses_at_lim_mag_par[ii][0][2]
        max_min_mass_par[ii][0] = masses_at_lim_mag_par[ii][0][1]               # for PLOTTING
        max_min_mass_par[ii][1] = masses_at_lim_mag_par[ii][0][1]
    else:    
        the_good_shit_par, the_good_shit_index_par = max_nested_list(masses_at_lim_mag_par[ii],1)
        the_not_quite_as_good_shit_par, the_not_quite_as_good_shit_index_par = min_nested_list(masses_at_lim_mag_par[ii],1)
        limiting_mass_par[ii] = the_good_shit_par[1] # store highest mass you have
        mag_of_limiting_mass_par[ii] = the_good_shit_par[0] # store highest mass you have
        z_of_limiting_mass_par[ii] = the_good_shit_par[2]
        max_min_mass_par[ii][0] = the_good_shit_par[1]               # for PLOTTING
        max_min_mass_par[ii][1] = the_not_quite_as_good_shit_par[1]
    if (diag_flag_1 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
        if project_diagnostic_flag == 0:
            pass
        else:
            avg_mag_temp_par = np.mean(store_mags_par)                        # report some useful info to the user
            max_mag_temp_par = max(store_mags_par)
            avg_mag_temp_member_par = np.mean(store_mags_member_par)
            max_mag_temp_member_par = max(store_mags_member_par)
            # now print the above
            print('\nDiag. 1: compare w/ "master_data*.py" result\nCluster %s'% (ii+1),': FULL IMAGE STATS\nMax. mag (i.e. dimmest source in ENTIRE image) %.3f'%max_mag_temp_par,', Avg. mag of IMAGE: %.3f'%avg_mag_temp_par,'\n# of objects w/ "NaN" mass estimates: %s'%count_nans_par[ii])
            print('Cluster %s'%(ii+1),': MEMBER-ONLY STATS\nMax. mag (i.e. dimmest source in CLUSTER only) %.3f'%max_mag_temp_member_par,', Avg. mag of MEMBERS: %.3f'%avg_mag_temp_member_par,'\n# of objects w/ "NaN" mass estimates: %s'%count_nans_member_par[ii])
            #print('\n"master_data*.py" reports (SF+Q): %s'%((mem[0][ii]+mem[1][ii])),' (spec) cluster members & %s'%((mem_phot[0][ii]+mem_phot[1][ii])),' (phot). Total: %s'%((mem[0][ii]+mem[1][ii])+(mem_phot[0][ii]+mem_phot[1][ii])),' for cluster %i'%(ii+1))
            #print('"master_data*.py" reports total: %s'%((mem[0][ii]+mem_phot[0][ii])),' (SF) cluster members & %s'%((mem[1][ii]+mem_phot[1][ii])),' (Q).')
            print('\nTHIS PROGRAM reports: %i'%(len(mag_by_cluster_member_par[ii])),' total cluster members for cluster %i'%(ii+1),'\n%i'%bad_flag_par[ii],' flagged objects\n%i'%bad_flux_par[ii],' objects with bad fluxes (<0)\n%i'%count_nans_member_par[ii],' NaNs\nAll this action happened in cluster %i'%(ii+1),', yo.\n')
#        
#
if diag_flag_2 == 1:
    print('Limiting magnitudes per Shipley et a. 2018 (PAR): %s'%limiting_mag_par)
    print('Limiting masses: %s'%limiting_mass_par)
    print('Magnitudes at limiting mass: %s'%mag_of_limiting_mass_par)
#
#
## The code now handles finding the highest mass at the limiting magnitude (or the mass of the object closest in mmagnitude to the limit).
#
#
#
## SECTION (1.1): SUMMARY TABLE
#
## Summary Table
#
#
if summary_flag == 1 or adams_flag == 1:
    ## Summarize limiting mass calculation in table for PARENT SAMPLE
    lim_mass_par_names = Column(['Full CATALOGUE ("master_parallel*.py")','Phot Members','Spec Members','Non-members','Bad flag (!=0)','Bad flux (F160W<0)','Mass = NaN','SUM'],name='Property')
    col_names = cluster_names_par
    lim_mass_par0 = Column([np.sum([phot_only_par,both_par,spec_only_par,no_data_par,stars_sub_par]),np.sum(mag_phot1_par),np.sum(mag_spec1_par),np.sum(count_non_member_par),np.sum(bad_flag1_par),np.sum(bad_flux1_par),np.sum(count_nans1_par),np.sum([mag_phot1_par,mag_spec1_par,bad_flag1_par,bad_flux1_par,count_nans1_par,count_non_member_par])],name='Total')  # total column
    lim_mass_par_stats = Table([lim_mass_par_names,lim_mass_par0])
    for ii in range(len(mag_phot1_par)):
        lim_mass_par_col = Column([np.sum([phot_only_par[ii],both_par[ii],spec_only_par[ii],no_data_par[ii],stars_sub_par[ii]]),mag_phot1_par[ii],mag_spec1_par[ii],count_non_member_par[ii],bad_flag1_par[ii],bad_flux1_par[ii],count_nans1_par[ii],np.sum([mag_phot1_par[ii],mag_spec1_par[ii],bad_flag1_par[ii],bad_flux1_par[ii],count_nans1_par[ii],count_non_member_par[ii]])],name=col_names[ii])               # cluster columns
        lim_mass_par_stats.add_column(lim_mass_par_col) 
    #
    #
    #
    ## Summarize limiting mass calculation in table for CLUSTER MEMBERS
    lim_mass_member_par_names = Column(['Members ("master_parallel*.py")','Phot Members','Spec Members','Outside Phot. Parent Sample','Bad flag (!=0)','Bad flux (F160W<0)','Mass = NaN','SUM'],name='Property')
    col_names = cluster_names_par
    lim_member_mass_par0 = Column([np.sum([np.sum(count_field_sample)]),np.sum(mag_phot2_par),np.sum(mag_spec2_par),np.sum(outside_parent_sample_par),np.sum(bad_flag2_par),np.sum(bad_flux2_par),np.sum(count_nans2_par),np.sum([mag_phot2_par,mag_spec2_par,bad_flag2_par,bad_flux2_par,count_nans2_par])],name='Total')  # total column
    lim_mass_member_par_stats = Table([lim_mass_member_par_names,lim_member_mass_par0])
    for ii in range(len(mag_phot1_par)):
        lim_mass_member_par_col = Column([np.sum([np.sum(count_field_sample_type[0][ii]),np.sum(count_field_sample_type[1][ii])]),mag_phot2_par[ii],mag_spec2_par[ii],outside_parent_sample_par[ii],bad_flag2_par[ii],bad_flux2_par[ii],count_nans2_par[ii],np.sum([mag_phot2_par[ii],mag_spec2_par[ii],bad_flag2_par[ii],bad_flux2_par[ii],count_nans2_par[ii],outside_parent_sample_par[ii]])],name=col_names[ii])               # cluster columns
        lim_mass_member_par_stats.add_column(lim_mass_member_par_col) 
    #
#
    print('\nFull CATALOGUE - PAR - breakdown\n%s'%lim_mass_par_stats)
    print('\nMember sample - PAR - breakdown\n%s'%lim_mass_member_par_stats)
#
#
#
#
## SECTION (2) VISUALIZATION
#
## define range of all histograms to come!
#
max_mass_par_array = [ [], [], [], [], [], [] ]
max_mass_index_par = np.array([0]*6)
max_mass_par = np.array([0]*1)
#
for ii in range(len(z_cluster)):
    max_mass_par_array[ii], max_mass_index_par[ii] = max_nested_list(mag_by_cluster_member_par[ii],1)
max_mass_par, max_mass_par_index = max_nested_list(max_mass_par_array,1)
#
range2_par=[min(limiting_mass_par),max_mass_par[1]]
#
## Now I have lists of corresponding magnitudes and masses for each cluster. make a single plot for the first cluster (as a clear example of the method), then a subplot (w/ 6 plots) all together. add lines to the plot for (A) the limiting magnitude (vertical); and (B) range of masses at that magnitude (horizontal), for clusters w/ many such objects, or circle/point to/identify in some way on the plot the galaxy closest to the limiting magnitude which was chosen as the standard for limiting mass
#
## Fig 1: mass vs magnitude scatter plot, colourbar = redshift (z_phot)
## Now make the single scatter plot for M0416 (i.e. cluster stored in the first position of the above lists)
#
## 
## setup dedicated arrays for plotting, since plotting them from my nested list of lists of lists doesn't seem to be working...
#
if plot_flag_1 == 1:
    for cluster in range(len(limiting_mag_par)):
        plotting_array_temp_par = np.array([[0]*len(mag_by_cluster_member_par[cluster])]*3,dtype='float32') # 0=mag;1=mass;2=z
        for ii in range(len(mag_by_cluster_member_par[cluster])):
            plotting_array_temp_par[0][ii] = mag_by_cluster_member_par[cluster][ii][0]   # store mass/magnitudes for the first cluster
            plotting_array_temp_par[1][ii] = mag_by_cluster_member_par[cluster][ii][1]   # index order: [cluster][object][mag/mass,z_phot]([cluster,row,col])
            plotting_array_temp_par[2][ii] = mag_by_cluster_member_par[cluster][ii][2]
        #
        ## calculate cluster members; display as "Total(Shown)"
        total_mem_par = np.sum([count_field_sample[0][cluster],count_field_sample[1][cluster]])
        mem_shown_par = len(plotting_array_temp_par[0])
        #
        members_string_par = '# members(F160W): %s'%total_mem_par+'(%s'%mem_shown_par+')'
        ## Compute std. dev. of cluster objects
        std_dev_par = np.std(plotting_array_temp_par[1])
        #
        ## The figure:
        #
        c_map = plt.get_cmap("spring")        # set colorbar map
        #
        #    plt.close()
        fig, ax = plt.subplots()
        ax.scatter(plotting_array_temp_par[0],plotting_array_temp_par[1], marker='o', s=10, c=plotting_array_temp_par[2], cmap= c_map, vmin=min(plotting_array_temp_par[2]), vmax=max(plotting_array_temp_par[2]))#, edgecolors='k')#,facecolors='none', edgecolors='k')
        #ax.scatter(SF_masses[0][0],SF_masses[0][1], marker='*', s=10, facecolors='none', edgecolors='b')
        ax.plot([limiting_mag_par[cluster],limiting_mag_par[cluster]],[0,15],'--r', linewidth=1.4)
        ax.plot([0,35],[max_min_mass_par[cluster][0],max_min_mass_par[cluster][0]], color='maroon', linestyle='--',linewidth=1.2)    #chang indexing of max_min_mass to [ii][0] & [ii][1] for looped subplots
        ax.plot([0,35],[max_min_mass_par[cluster][1],max_min_mass_par[cluster][1]], color='maroon', linestyle='--',linewidth=1.2)
        sm =  ScalarMappable(cmap=c_map)
        sm.set_array([])
        sm.set_clim([0,z_cutoff[1]])
        cbar = fig.colorbar(sm, ax=ax)
        cbar.ax.set_title("$z$_{phot,F160W}$")
        #ax.pcolormesh(plotting_array_temp[0], plotting_array_temp[1], plotting_array_temp[2], cmap=c_map)
        #ax.colorbar()

        #plt.plot([0,35],[max_mass,max_mass], '-.k', linewidth=0.8)
        #plt.plot([0,35],[min_mass,min_mass], '-.k', linewidth=0.8)
        ax.set_xlabel('$m_{F160W}$')
        ax.set_xlim(17,30)
        ax.set_ylabel('$log(M/M_{\odot})$')
        ax.set_ylim(5,13)
        ax.grid(b=False, which='major', axis='both', color = 'k', linestyle = ':')
        ax.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=False, labelleft=True)
        ax.minorticks_on()
        #
        plt.text(21.6,12.1,cluster_names_par[cluster],fontsize=12)
        plt.text(21.6,11.5,'$M/L_{max}$: %s'%limiting_mass_par[cluster],fontsize=10)
        plt.text(21.6,10.5,members_string_par,fontsize=10)
        plt.text(18.1,6.1,'0.3 < z < 0.55',fontsize=10)
        plt.text(18.1,5.7,'$\sigma_{mass}$ = %s'%std_dev_par,fontsize=10)
        #
        plt.show()
    #
#
#
#
## Visualize by cluster
#
if plot_flag_2 == 1:
    #
    fig, axs = plt.subplots(nrows=2,ncols=3,sharex=True,sharey=True)#,tight_layout=True)
    fig.subplots_adjust(wspace=0,hspace=0)
    #
    for cluster in range(len(limiting_mag_par)):
        plotting_array_temp_par = np.array([[0]*len(mag_by_cluster_member_par[cluster])]*3,dtype='float32') # 0=mag;1=mass;2=z
        for ii in range(len(mag_by_cluster_member_par[cluster])):
            plotting_array_temp_par[0][ii] = mag_by_cluster_member_par[cluster][ii][0]   # store mass/magnitudes for the first cluster
            plotting_array_temp_par[1][ii] = mag_by_cluster_member_par[cluster][ii][1]   # index order: [cluster][object][mag/mass,z_phot]([cluster,row,col])
            plotting_array_temp_par[2][ii] = mag_by_cluster_member_par[cluster][ii][2]
        #
        ## calculate cluster members; display as "Total(Shown)"
        total_mem_par = np.sum([mem_phot_par[0][cluster],mem_phot_par[1][cluster],mem_spec_par[0][cluster],mem_spec_par[1][cluster]])
        mem_shown_par = len(plotting_array_temp_par[0])
        #
        members_string_par = '# cluster members(F814W): %s'%total_mem_par+'(%s'%mem_shown_par+')'
        ## Compute std. dev. of cluster objects
        std_dev_par = np.std(plotting_array_temp_par[1])
        #
        ## The figure:
        #
        c_map = plt.get_cmap("spring")        # set colorbar map
        #
        lim_mag_minus_par = (limiting_mag_par[cluster]-TOL)
        lim_mag_plus_par = (limiting_mag_par[cluster]+TOL)
        #
        #
        ax = axs.flat[cluster]
        #
        ax.scatter(plotting_array_temp_par[0],plotting_array_temp_par[1], marker='o', s=30, c=plotting_array_temp_par[2], cmap= c_map, vmin=min(plotting_array_temp_par[2]), vmax=max(plotting_array_temp_par[2]))
        ax.plot([limiting_mag_par[cluster],limiting_mag_par[cluster]],[0,15],'--k', linewidth=1.4)
        ax.plot([lim_mag_minus_par,lim_mag_minus_par],[0,15],':r', linewidth=1.2)
        ax.plot([lim_mag_plus_par,lim_mag_plus_par],[0,15],':r', linewidth=1.2)
        ax.plot([0,35],[max_min_mass_par[cluster][0],max_min_mass_par[cluster][0]], color='maroon', linestyle='--',linewidth=1.2)    #chang indexing of max_min_mass to [ii][0] & [ii][1] for looped subplots
        ax.plot([0,35],[max_min_mass_par[cluster][1],max_min_mass_par[cluster][1]], color='maroon', linestyle='--',linewidth=1.2)
        sm =  ScalarMappable(cmap=c_map)
        sm.set_array([])
        sm.set_clim([0,z_cutoff[1]])
        ax.grid(axis='both', alpha=0.75)
        ax.set_xlim([25,27])
        ax.set_ylim([5.5,8.5])
        #
        #plt.text(21.6,11.5,'$M/L_{max}$: %s'%limiting_mass[cluster],fontsize=10)
    ## label locations
    axs.flat[0].text(25.15,6.05,cluster_names_par[0],fontsize=12)
    axs.flat[0].text(25.15,5.75,'$M/L_{max,F160W}$: %s'%limiting_mass_par[0],fontsize=12)
    axs.flat[1].text(25.15,6.05,cluster_names_par[1],fontsize=12)
    axs.flat[1].text(25.15,5.75,'$M/L_{max,F160W}$: %s'%limiting_mass_par[1],fontsize=12)
    axs.flat[2].text(25.15,6.05,cluster_names_par[2],fontsize=12)
    axs.flat[2].text(25.15,5.75,'$M/L_{max,F160W}$: %s'%limiting_mass_par[2],fontsize=12)
    axs.flat[3].text(25.15,6.05,cluster_names_par[3],fontsize=12)
    axs.flat[3].text(25.15,5.75,'$M/L_{max,F160W}$: %s'%limiting_mass_par[3],fontsize=12)
    axs.flat[4].text(25.15,6.05,cluster_names_par[4],fontsize=12)
    axs.flat[4].text(25.15,5.75,'$M/L_{max,F160W}$: %s'%limiting_mass_par[4],fontsize=12)
    axs.flat[5].text(25.15,6.05,cluster_names_par[5],fontsize=12)
    axs.flat[5].text(25.15,5.75,'$M/L_{max,F160W}$: %s'%limiting_mass_par[5],fontsize=12)
    #axs.clim(0,z_cutoff[1])
    cbar = fig.colorbar(sm, ax=axs[:, 2],location='right')#, vmin=0.0, vmax=z_cutoff[1])
    cbar.ax.set_title("$z_{phot,F160W}$")
    #            
    #            
    plt.show()
#
#
#
#
## TIME_FLAG END
#
if time_flag == 1:
    if project_time_flag == 1:
        pass
    else:
        print('Program "data_mass_completeness_F160W_par.py" took: %s seconds to run.\n\n' % (time.time() - start_time))
#
#
#
#
print('\n\n"data_mass_completeness_F160W_par.py"  terminated successfully.\n')
#
#                        
###################     PROGRAM END     ###################