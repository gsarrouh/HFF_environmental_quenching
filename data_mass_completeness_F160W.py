#Created on Tue Jun 30 13:58:21 2020
#
#####################  data_mass_completeness_F160W.py  #####################
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
### (2)    determine LIMITING MASS
### (3)    VISUALIZATION; individually (PLOT_FLAG_1) and tiled subplot (PLOT_FLAG_2) 
### (4)    REMOVE all galaxies below limiting mass from prelimenary cluster member sample (i.e. reclassify 
### (4.1)   prepare SUMMARY TABLES
#
### PROGRAM END
#
#
### NOTE: if DIAG_FLAG_1 IS NOT TURNED ON, THIS PROGRAM WILL NOT COMPUTE LIMITING MASSES. This is a pretty big design flaw on my part. 
#
#
#
###################     PROGRAM START     ###################
#
#
time_flag = 1     # 0= all timers off;   1=on, time entire program;    2=off, allow individual sections to be timed
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
## FLAGS - search flag name to find it, or search "MAY NEED TO EDIT"
#
summary_flag = 1          # print summary table comparing total & cluster members here to those in 'master_data*.py'
#
diag_flag = 1
diag_flag_1 = 1           # display cluster stats and comparison w/ "master_data*.py" analysis;
diag_flag_2 = 1           # display result: limiting mag, limiting mass, magnitude of object chosen as limiting mass
diag_flag_3 = 0           # prints the indices of objects w/ "NaN" mass estimates;  
diag_flag_4 = 1           # Count galaxies sorted through at each step of initial 'for' loop
diag_flag_5 = 0           # AVAILABLE
#
plot_flag_1 = 0           # produce Fig.1: mass vs mag for 1 cluster at a time (plots 6 figures total);
plot_flag_2 = 0           # produce Fig.2: tiled subplot of all 6 clusters zoomed in (single figure);   
plot_flag_3 = 0           # produce Fig.3: PUBLICATION FIGURE: left panel: zoomed out; right lanel: zoomed in;   
#
#
#
## define "global" arrays to be used ubiquitouslsy throughout program
limiting_mag = np.array([26.2,26.5,25.5,26.1,26.5,26.7])    # NOTE: order is [M0416,M1149,M0717,A370,A1063,A274]
limiting_mag_814 = np.array([27.5,27.4,26.9,27.1,27.4,27.2])
cluster_names = ['M0416','M1149','M0717','A370','A1063','A2744']   # in the order corresponding to limiting mags above
#
#
#
## SECTION (1): compute MAGNITUDES & perform accounting of all objects in catalogue, specifically selected members
#
## start by adding a column to master_cat store magnitudes, if it doesn't already exist (add try/except for repeatablility)
#
try: 
    F = Column([-99]*len(master_cat),name='F160W_mag',dtype='float32')   # add a column in which to store magnitudes
    master_cat.add_column(F) 
    G = Column([-99]*len(master_cat),name='F814W_mag',dtype='float32')   # add a column in which to store magnitudes
    master_cat.add_column(G)
    print("This is the first time you run this program today, isn't it ya ghassan?")
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
#
if (diag_flag == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        print('F160W flag diagnostic\n\na = [flag=0, =1, =2, =3, =4, =-1]\nb=[flux>0,flux<0] for each flag value')
#
## Compute the F160W - F814W offset as the difference between mdeian magnitudes for members which have both flux measurements
#
mag_F160W = [ [], [], [], [], [], [] ]
mag_F814W = [ [], [], [], [], [], [] ]
for counter in range(len(master_cat)):
    for cluster in range(len(limiting_mag)):
        if master_cat['cluster'][counter] == (cluster+1):                 # track by cluster
            if (master_cat['flag_F160W'][counter] == 0 and master_cat['f_F160W'][counter] > 0) and (master_cat['flag_F814W'][counter] == 0 and master_cat['f_F814W'][counter] > 0):
                master_cat['F160W_mag'][counter] = ((-2.5*np.log10(master_cat['f_F160W'][counter]))+25)
                master_cat['F814W_mag'][counter] = ((-2.5*np.log10(master_cat['f_F814W'][counter]))+25)
                mag_F160W[cluster].append(master_cat['F160W_mag'][counter])
                mag_F814W[cluster].append(master_cat['F814W_mag'][counter])
#
## compute medians and offset ( offset = median(F160W_mag) - median(F814W_mag) )

offset = np.array([0.0]*6)
#
for ii in range(len(limiting_mag)):
    temp_F160W = np.array(mag_F160W[ii])
    temp_F814W = np.array(mag_F814W[ii])
    medians = np.array([0.0]*2)
    medians[0] = np.median(temp_F160W)
    medians[1] = np.median(temp_F814W)
#
#    offset[ii] = medians[0] - medians[1]
#
# store the offset in magnitudes to be used between magnitudes measured in F160W filter and F814W
#offset = limiting_mag - limiting_mag_814    # [limiting mag F160] - [limiting mag F814W]
#
#
#
## to store RESULT of program
limiting_mass = np.array([0]*6,dtype='float32')                    
mag_of_limiting_mass = np.array([0]*6,dtype='float32')             # to store magnitude of the object which sets the limiting mass
z_of_limiting_mass =  np.array([0]*6,dtype='float32')             # to store del_z photometric redshift difference ('z_clusterphot') of the object which sets the limiting mass (for comparison w/ redshift cuts)
#
masses_at_lim_mag =  [ [], [], [], [], [], [] ]    # array to store masses at (or near) the limiting magnitude, to determine a range of masses corresponding to the cluster's limiting mag (along w/ mag for a diagnostic check)
#
## MAY NEED TO EDIT: tolerance ("TOL")
TOL = 0.05                                         # this is how near the limiting mag. i'm willing to accept an object
#
mag_by_cluster = [ [], [], [], [], [], [] ]        # two lists: 1-all objects in the ii'th image; 1-cluster members only
mag_by_cluster_member = [ [], [], [], [], [], [] ] # each list (col) contains lists with 3 entries: [mag,mass,del_z(phot)]
# for storing info on stats (diagnostic)
avg_mag_by_cluster = np.array([0.0]*6)               # diagnostic
max_mag_by_cluster = np.array([0.0]*6)               # i.e. store the dimmest (least bright) source
# for storing the max/min masses at the limiting magnitude of each cluster (for plotting)
max_min_mass = np.array([[0]*2]*6,dtype='float32')       # [min,max]
#
mag_spec1 = np.array([0]*6) 
mag_phot1 = np.array([0]*6) 
bad_flag1 = np.array([0]*6)                         # track how many objects are flagged by cluster
bad_flux1 = np.array([0]*6)                         # track # of objects w/ bad fluxes not picked up by 'flag_F160W'
count_nans1 = np.array([0]*6)                       # keep track of objects w/ NaN mass estimates by cluster (entire image)
count_nans_mem = np.array([0]*6) 
count_non_member = np.array([0]*6) 
#
counting_array1 = np.array([0]*9)
#
#
for counter in range(len(master_cat)):             # loop through catalogue
    for cluster in range(len(limiting_mag)):
        if master_cat['cluster'][counter] == (cluster+1):       # start by looking at clusters, one at a time
            if master_cat['member'][counter] == 0:    # now look only at cluster MEMBERS only  
                if master_cat['flag_F160W'][counter] != -1 and master_cat['f_F160W'][counter] > 0:         # look at good objects w/ good fluxes
                    master_cat['F160W_mag'][counter] = ((-2.5*np.log10(master_cat['f_F160W'][counter]))+25)     # AB_mag system zero point = 25, confirmed in Shipley et al. 
                    if (diag_flag_3 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:   # diagnostic check
                        if project_diagnostic_flag == 0:
                            pass
                        else:
                            if np.isnan(master_cat['F160W_mag'][counter]) == 1:      # to identify bad mass estimates (i.e NaNs)
                                print('Diag. 3: NaNs in mag estimates\nNan index in master_cat (mag) = %s'%counter)
                            elif np.isnan(master_cat['lmass'][counter]) == 1:
                                print('Diag. 3: NaNs in mass estimates\nNan index in master_cat (mass) = %s'%counter)
                                #
                    #
                    if np.isnan(master_cat['lmass'][counter]) == 1:          # track all nans in IMAGE
                        count_nans1[cluster]+=1
                        counting_array1[2]+=1 # galaxies w/ bad mass estimates
                    else:
                        #masses[ii].append(master_cat['lmass'][counter])      # store all the masses (regardless of type)
                        mag_by_cluster[cluster].append([master_cat['F160W_mag'][counter],master_cat['lmass'][counter],master_cat['z_clusterphot'][counter]])     # look at all objects in the ii'th image
                    #
                    if np.isnan(master_cat['lmass'][counter]) == 1:          # track all nans in IMAGE
                        count_nans_mem[cluster]+=1
                    elif master_cat['sub'][counter] == 1:     # (spec+phot)
                        mag_spec1[cluster]+=1
                        mag_by_cluster_member[cluster].append([master_cat['F160W_mag'][counter],master_cat['lmass'][counter],master_cat['z_clusterphot'][counter]])
                        if master_cat['F160W_mag'][counter] < (limiting_mag[cluster]+TOL) and master_cat['F160W_mag'][counter] > (limiting_mag[cluster]-TOL):           # if flux is equal to limiting magnitude, +/- 'TOL', store mass value...
                            masses_at_lim_mag[cluster].append([master_cat['F160W_mag'][counter],master_cat['lmass'][counter],master_cat['z_clusterphot'][counter]])
                    elif master_cat['sub'][counter] == 2:   # phot only
                        mag_phot1[cluster]+=1
                        mag_by_cluster_member[cluster].append([master_cat['F160W_mag'][counter],master_cat['lmass'][counter],master_cat['z_clusterphot'][counter]])
                        if master_cat['F160W_mag'][counter] < (limiting_mag[cluster]+TOL) and master_cat['F160W_mag'][counter] > (limiting_mag[cluster]-TOL):           # if flux is equal to limiting magnitude, +/- 'TOL', store mass value...
                            masses_at_lim_mag[cluster].append([master_cat['F160W_mag'][counter],master_cat['lmass'][counter],master_cat['z_clusterphot'][counter]])
                        #
                        #
                        ###  if you want to track SF/Q TYPE galaxies, add a list here
                        #
                        #
                        ## now track objects at the limiting magnitude of the cluster
                        #
                    #
                elif master_cat['flag_F814W'][counter] != -1 and master_cat['f_F814W'][counter] > 0:   #use F814W where F160W unavailable
                    master_cat['F814W_mag'][counter] = ((-2.5*np.log10(master_cat['f_F814W'][counter]))+25)
                    master_cat['F160W_mag'][counter] = master_cat['F814W_mag'][counter] + offset[cluster]    # AB_mag system zero point = 25, confirmed in Shipley et al.; "offset" refers to offset b/w F160W & F814W bands 
                    #
                    if (diag_flag_3 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:   # diagnostic check
                        if project_diagnostic_flag == 0:
                            pass
                        else:
                            if np.isnan(master_cat['F160W_mag'][counter]) == 1:      # to identify bad mass estimates (i.e NaNs)
                                print('Diag. 3: NaNs in mag estimates\nNan index in master_cat (mag) = %s'%counter)
                            elif np.isnan(master_cat['lmass'][counter]) == 1:
                                print('Diag. 3: NaNs in mass estimates\nNan index in master_cat (mass) = %s'%counter)
                                #
                        #
                    #
                    #
                    if np.isnan(master_cat['lmass'][counter]) == 1:          # track all nans in IMAGE
                        count_nans1[cluster]+=1
                        counting_array1[2]+=1 # galaxies w/ bad mass estimates
                    else:
                        mag_by_cluster_member[cluster].append([master_cat['F160W_mag'][counter],master_cat['lmass'][counter],master_cat['z_clusterphot'][counter]])
                        if master_cat['F160W_mag'][counter] < (limiting_mag[cluster]+TOL) and master_cat['F160W_mag'][counter] > (limiting_mag[cluster]-TOL):           # if flux is equal to limiting magnitude, +/- 'TOL', store mass value...
                            masses_at_lim_mag[cluster].append([master_cat['F160W_mag'][counter],master_cat['lmass'][counter],master_cat['z_clusterphot'][counter]])
                        if master_cat['sub'][counter] == 1:     # (spec+phot)
                            mag_spec1[cluster]+=1

                        elif master_cat['sub'][counter] == 2:   # phot only
                            mag_phot1[cluster]+=1
                        #
                        #
                        ###  if you want to track SF/Q TYPE galaxies, add a list here
                        #
                        #
                        ## now track objects at the limiting magnitude of the cluster
                        #
                    #
                else:
                    bad_flux1[cluster]+=1
                #
            else:
                counting_array1[8]+=1   # all non-member galaxies            
                for ii in range(len(limiting_mag)):
                    if master_cat['cluster'][counter] == (ii+1):
                        count_non_member[ii]+=1
            #
        #
    #
    #                         
#
## DIAG_FLAG_4: ACCOUNTING for each galaxy at each condition imposed within the above for loop, which sets up the lists of galaxies to be plotted
#
## SECTION (1.1): SUMMARY TABLE
#
## Summary Table
#
#
if summary_flag == 1 or adams_flag == 1:
    ## Summarize limiting mass calculation in table for PARENT SAMPLE
    lim_mass_names = Column(['Full CATALOGUE ("master_data*.py")','Phot Members','Spec Members','Non-members','Bad flux (F160W<0 & F814W<0 & flag == -1)','Mass = NaN (members)','SUM'],name='Property')
    col_names = cluster_names
    lim_mass0 = Column([np.sum([phot_only,both,spec_only,no_data,stars_sub]),np.sum(mag_phot1),np.sum(mag_spec1),np.sum(count_non_member),np.sum(bad_flux1),np.sum(count_nans_mem),np.sum([mag_phot1,mag_spec1,bad_flux1,count_nans_mem,count_non_member])],name='Total')  # total column
    lim_mass_stats = Table([lim_mass_names,lim_mass0])
    for ii in range(len(mag_phot1)):
        lim_mass_col = Column([np.sum([phot_only[ii],both[ii],spec_only[ii],no_data[ii],stars_sub[ii]]),mag_phot1[ii],mag_spec1[ii],count_non_member[ii],bad_flux1[ii],count_nans_mem[ii],np.sum([mag_phot1[ii],mag_spec1[ii],bad_flux1[ii],count_nans_mem[ii],count_non_member[ii]])],name=col_names[ii])               # cluster columns
        lim_mass_stats.add_column(lim_mass_col) 
    #
    #
    print('\n"data_mass_completeness*.py" Full CATALOGUE breakdown\n%s'%lim_mass_stats)
#
#
#
#
## SECTION (2): determine limiting mass
#
#
for ii in range(len(limiting_mag)):
    mag_by_cluster[ii].sort(key=lambda x: x)                # sort data by magnitude within each clust
    mag_by_cluster_member[ii].sort(key=lambda x: x)
    # temp_*_array(s) store the contents of "mag_by_cluster_*" for a single cluster
    temp_mag_array = np.array([[0]*3]*len(mag_by_cluster[ii]),dtype='float32')        # find max magnitude of image & cluser
    temp_mag_member_array = np.array([[0]*3]*len(mag_by_cluster_member[ii]),dtype='float32')
    # similar to temp*arrays, but dedicated for magnitude
    store_mags = np.array([0]*len(mag_by_cluster[ii]),dtype='float32') 
    store_mags_member = np.array([0]*len(mag_by_cluster_member[ii]),dtype='float32') 
    #
    for jj in range(len(mag_by_cluster[ii])):
        temp_mag_array[jj][0] = mag_by_cluster[ii][jj][0]            # assign fluxes from a given cluster to a dedicated temporary list
        temp_mag_array[jj][1] = mag_by_cluster[ii][jj][1]
        temp_mag_array[jj][2] = mag_by_cluster[ii][jj][2]
        store_mags[jj] = mag_by_cluster[ii][jj][0]                   # magnitude
    for jj in range(len(mag_by_cluster_member[ii])):
        temp_mag_member_array[jj][0] = mag_by_cluster_member[ii][jj][0]  # do the same for cluster members - MAG
        temp_mag_member_array[jj][1] = mag_by_cluster_member[ii][jj][1]  # do the same for cluster members - MASS
        temp_mag_member_array[jj][2] = mag_by_cluster_member[ii][jj][2]  # do the same for cluster members - z_peak (i.e. z_phot)
        store_mags_member[jj] = mag_by_cluster_member[ii][jj][0]                   # magnitude
    #
    if len(masses_at_lim_mag[ii]) == 0:                     # if no objects were found within the tolerance range...
        temp_list = []
        for jj in range(len(temp_mag_member_array)):
            if temp_mag_member_array[jj][0] <= limiting_mag[ii]:
                temp_list.append([temp_mag_member_array[jj][0],temp_mag_member_array[jj][1],temp_mag_member_array[jj][2]])   #[mag,mass,redshift]
        temp_limit = min(temp_list, key=lambda x: abs(x[0]-limiting_mag[ii])) 
        while temp_limit[0] > limiting_mag[ii]:
            print('Closest limiting mass found was above limiting mag. Looking for a new one...')
            index = np.where(temp_mag_member_array == temp_limit[0])
            temp_mag_member_array = np.delete(temp_mag_member_array, (index), axis=0)
            temp_limit = min(temp_mag_member_array, key=lambda x: abs(x[0]-limiting_mag[ii]))
        print('WARNING: NO GALAXIES FOUND AT CLUSTER %s'%(ii+1),' LIMITING MAG. of %s'%limiting_mag[ii],' +/- %.2f'%TOL,'\nClosest thing Old Faithful could find was mass of %.2f'%temp_limit[1],' w/ mag %.2f'%temp_limit[0])
        limiting_mass[ii] = temp_limit[1]  
        mag_of_limiting_mass[ii] = temp_limit[0]
        z_of_limiting_mass[ii] = temp_limit[2]
        max_min_mass[ii][0] = temp_limit[1]               # for PLOTTING
        max_min_mass[ii][1] = temp_limit[1]
    elif len(masses_at_lim_mag[ii]) == 1:
        limiting_mass[ii] = masses_at_lim_mag[ii][0][1]
        mag_of_limiting_mass[ii] = masses_at_lim_mag[ii][0][0]
        z_of_limiting_mass[ii] = masses_at_lim_mag[ii][0][2]
        max_min_mass[ii][0] = masses_at_lim_mag[ii][0][1]               # for PLOTTING
        max_min_mass[ii][1] = masses_at_lim_mag[ii][0][1]
    else:    
        the_good_shit, the_good_shit_index = max_nested_list(masses_at_lim_mag[ii],1)
        the_not_quite_as_good_shit, the_not_quite_as_good_shit_index = min_nested_list(masses_at_lim_mag[ii],1)
        limiting_mass[ii] = the_good_shit[1] # store highest mass you have
        mag_of_limiting_mass[ii] = the_good_shit[0] # store highest mass you have
        z_of_limiting_mass[ii] = the_good_shit[2]
        max_min_mass[ii][0] = the_good_shit[1]               # for PLOTTING
        max_min_mass[ii][1] = the_not_quite_as_good_shit[1]
    if (diag_flag_1 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
        if project_diagnostic_flag == 0:
            pass
        else:
            avg_mag_temp = np.mean(store_mags)                        # report some useful info to the user
            max_mag_temp = max(store_mags)
            avg_mag_temp_member = np.mean(store_mags_member)
            max_mag_temp_member = max(store_mags_member)
            # now print the above
            print('\nDiag. 1: compare w/ "master_data*.py" result\nCluster %s'% (ii+1),': FULL IMAGE STATS\nMax. mag (i.e. dimmest source in ENTIRE image) %.3f'%max_mag_temp,', Avg. mag of IMAGE: %.3f'%avg_mag_temp,'\n# of objects w/ "NaN" mass estimates: %s'%count_nans1[ii])
            print('Cluster %s'%(ii+1),': MEMBER-ONLY STATS\nMax. mag (i.e. dimmest source in CLUSTER only) %.3f'%max_mag_temp_member,', Avg. mag of MEMBERS: %.3f'%avg_mag_temp_member,'\n# of objects w/ "NaN" mass estimates: %s'%count_nans_mem[ii])
            #print('\n"master_data*.py" reports (SF+Q): %s'%((mem[0][ii]+mem[1][ii])),' (spec) cluster members & %s'%((mem_phot[0][ii]+mem_phot[1][ii])),' (phot). Total: %s'%((mem[0][ii]+mem[1][ii])+(mem_phot[0][ii]+mem_phot[1][ii])),' for cluster %i'%(ii+1))
            #print('"master_data*.py" reports total: %s'%((mem[0][ii]+mem_phot[0][ii])),' (SF) cluster members & %s'%((mem[1][ii]+mem_phot[1][ii])),' (Q).')
            print('\nTHIS PROGRAM reports: %i'%(len(mag_by_cluster_member[ii])),' total cluster members for cluster %i'%(ii+1),'\n%i'%bad_flag1[ii],' flagged objects\n%i'%bad_flux1[ii],' objects with bad fluxes (<0)\n%i'%count_nans_mem[ii],' NaNs\nAll this action happened in cluster %i'%(ii+1),', yo.\n')
#        
#
if (diag_flag_2 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        print('\nMembership definition [spec,phot]: %s'%z_cutoff,'\nLimiting magnitudes per Shipley et al. 2018: %s'%limiting_mag)
        print('Limiting masses: %s'%limiting_mass)
        print('Magnitudes at limiting mass: %s'%mag_of_limiting_mass)
#
#
## The code now handles finding the highest mass at the limiting magnitude (or the mass of the object closest in mmagnitude to the limit).
#
#
#
#
## SECTION (3) VISUALIZATION
#
## define range of all histograms to come!
#
max_mass_array = [ [], [], [], [], [], [] ]
max_mass_index = np.array([0]*6)
max_mass = np.array([0]*1)
#
for ii in range(len(z_cluster)):
    max_mass_array[ii], max_mass_index[ii] = max_nested_list(mag_by_cluster_member[ii],1)
max_mass, max_mass_index = max_nested_list(max_mass_array,1)
#
range2=[min(limiting_mass),max_mass[1]]
#
## Now I have lists of corresponding magnitudes and masses for each cluster. make a single plot for the first cluster (as a clear example of the method), then a subplot (w/ 6 plots) all together. add lines to the plot for (A) the limiting magnitude (vertical); and (B) range of masses at that magnitude (horizontal), for clusters w/ many such objects, or circle/point to/identify in some way on the plot the galaxy closest to the limiting magnitude which was chosen as the standard for limiting mass
#
## Fig 1: mass vs magnitude scatter plot, colourbar = redshift (z_phot)
## Now make the single scatter plot for M0416 (i.e. cluster stored in the first position of the above lists)
#
## 
## setup dedicated arrays for plotting, since plotting them from my nested list of lists of lists doesn't seem to be working...
#
if (plot_flag_1 == 1 and project_plot_flag ==2) or project_plot_flag == 1:
    if project_plot_flag == 0:
        pass
    else:
        for cluster in range(len(limiting_mag)):
            plotting_array_temp = np.array([[0]*len(mag_by_cluster_member[cluster])]*3,dtype='float32') # 0=mag;1=mass;2=z
            for ii in range(len(mag_by_cluster_member[cluster])):
                plotting_array_temp[0][ii] = mag_by_cluster_member[cluster][ii][0]   # store mass/magnitudes for the first cluster
                plotting_array_temp[1][ii] = mag_by_cluster_member[cluster][ii][1]   # index order: [cluster][object][mag/mass,z_phot]([cluster,row,col])
                plotting_array_temp[2][ii] = mag_by_cluster_member[cluster][ii][2]
            #
            ## calculate cluster members; display as "Total(Shown)"
            total_mem = np.sum([mem_phot[0][cluster],mem_phot[1][cluster],mem_spec[0][cluster],mem_spec[1][cluster]])
            mem_shown = len(plotting_array_temp[0])
            #
            members_string = '# members(F160W+F814W): %s'%total_mem+'(%s'%mem_shown+')'
            ## Compute std. dev. of cluster objects
            std_dev = np.std(plotting_array_temp[1])
            #
            ## The figure:
            #
            c_map = plt.get_cmap("spring")        # set colorbar map
            #
            #    plt.close()
            fig, ax = plt.subplots()
            ax.scatter(plotting_array_temp[0],plotting_array_temp[1], marker='o', s=10, c=plotting_array_temp[2], cmap= c_map, vmin=min(plotting_array_temp[2]), vmax=max(plotting_array_temp[2]))#, edgecolors='k')#,facecolors='none', edgecolors='k')
            #ax.scatter(SF_masses[0][0],SF_masses[0][1], marker='*', s=10, facecolors='none', edgecolors='b')
            ax.plot([limiting_mag[cluster],limiting_mag[cluster]],[0,15],'--r', linewidth=1.4)
            ax.plot([0,35],[max_min_mass[cluster][0],max_min_mass[cluster][0]], color='maroon', linestyle='--',linewidth=1.2)    #chang indexing of max_min_mass to [ii][0] & [ii][1] for looped subplots
            ax.plot([0,35],[max_min_mass[cluster][1],max_min_mass[cluster][1]], color='maroon', linestyle='--',linewidth=1.2)
            sm =  ScalarMappable(cmap=c_map)
            sm.set_array([])
            sm.set_clim([0,z_cutoff[1]])
            cbar = fig.colorbar(sm, ax=ax)
            cbar.ax.set_title("|${\Delta}$z|$_{phot}$")
            #ax.pcolormesh(plotting_array_temp[0], plotting_array_temp[1], plotting_array_temp[2], cmap=c_map)
            #ax.colorbar()

            #plt.plot([0,35],[max_mass,max_mass], '-.k', linewidth=0.8)
            #plt.plot([0,35],[min_mass,min_mass], '-.k', linewidth=0.8)
            ax.set_xlabel('$m_{F160W}$')
            ax.set_xlim(17,30)
            ax.set_ylabel('$log_{10}$(M/M$_{\odot})$')
            ax.set_ylim(5,13)
            ax.grid(b=False, which='major', axis='both', color = 'k', linestyle = ':')
            ax.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=False, labelleft=True)
            ax.minorticks_on()
            #
            plt.text(20.1,12.1,cluster_names[cluster],fontsize=12)
            plt.text(20.1,11.5,'$M/L_{max}$: %s'%limiting_mass[cluster],fontsize=10)
            plt.text(20.1,10.5,members_string,fontsize=10)
            plt.text(18.1,6.1,'z = %s'%z_cluster[cluster],fontsize=10)
            plt.text(18.1,5.7,'$\sigma_{mass}$ = %s'%np.round(std_dev,decimals=3),fontsize=10)
            #
            plt.show()
            #
            #
            lim_mag_minus = (limiting_mag[cluster]-TOL)
            lim_mag_plus = (limiting_mag[cluster]+TOL)
        #
        fig, axs = plt.subplots(nrows=2,ncols=3,sharex=True,sharey=True)#,tight_layout=True)
        fig.subplots_adjust(wspace=0,hspace=0)
        #
        for cluster in range(len(limiting_mag)):
            plotting_array_temp = np.array([[0]*len(mag_by_cluster_member[cluster])]*3,dtype='float32') # 0=mag;1=mass;2=z
            for ii in range(len(mag_by_cluster_member[cluster])):
                plotting_array_temp[0][ii] = mag_by_cluster_member[cluster][ii][0]   # store mass/magnitudes for the first cluster
                plotting_array_temp[1][ii] = mag_by_cluster_member[cluster][ii][1]   # index order: [cluster][object][mag/mass,z_phot]([cluster,row,col])
                plotting_array_temp[2][ii] = mag_by_cluster_member[cluster][ii][2]
            #
            ## calculate cluster members; display as "Total(Shown)"
            total_mem = np.sum([mem_phot[0][cluster],mem_phot[1][cluster],mem_spec[0][cluster],mem_spec[1][cluster]])
            mem_shown = len(plotting_array_temp[0])
            #
            members_string = '# cluster members(F160W+F814W): %s'%total_mem+'(%s'%mem_shown+')'
            ## Compute std. dev. of cluster objects
            std_dev = np.std(plotting_array_temp[1])
            #
            ## The figure:
            #
            c_map = plt.get_cmap("spring")        # set colorbar map
            #
            lim_mag_minus = (limiting_mag[cluster]-TOL)
            lim_mag_plus = (limiting_mag[cluster]+TOL)
            #
            #
            ax = axs.flat[cluster]
            #
            ax.scatter(plotting_array_temp[0],plotting_array_temp[1], marker='o', s=30, c=plotting_array_temp[2], cmap= c_map, vmin=min(plotting_array_temp[2]), vmax=max(plotting_array_temp[2]))
            ax.plot([limiting_mag[cluster],limiting_mag[cluster]],[0,15],'--k', linewidth=1.4)
            #ax.plot([lim_mag_minus,lim_mag_minus],[0,15],':r', linewidth=1.2)
            #ax.plot([lim_mag_plus,lim_mag_plus],[0,15],':r', linewidth=1.2)
            ax.plot([0,35],[max_min_mass[cluster][0],max_min_mass[cluster][0]], color='maroon', linestyle='--',linewidth=1.2)    #chang indexing of max_min_mass to [ii][0] & [ii][1] for looped subplots
            ax.plot([0,35],[max_min_mass[cluster][1],max_min_mass[cluster][1]], color='maroon', linestyle='--',linewidth=1.2)
            sm =  ScalarMappable(cmap=c_map)
            sm.set_array([])
            sm.set_clim([0,z_cutoff[1]])
            ax.grid(axis='both', alpha=0.75)
            ax.set_xlim(17,30)
            ax.set_ylim(5,13)
            ## label locations
            axs.flat[cluster].text(20.1,12.1,cluster_names[cluster],fontsize=12)
            axs.flat[cluster].text(20.1,11.5,'$M/L_{max}$: %s'%limiting_mass[cluster],fontsize=12)
            axs.flat[cluster].text(20.1,10.5,members_string,fontsize=10)
            axs.flat[cluster].text(18.1,6.1,'z = %s'%z_cluster[cluster],fontsize=10)
            axs.flat[cluster].text(18.1,5.7,'$\sigma_{mass}$ = %s'%std_dev,fontsize=10)
            #
            #plt.text(21.6,11.5,'$M/L_{max}$: %s'%limiting_mass[cluster],fontsize=10)
        #
        cbar = fig.colorbar(sm, ax=axs[:, 2],location='right')#, vmin=0.0, vmax=z_cutoff[1])
        cbar.ax.set_title("|${\Delta}$z|$_{phot}$")
        #            
        #            
        plt.show()
        #
        #
    #
#
#
## Visualize by cluster
#
if (plot_flag_2 == 1 and project_plot_flag ==2) or project_plot_flag == 1:
    if project_plot_flag == 0:
        pass
    else:
        #
        fig, axs = plt.subplots(nrows=2,ncols=3,sharex=True,sharey=True)#,tight_layout=True)
        fig.subplots_adjust(wspace=0,hspace=0)
        #
        for cluster in range(len(limiting_mag)):
            plotting_array_temp = np.array([[0]*len(mag_by_cluster_member[cluster])]*3,dtype='float32') # 0=mag;1=mass;2=z
            for ii in range(len(mag_by_cluster_member[cluster])):
                plotting_array_temp[0][ii] = mag_by_cluster_member[cluster][ii][0]   # store mass/magnitudes for the first cluster
                plotting_array_temp[1][ii] = mag_by_cluster_member[cluster][ii][1]   # index order: [cluster][object][mag/mass,z_phot]([cluster,row,col])
                plotting_array_temp[2][ii] = mag_by_cluster_member[cluster][ii][2]
            #
            ## calculate cluster members; display as "Total(Shown)"
            total_mem = np.sum([mem_phot[0][cluster],mem_phot[1][cluster],mem_spec[0][cluster],mem_spec[1][cluster]])
            mem_shown = len(plotting_array_temp[0])
            #
            members_string = '# cluster members(F160W+F814W): %s'%total_mem+'(%s'%mem_shown+')'
            ## Compute std. dev. of cluster objects
            std_dev = np.std(plotting_array_temp[1])
            #
            ## The figure:
            #
            c_map = plt.get_cmap("spring")        # set colorbar map
            #
            lim_mag_minus = (limiting_mag[cluster]-TOL)
            lim_mag_plus = (limiting_mag[cluster]+TOL)
            #
            #
            ax = axs.flat[cluster]
            #
            ax.scatter(plotting_array_temp[0],plotting_array_temp[1], marker='o', s=30, c=plotting_array_temp[2], cmap= c_map, vmin=min(plotting_array_temp[2]), vmax=max(plotting_array_temp[2]))
            ax.plot([limiting_mag[cluster],limiting_mag[cluster]],[0,15],'--k', linewidth=1.4)
            ax.plot([lim_mag_minus,lim_mag_minus],[0,15],':r', linewidth=1.2)
            ax.plot([lim_mag_plus,lim_mag_plus],[0,15],':r', linewidth=1.2)
            ax.plot([0,35],[max_min_mass[cluster][0],max_min_mass[cluster][0]], color='maroon', linestyle='--',linewidth=1.2)    #chang indexing of max_min_mass to [ii][0] & [ii][1] for looped subplots
            ax.plot([0,35],[max_min_mass[cluster][1],max_min_mass[cluster][1]], color='maroon', linestyle='--',linewidth=1.2)
            sm =  ScalarMappable(cmap=c_map)
            sm.set_array([])
            sm.set_clim([0,z_cutoff[1]])
            ax.grid(axis='both', alpha=0.75)
            ax.set_xlim([25,27])
            ax.set_ylim([5.5,8.5])
            #
            #plt.text(21.6,11.5,'$M/L_{max}$: %s'%limiting_mass[cluster],fontsize=10)
        ## label locations
        axs.flat[0].text(25.15,6.05,cluster_names[0],fontsize=12)
        axs.flat[0].text(25.15,5.75,'$M/L_{max}$: %s'%limiting_mass[0],fontsize=12)
        axs.flat[1].text(25.15,6.05,cluster_names[1],fontsize=12)
        axs.flat[1].text(25.15,5.75,'$M/L_{max}$: %s'%limiting_mass[1],fontsize=12)
        axs.flat[2].text(25.15,6.05,cluster_names[2],fontsize=12)
        axs.flat[2].text(25.15,5.75,'$M/L_{max}$: %s'%limiting_mass[2],fontsize=12)
        axs.flat[3].text(25.15,6.05,cluster_names[3],fontsize=12)
        axs.flat[3].text(25.15,5.75,'$M/L_{max}$: %s'%limiting_mass[3],fontsize=12)
        axs.flat[4].text(25.15,6.05,cluster_names[4],fontsize=12)
        axs.flat[4].text(25.15,5.75,'$M/L_{max}$: %s'%limiting_mass[4],fontsize=12)
        axs.flat[5].text(25.15,6.05,cluster_names[5],fontsize=12)
        axs.flat[5].text(25.15,5.75,'$M/L_{max}$: %s'%limiting_mass[5],fontsize=12)
        #axs.clim(0,z_cutoff[1])
        cbar = fig.colorbar(sm, ax=axs[:, 2],location='right')#, vmin=0.0, vmax=z_cutoff[1])
        cbar.ax.set_title("|${\Delta}$z|$_{phot}$")
        #            
        #            
        plt.show()
    #
#
#
if (plot_flag_3 == 1 and project_plot_flag ==2) or project_plot_flag == 1:
    if project_plot_flag == 0:
        pass
    else:
        for cluster in range(len(limiting_mag)): 
            plotting_array_temp = np.array([[0]*len(mag_by_cluster_member[cluster])]*3,dtype='float32') # 0=mag;1=mass;2=z
            for ii in range(len(mag_by_cluster_member[cluster])):
                plotting_array_temp[0][ii] = mag_by_cluster_member[cluster][ii][0]   # store mass/magnitudes for the first cluster
                plotting_array_temp[1][ii] = mag_by_cluster_member[cluster][ii][1]   # index order: [cluster (ii)][object][mag/mass,z_phot]([cluster,row,col])
                plotting_array_temp[2][ii] = mag_by_cluster_member[cluster][ii][2]
            #
            ## calculate cluster members; display as "Total(Shown)"
            total_mem = np.sum([mem_phot[0][cluster],mem_phot[1][cluster],mem_spec[0][cluster],mem_spec[1][cluster]])
            mem_shown = len(plotting_array_temp[0])
            #
            members_string = '$M^{*}_{lim}$: %s'%limiting_mass[cluster]+'\n# cluster members(F160W): %s'%total_mem+'(%s'%mem_shown+')\nz$_{cluster}$ = %s'%z_cluster[cluster]
            ## Compute std. dev. of cluster objects
            std_dev = np.std(plotting_array_temp[1])
            #
            ## The figure:
            #
            fig, axs = plt.subplots(nrows=1,ncols=2,sharex=False,sharey=False)#,tight_layout=True)
            fig.subplots_adjust(wspace=0,hspace=0)
            fig.set_size_inches(20,8)
            #
            #
            c_map = plt.get_cmap("autumn")        # set colorbar map
            #
            lim_mag_minus = (limiting_mag[cluster]-TOL)
            lim_mag_plus = (limiting_mag[cluster]+TOL)
            #
            #
            ax1 = axs.flat[0]
            #
            ax1.scatter(plotting_array_temp[0],plotting_array_temp[1], marker='o', s=150, c=plotting_array_temp[2], cmap= c_map, vmin=min(plotting_array_temp[2]), vmax=max(plotting_array_temp[2]))
            ax1.plot([limiting_mag[cluster],limiting_mag[cluster]],[0,15],'--k', linewidth=2.5)
            ax1.plot([lim_mag_minus,lim_mag_minus],[0,15],':r', linewidth=1.5)
            ax1.plot([lim_mag_plus,lim_mag_plus],[0,15],':r', linewidth=1.5)
            ax1.plot([0,35],[max_min_mass[cluster][0],max_min_mass[cluster][0]], color='k', linestyle='--',linewidth=2.5)    #chang indexing of max_min_mass to [ii][0] & [ii][1] for looped subplots
            ax1.plot([0,35],[max_min_mass[cluster][1],max_min_mass[cluster][1]], color='k', linestyle='--',linewidth=2.5)
            #
            ax1.grid(axis='both', alpha=0.75)
            ax1.set_xlim([17,30])
            ax1.set_xlabel('m$_{F160W}$',fontsize=20)
            ax1.set_ylim([5,13])
            ax1.set_ylabel('$log_{10}$(M/M$_{\odot})$',fontsize=20)
            ax1.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=False,labelright=False, labelleft=True,labelsize=20)
            #
            ## label locations
            ax1.text(20.1,12.1,cluster_names[cluster],fontsize=25)
            ax1.text(20.1,11.0,members_string,fontsize=20)
            #
            #
            # zoom in
            ax2 = axs.flat[1]
            #
            #
            ## define function for drawing a circle on plots
            def circle(x, y, radius=0.15):
                from matplotlib.patches import Circle
                from matplotlib.patheffects import withStroke
                circle = Circle((x, y), radius, clip_on=False, zorder=10, linewidth=1,
                                edgecolor='black', facecolor=(0, 0, 0, .0125),
                                path_effects=[withStroke(linewidth=5, foreground='w')])
                ax2.add_artist(circle)
            #
            #
            ax2.scatter(plotting_array_temp[0],plotting_array_temp[1], marker='o', s=150, c=plotting_array_temp[2], cmap= c_map, vmin=min(plotting_array_temp[2]), vmax=max(plotting_array_temp[2]))
            ax2.plot([limiting_mag[cluster],limiting_mag[cluster]],[0,15],'--k', linewidth=2.0)
            ax2.plot([lim_mag_minus,lim_mag_minus],[0,15],':r', linewidth=2.0)
            ax2.plot([lim_mag_plus,lim_mag_plus],[0,15],':r', linewidth=2.0)
            ax2.plot([0,35],[max_min_mass[cluster][0],max_min_mass[cluster][0]], color='k', linestyle='--',linewidth=2.5)    #chang indexing of max_min_mass to [ii][0] & [ii][1] for looped subplots
            ax2.plot([0,35],[max_min_mass[cluster][1],max_min_mass[cluster][1]], color='k', linestyle='--',linewidth=2.5)
            sm =  ScalarMappable(cmap=c_map)
            sm.set_array([])
            sm.set_clim([0,z_cutoff[1]])
            ax2.grid(axis='both', alpha=0.75)
            ax2.set_xlim([25.5,27.5])
            ax2.set_xlabel('m$_{F160W}$',fontsize=20)
            ax2.set_ylim([(max_min_mass[cluster][1]-0.2),9])
            circle(mag_of_limiting_mass[cluster],limiting_mass[cluster])
            ax2.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=True, labelleft=False,labelsize=20)
            cbar = fig.colorbar(sm, ax=ax2)#,location='right')#, vmin=0.0, vmax=z_cutoff[1])
            cbar.ax.set_title("|${\Delta}$z|$_{phot}$",fontsize=20)
            #     
            plt.show()

#
#
#
## SECTION (4): REMOVE all galaxies from cluster member sample which are below limiting mass (filter: from member=0(1) to member=6(7) for cluster(field) galaxies)
#
##
above_lim_mass_mem = np.array([0]*6)
above_lim_mass_field = np.array([0]*6)
below_lim_mass_mem = np.array([0]*6)
below_lim_mass_field = np.array([0]*6)
other = np.array([0]*6)
#
for counter in range(len(master_cat)):
    for cluster in range(len(below_lim_mass_mem)):
        if master_cat['cluster'][counter] == (cluster+1):              # track galaxies by CLUSTER
            if master_cat['member'][counter] == 0:                     # deal with CLUSTER MEMBERS
                if master_cat['lmass'][counter] < limiting_mass[cluster]:    # isolate galaxies below LIMITING MASS by cluster
                    master_cat['member'][counter] = 6                  # RECLASSIFY
                    below_lim_mass_mem[cluster]+=1                     # keep count of numbers
                else:
                    above_lim_mass_mem[cluster]+=1
            elif master_cat['member'][counter] == 1:                     # deal with FIELD MEMBERS
                if master_cat['lmass'][counter] < limiting_mass[cluster]:    # isolate galaxies below LIMITING MASS by cluster
                    master_cat['member'][counter] = 7                  # RECLASSIFY
                    below_lim_mass_field[cluster]+=1                     # keep count of numbers
                else:
                    above_lim_mass_field[cluster]+=1
            else: 
                other[cluster]+=1                                      # all others (i.e. not cluster members or field members)
        #
    #
#
#
#
## SECTION (4.1): Summary Table for Cluster/Field MEMBERS
#
#
if summary_flag == 1 or adams_flag == 1:
    ## Summarize limiting mass calculation in table for PARENT SAMPLE
    below_lim_mass_names = Column(['CLUSTER members ("master_data*.py")','FIELD members ("master_data*.py")','Cluster above lim. mass','Cluster below lim. mass','SUM (Cluster)','Field above lim. mass','Field below lim. mass','SUM (Field)'],name='Property')
    col_names = cluster_names
    below_lim_mass0 = Column([np.sum([mem_phot,mem_spec]),np.sum([field_phot,field_spec]),np.sum(above_lim_mass_mem),np.sum(below_lim_mass_mem),np.sum([above_lim_mass_mem,below_lim_mass_mem]),np.sum(above_lim_mass_field),np.sum(below_lim_mass_field),np.sum([above_lim_mass_field,below_lim_mass_field])],name='Total')  # total column
    below_lim_mass_stats = Table([below_lim_mass_names,below_lim_mass0])
    for ii in range(len(mag_phot1)):
        below_lim_mass_col = Column([np.sum([mem_phot[0][ii],mem_phot[1][ii],mem_spec[0][ii],mem_spec[1][ii]]),np.sum([field_phot[0][ii],field_phot[1][ii],field_spec[0][ii],field_spec[1][ii]]),above_lim_mass_mem[ii],below_lim_mass_mem[ii],np.sum([above_lim_mass_mem[ii],below_lim_mass_mem[ii]]),above_lim_mass_field[ii],below_lim_mass_field[ii],np.sum([above_lim_mass_field[ii],below_lim_mass_field[ii]])],name=col_names[ii])               # cluster columns
        below_lim_mass_stats.add_column(below_lim_mass_col) 
    #
    #
    #
#
    print('\n"data_mass_completeness*.py" Full CATALOGUE breakdown - POST limiting mass calculation\n%s'%below_lim_mass_stats)
#
#
#
## TIME_FLAG END
#
if time_flag == 1:
    if project_time_flag == 1:
        pass
    else:
        print('Program "data_mass_completeness*.py" took: %s seconds to run.\n\n' % (time.time() - start_time))
#
#
#
#
print('\n\n"data_mass_completeness_F160W.py"  terminated successfully.\n')
#
#                        
###################     PROGRAM END     ###################