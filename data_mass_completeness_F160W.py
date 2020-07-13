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
#
limiting_mag = np.array([26.2,26.5,25.5,26.1,26.5,26.7])    # NOTE: order is [M0416,M1149,M0717,A370,A1063,A274]
cluster_names = ['M0416','M1149','M0717','A370','A1063','A2744']   # in the order corresponding to limiting mags above
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
    F = Column([-99]*len(master_cat),name='F160W_mag',dtype='float32')   # add a column in which to store magnitudes
    master_cat.add_column(F)
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
print('F160W flag diagnostic\n\na = [flag=0, =1, =2, =3, =4, =-1]\nb=[flux>0,flux<0] for each flag value')
#
a = np.array([0]*6)
b = np.array([[0]*2]*6)
for counter in range(len(master_cat)):
    for ii in range(-1,5):
        if master_cat['flag_F160W'][counter] == ii:
            a[ii]+=1
            if master_cat['f_F160W'][counter] > 0:
                b[ii][0]+=1
            else: b[ii][1]+=1
print('a = %s'%a)
print('total count: %s'%np.sum(a))
print('length of catalogue: %s'%len(master_cat))
ratios = []
for ii in range(len(a)):
    ratios.append(a[ii]/len(master_cat))
print(ratios)
print(np.sum(ratios))
print('b = %s'%b)
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
count_non_member = np.array([0]*6) 
#
counting_array1 = np.array([0]*9)
#
#
for counter in range(len(master_cat)):             # loop through catalogue
    if master_cat['flag_F160W'][counter] == 0:         # avoid objects w/ bad phot; good galaxies marked by flag_band == 0
        counting_array1[0]+=1  # good flag galaxies
        if master_cat['f_F160W'][counter] > 0:         # double check, since bad phots are set to == -99
            counting_array1[1]+=1  # good flux galaxies
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
            for ii in range(len(limiting_mag)):
                if master_cat['cluster'][counter] == (ii+1):                 # track objects by cluster
                    if np.isnan(master_cat['lmass'][counter]) == 1:          # track all nans in IMAGE
                        count_nans1[ii]+=1
                        counting_array1[2]+=1 # galaxies w/ bad mass estimates
                    else:
                        #masses[ii].append(master_cat['lmass'][counter])      # store all the masses (regardless of type)
                        mag_by_cluster[ii].append([master_cat['F160W_mag'][counter],master_cat['lmass'][counter],master_cat['z_clusterphot'][counter]])     # look at all objects in the ii'th image
                        if np.isnan(master_cat['lmass'][counter]) == 1:      # track all nans in CLUSTER
                            counting_array1[3]+=1
                        else:
                            if master_cat['member'][counter] == 0:    # now look only at cluster MEMBERS only  
                                counting_array1[4]+=1
                                if master_cat['sub'][counter] == 1:     # (spec+phot)
                                    mag_spec1[ii]+=1
                                elif master_cat['sub'][counter] == 2:   # phot only
                                    mag_phot1[ii]+=1
                                #
                                #
                                ###  if you want to track SF/Q TYPE galaxies, add a list here
                                #
                                #
                                ## now track objects at the limiting magnitude of the cluster
                                #
                            else:
                                counting_array1[8]+=1   # all non-member galaxies            
                                for ii in range(len(limiting_mag)):
                                    if master_cat['cluster'][counter] == (ii+1):
                                        count_non_member[ii]+=1
                            #  
        # 
        else:
            for cluster in range(len(limiting_mag)):
                if master_cat['cluster'][counter] == (cluster+1):
                    bad_flux1[cluster]+=1
                    counting_array1[5]+=1    # bad flux (<0)
    else:
        counting_array1[6]+=1     # bad flag ALL
        for cluster in range(len(limiting_mag)):
            if master_cat['cluster'][counter] == (cluster+1):
                bad_flag1[cluster]+=1
        for cluster in range(len(limiting_mag)):
            if master_cat['cluster'][counter] == (cluster+1):
                if master_cat['member'][counter] == 0:     # we're only interesting in comparing for cluster members
                    
                    counting_array1[7]+=1     # bad flag member
            #
#
if (diag_flag_4 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    print('\n[good_flag,good_flux,NaN-mass,NaN-mass(member),member,bad_flux,bad_flag,bad_flag_member,non_member]:\n%s'%counting_array1)
#
#
mag_spec2 = np.array([0]*6) 
mag_phot2 = np.array([0]*6) 
bad_flag2 = np.array([0]*6) 
bad_flux2 = np.array([0]*6) 
count_nans2 = np.array([0]*6)                       # keep track of objects w/ NaN mass estimates by cluster (entire image)
counting_array2 = np.array([0]*6)
#
for counter in range(len(master_cat)):
    if master_cat['member'][counter] == 0:
        counting_array2[0]+=1                                     # number of members
        if master_cat['flag_F160W'][counter] == 0:                # good flag
            counting_array2[1]+=1                       
            if master_cat['f_F160W'][counter] > 0:                # good flux
                counting_array2[2]+=1
                if np.isnan(master_cat['lmass'][counter]) == 1:
                    for ii in range(len(limiting_mag)):
                        if master_cat['cluster'][counter] == (ii+1):
                            count_nans2[ii]+=1
                elif master_cat['sub'][counter] == 1:     # (spec+phot)
                    for ii in range(len(limiting_mag)):
                        if master_cat['cluster'][counter] == (ii+1):
                            mag_by_cluster_member[ii].append([master_cat['F160W_mag'][counter],master_cat['lmass'][counter],master_cat['z_clusterphot'][counter]])
                            mag_spec2[ii]+=1
                            if master_cat['F160W_mag'][counter] < (limiting_mag[ii]+TOL) and master_cat['F160W_mag'][counter] > (limiting_mag[ii]-TOL):           # if flux is equal to limiting magnitude, +/- 'TOL', store mass value...
                                masses_at_lim_mag[ii].append([master_cat['F160W_mag'][counter],master_cat['lmass'][counter],master_cat['z_clusterphot'][counter]])
                                    #
                elif master_cat['sub'][counter] == 2:   # phot only
                    for ii in range(len(limiting_mag)):
                        if master_cat['cluster'][counter] == (ii+1):
                            mag_by_cluster_member[ii].append([master_cat['F160W_mag'][counter],master_cat['lmass'][counter],master_cat['z_clusterphot'][counter]])
                            mag_phot2[ii]+=1
                            if master_cat['F160W_mag'][counter] < (limiting_mag[ii]+TOL) and master_cat['F160W_mag'][counter] > (limiting_mag[ii]-TOL):           # if flux is equal to limiting magnitude, +/- 'TOL', store mass value...
                                masses_at_lim_mag[ii].append([master_cat['F160W_mag'][counter],master_cat['lmass'][counter],master_cat['z_clusterphot'][counter]])
                                #
                        #
                else:
                    counting_array2[3]+=1               # outside parent sample
            else:
                for ii in range(len(limiting_mag)):
                    if master_cat['cluster'][counter] == (ii+1):  # bad flux
                        bad_flux2[ii]+=1
                counting_array2[4]+=1                            
        else:
            for ii in range(len(limiting_mag)):
                if master_cat['cluster'][counter] == (ii+1):      # bad flag
                    bad_flag2[ii]+=1
            counting_array2[5]+=1                           
            
#
if (diag_flag_4 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    print('\n[members,good_flag,good_flux,not_in_parent,bad_flux,bad_flag]:\n%s'%counting_array2)
#
## DIAG_FLAG_4: ACCOUNTING for each galaxy at each condition imposed within the above for loop, which sets up the lists of galaxies to be plotted
if (diag_flag_4 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    print('\n"flag_F160W==0" OK: %s'%counting_array1[0],'\n"flag_F160W!=0" BAD: %s'%counting_array1[6],'\nGood flag + Bad flag = %s'%counting_array1[0],' + %s'%counting_array1[6],' = %s'%np.sum([counting_array1[0],counting_array1[6]]))
    print('\nOf the "Good flag" galaxies\nflux_F160W>0: %s'%counting_array1[1],'\nflux_F160W<0: %s'%counting_array1[5],'\n(flux>0) + (flux<0) = %s'%counting_array1[1],' + %s'%counting_array1[5],' = %s'%(np.sum([counting_array1[1],counting_array1[5]])))
    print('\nlmass = NaN (all): %s'%counting_array1[2],'\nCluster members selected: %s'%counting_array1[4],'\nlmass = NaN (members): %s'%counting_array1[3],'\nBad flag: %s'%counting_array1[7],'\nBad flux (<0): %s'%counting_array1[5])

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
        print('WARNING: NO GALAXIES FOUND AT CLUSTER %s'%(ii+1),' LIMITING MAG. of %s'%limiting_mag[ii],' +/- %.2f'%TOL,'\nClosest thing Old Faithful could find was %.2f'%temp_limit[1],' w/ mag %.2f'%temp_limit[0])
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
            print('\nDiag. 1: compare w/ "master_data*.py" result\nCluster %s'% (ii+1),': FULL IMAGE STATS\nMax. mag (i.e. dimmest source in ENTIRE image) %.3f'%max_mag_temp,', Avg. mag of IMAGE: %.3f'%avg_mag_temp,'\n# of objects w/ "NaN" mass estimates: %s'%count_nans[ii])
            print('Cluster %s'%(ii+1),': MEMBER-ONLY STATS\nMax. mag (i.e. dimmest source in CLUSTER only) %.3f'%max_mag_temp_member,', Avg. mag of MEMBERS: %.3f'%avg_mag_temp_member,'\n# of objects w/ "NaN" mass estimates: %s'%count_nans_member[ii])
            #print('\n"master_data*.py" reports (SF+Q): %s'%((mem[0][ii]+mem[1][ii])),' (spec) cluster members & %s'%((mem_phot[0][ii]+mem_phot[1][ii])),' (phot). Total: %s'%((mem[0][ii]+mem[1][ii])+(mem_phot[0][ii]+mem_phot[1][ii])),' for cluster %i'%(ii+1))
            #print('"master_data*.py" reports total: %s'%((mem[0][ii]+mem_phot[0][ii])),' (SF) cluster members & %s'%((mem[1][ii]+mem_phot[1][ii])),' (Q).')
            print('\nTHIS PROGRAM reports: %i'%(len(mag_by_cluster_member[ii])),' total cluster members for cluster %i'%(ii+1),'\n%i'%bad_flag[ii],' flagged objects\n%i'%bad_flux[ii],' objects with bad fluxes (<0)\n%i'%count_nans_member[ii],' NaNs\nAll this action happened in cluster %i'%(ii+1),', yo.\n')
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
## SECTION (1.1): SUMMARY TABLE
#
## Summary Table
#
#
if summary_flag == 1 or adams_flag == 1:
    ## Summarize limiting mass calculation in table for PARENT SAMPLE
    lim_mass_names = Column(['Full CATALOGUE ("master_data*.py")','Phot Members','Spec Members','Non-members','Bad flag (!=0)','Bad flux (F160W<0)','Mass = NaN','SUM'],name='Property')
    col_names = cluster_names
    lim_mass0 = Column([np.sum([phot_only,both,spec_only,no_data,stars_sub]),np.sum(mag_phot1),np.sum(mag_spec1),np.sum(count_non_member),np.sum(bad_flag1),np.sum(bad_flux1),np.sum(count_nans1),np.sum([mag_phot1,mag_spec1,bad_flag1,bad_flux1,count_nans1,count_non_member])],name='Total')  # total column
    lim_mass_stats = Table([lim_mass_names,lim_mass0])
    for ii in range(len(mag_phot1)):
        lim_mass_col = Column([np.sum([phot_only[ii],both[ii],spec_only[ii],no_data[ii],stars_sub[ii]]),mag_phot1[ii],mag_spec1[ii],count_non_member[ii],bad_flag1[ii],bad_flux1[ii],count_nans1[ii],np.sum([mag_phot1[ii],mag_spec1[ii],bad_flag1[ii],bad_flux1[ii],count_nans1[ii],count_non_member[ii]])],name=col_names[ii])               # cluster columns
        lim_mass_stats.add_column(lim_mass_col) 
    #
    #
    #
    ## Summarize limiting mass calculation in table for CLUSTER MEMBERS
    lim_mass_member_names = Column(['Members ("master_data*.py")','Phot Members','Spec Members','Bad flag (!=0)','Bad flux (F160W<0)','Mass = NaN','SUM'],name='Property')
    col_names = cluster_names
    lim_member_mass0 = Column([np.sum([mem_phot,mem_spec]),np.sum(mag_phot2),np.sum(mag_spec2),np.sum(bad_flag2),np.sum(bad_flux2),np.sum(count_nans2),np.sum([mag_phot2,mag_spec2,bad_flag2,bad_flux2,count_nans2])],name='Total')  # total column
    lim_mass_member_stats = Table([lim_mass_member_names,lim_member_mass0])
    for ii in range(len(mag_phot1)):
        lim_mass_member_col = Column([np.sum([mem_phot[0][ii],mem_phot[1][ii],mem_spec[0][ii],mem_spec[1][ii]]),mag_phot2[ii],mag_spec2[ii],bad_flag2[ii],bad_flux2[ii],count_nans2[ii],np.sum([mag_phot2[ii],mag_spec2[ii],bad_flag2[ii],bad_flux2[ii],count_nans2[ii]])],name=col_names[ii])               # cluster columns
        lim_mass_member_stats.add_column(lim_mass_member_col) 
    #
#
    print('\n"data_mass_completeness*.py" Full CATALOGUE breakdown\n%s'%lim_mass_stats)
    print('\nMember sample breakdown\n%s'%lim_mass_member_stats)
#
#
#
#
## SECTION (2) VISUALIZATION
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
if plot_flag_1 == 1:
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
        members_string = '# members(F160W): %s'%total_mem+'(%s'%mem_shown+')'
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
        ax.set_ylabel('$log(M/M_{\odot})$')
        ax.set_ylim(5,13)
        ax.grid(b=False, which='major', axis='both', color = 'k', linestyle = ':')
        ax.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=False, labelleft=True)
        ax.minorticks_on()
        #
        plt.text(21.6,12.1,cluster_names[cluster],fontsize=12)
        plt.text(21.6,11.5,'$M/L_{max}$: %s'%limiting_mass[cluster],fontsize=10)
        plt.text(21.6,10.5,members_string,fontsize=10)
        plt.text(18.1,6.1,'z = %s'%z_cluster[cluster],fontsize=10)
        plt.text(18.1,5.7,'$\sigma_{mass}$ = %s'%std_dev,fontsize=10)
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
        members_string = '# cluster members(F814W): %s'%total_mem+'(%s'%mem_shown+')'
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
    axs.flat[0].text(25.15,5.75,'$M/L_{max,F160W}$: %s'%limiting_mass[0],fontsize=12)
    axs.flat[1].text(25.15,6.05,cluster_names[1],fontsize=12)
    axs.flat[1].text(25.15,5.75,'$M/L_{max,F160W}$: %s'%limiting_mass[1],fontsize=12)
    axs.flat[2].text(25.15,6.05,cluster_names[2],fontsize=12)
    axs.flat[2].text(25.15,5.75,'$M/L_{max,F160W}$: %s'%limiting_mass[2],fontsize=12)
    axs.flat[3].text(25.15,6.05,cluster_names[3],fontsize=12)
    axs.flat[3].text(25.15,5.75,'$M/L_{max,F160W}$: %s'%limiting_mass[3],fontsize=12)
    axs.flat[4].text(25.15,6.05,cluster_names[4],fontsize=12)
    axs.flat[4].text(25.15,5.75,'$M/L_{max,F160W}$: %s'%limiting_mass[4],fontsize=12)
    axs.flat[5].text(25.15,6.05,cluster_names[5],fontsize=12)
    axs.flat[5].text(25.15,5.75,'$M/L_{max,F160W}$: %s'%limiting_mass[5],fontsize=12)
    #axs.clim(0,z_cutoff[1])
    cbar = fig.colorbar(sm, ax=axs[:, 2],location='right')#, vmin=0.0, vmax=z_cutoff[1])
    cbar.ax.set_title("|${\Delta}$z|$_{phot,F160W}$")
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
        print('Program "data_mass_completeness*.py" took: %s seconds to run.\n\n' % (time.time() - start_time))
#
#
#
#
print('\n\n"data_mass_completeness_F160W.py"  terminated successfully.\n')
#
#                        
###################     PROGRAM END     ###################