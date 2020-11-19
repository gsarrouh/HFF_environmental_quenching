# Created on Tue Aug 18 06:36:54 2020
#
#
### WHAT THIS PROGRAM DOES:
### This script reads in all data for the UltraVISTA/COSMOS catalogue and prepares data for the SMF reference field. Key information is summarized in the tables enabled by activating "adams_flag == 1"
#
### Data is organized in a single catalogue ("master_cat_uvc"). galaxies are selected using the "USE" column of the main photometric catalogue, and from these galaxies for the field sample are identified with a designation column "member == 0", selecting all galaxies in the redshift range "z_cutoff" above mass of "limiting_mass_uvc" in units of log_10(M/M_sol).
#
#
## FILTER 1 - 'member': :   identifies type of data each object has
#
## 0: secure field member above limiting_mass_uvc M_sol
## -99: other
#
#
#
### Section summary:
#
### PROGRAM START
#
### (0)    import modules, define functions & limiting mass, FLAGS;
### (1)    import DATA
### (1.1)   calculate rest-fram UVJ, & U-V, V-J colors
### (1.2)   SUMMARY TABLE 1
### (2)    separate by TYPE (add filter)
### (3)    separate by MEMBER (i.e. select field sample above limiting mass, add filter)
#
### PROGRAM END
#
#
### NOTE: To find a section, search "SECTION (*)" in all caps, where * is the section number you want. To find fields which require user input, search "MAY NEED TO EDIT" in all caps. Some of these fields may direct you to manually enter input into a sub-program.
#
#
#
###################     PROGRAM START     ###################
#
#
## SECTION (0): modules, definitions, FLAGS
#
## MODULES
import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import astropy
from astropy.table import Table
from astropy.table import Column
import time
import pandas as pd
#
#
## TIME_FLAG: START
## superior time_flag which supercedes all others and times the entire program
time_flag = 1     # 0= all timers off;   1=on, time entire program;    2=off, allow individual sections to be timed
#
#
if 'project_time_flag' in locals():
    pass
else:
    project_time_flag = 0
#
if time_flag == 1:
    if project_time_flag == 1:
        pass
    else:
        start_time = time.time()
#
## DEFINITIONS
#
## MAY NEED TO EDIT
limiting_mass_uvc = 9.0
#
## FLAGS !!!
## MAY NEED TO EDIT
#
## section summary flags   (all are turned on if "adams_flag"==1)
## these flags print section summary tables for bookkeeping purposes
summary_flag_1 = 1          # S1.2: display diagnostic summary table, describes SUB-type filter
summary_flag_2 = 1          # S2: display outlier fractions & UVJ table
summary_flag_3 = 1          # S3: display TYPE filter summary table (SF/Q)
## diagnostic flags:
diag_flag_1 = 0             # S2: histograms of del_z spec/phot
diag_flag_2 = 0
diag_flag_3 = 0
#
#
if 'project_diagnostic_flag' in locals():
    pass
else:
    project_diagnostic_flag = 2
    project_plot_flag = 2
    project_master_variational_flag = 0
    z_cutoff = [0.012,0.055]     # [spec,phot] cutoffs for cluster membership
    z_cutoff_field = [0.05,0.10]
    z_cluster = [0.396,0.543,0.545,0.375,0.348,0.308]
    num_bins_SF_pos_neg = [7.36,8.93,10.16,12.3]
    num_bins_Q_pos_neg = [7.36,9.43,10.11,10.47,12.3]
    limiting_mass_flag = 1
    bin_width = 0.4
    diagnostic_round_flag = 2
    mcmc_flag = 0
#
if 'adams_flag' in locals():
    pass
else:
    adams_flag = 1
#
# Read in ALL data from WORKING DIRECTORY: NSERC17/HFFtoAdam/working_data/HFF_environmental_quenching
#
#
#
##SECTION 1: import all data from HFF team, convert flux to luminosity & gather full
#
print('"UVC_master_data.py" Section 1: import data beginning...')
#
## import catalogues into single table "master_cat_uvc"; separate objects for which there is no redshift data (both photo & spec) as "nodata"
#
##create table from EAZY output redshift ".zout" file
z_uvc = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/UltraVISTA_catalgoues/UVISTA_final_v4.1.zout',format='ascii')
#
#create table from FAST ".fout" file (contains mass estimates)
f_uvc = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/UltraVISTA_catalgoues/UVISTA_final_BC03_v4.1.fout',format='ascii')
#
# rename columns of the .fout files because for some reason the column names didn't register
col_names_old = ['col1','col2','col3','col4','col5','col6','col7','col8','col9','col10','col11']
col_names_new = ['id','z','ltau','metal','lage','Av','lmass','lsfr','lssfr','la2t','chi2']
for ii in range(len(col_names_new)):
    f_uvc.rename_column(col_names_old[ii],col_names_new[ii])
#
##read in the whole bloody catalogue
cat_uvc = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/UltraVISTA_catalgoues/UVISTA_final_v4.1.cat',format='ascii')
#
##creat table for colours
UV_uvc = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/UltraVISTA_catalgoues/UVISTA_final_v4.1.153-155.rf',format='ascii')
VJ_uvc = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/UltraVISTA_catalgoues/UVISTA_final_v4.1.155-161.rf',format='ascii')
##aggregate into a single table
#
global master_cat_uvc
master_cat_uvc = Table([z_uvc['id'],z_uvc['z_peak'],z_uvc['z_spec'],UV_uvc['L153'],UV_uvc['L155'],VJ_uvc['L161'],f_uvc['lmass'],cat_uvc['USE']], names=('id','z_peak','z_spec','u','v','j','lmass','USE'))
#
## add "empty" columns (with value set = -99) for the remaining sieves: type, member for,[SF or Q],[field sample] respectively
E2 = Column([-99]*len(master_cat_uvc), name='type', dtype=np.int8)
E3 = Column([-99]*len(master_cat_uvc), name='member', dtype=np.int8)
master_cat_uvc.add_columns([E2,E3],[-1,-1])                   # add columns to the end of table
#
# summarize dataset
counting_array1 = np.array([0]*2)   # [USE==1, galaxies in redshift range of clusters]
#
for counter in range(len(master_cat_uvc)):
    if master_cat_uvc['USE'][counter] == 1:  #non-contaminated GALAXIES w/ good photometric measurements, Ks_tot < 24.4
        counting_array1[0]+=1
        if master_cat_uvc['z_peak'][counter] > min(z_cluster) and master_cat_uvc['z_peak'][counter] < max(z_cluster):     # identify galaxies broadly in the cluster redshift range
            counting_array1[1]+=1
#
#
#
## SECTION (1.1): convert FLUX TO MAGNITUDE; using well-known mag = -2.5*log_10(flux) + zero_pt. zero_pt = 25
#
# add columns for luminosity calculations
empty_u = Column([-99]*len(master_cat_uvc), name='L_u', dtype=np.float64)
empty_v = Column([-99]*len(master_cat_uvc), name='L_v', dtype=np.float64)
empty_j = Column([-99]*len(master_cat_uvc), name='L_j', dtype=np.float64)
empty_uv = Column([-99]*len(master_cat_uvc), name='uv', dtype=np.float64)
empty_vj = Column([-99]*len(master_cat_uvc), name='vj', dtype=np.float64)
master_cat_uvc.add_columns([empty_u,empty_v,empty_j,empty_uv,empty_vj])
#
## convert flux to magnitude (erroneously labelled as luminosity, e.g. L_u for magnitude in UV), get color indices U-V, V-J, add to table.
## the number of objects this calculation is performed on should be equal to the number of (phot_only + (phot+spec)) objects identified above
#
counting_array2 = np.array([0]*6) #[USE==1,]
#
for counter in range(len(master_cat_uvc)):
    if master_cat_uvc['USE'][counter] == 1:
        master_cat_uvc['L_u'][counter] = -2.5*np.log10(master_cat_uvc['u'][counter]) + 25
        master_cat_uvc['L_v'][counter] = -2.5*np.log10(master_cat_uvc['v'][counter]) + 25
        master_cat_uvc['L_j'][counter] = -2.5*np.log10(master_cat_uvc['j'][counter]) + 25
        master_cat_uvc['uv'][counter] = master_cat_uvc['L_u'][counter] - master_cat_uvc['L_v'][counter]
        master_cat_uvc['vj'][counter] = master_cat_uvc['L_v'][counter] - master_cat_uvc['L_j'][counter]
        counting_array2[0]+=1
#
#
## SECTION (1.2): SUMMARY table
#
if summary_flag_1 == 1 or adams_flag == 1:
        ## Summarize initial data stats in table
        uvc_names = Column(['Total','USE==1','~0.3<z<0.55'],name='Property')
        uvc0 = Column([len(master_cat_uvc),counting_array1[0],counting_array1[1]],name='Total')  # total column
        uvc_stats = Table([uvc_names,uvc0])
        #
        print('\nSummary Table 1 - Catalogue summary (UltraVISTA):')
        print(uvc_stats)
        #
#
#
print('"UVC_master_data.py" Section 1 complete.\n')
#
#
#
#
#
## SECTION (2): add FILTER TYPE: separate SF/Q for both subsamples;    filter name: 'type'  IDs below; SELECTION CRITERIA from Shipley et al. 2018, Section (5.3)
##   0 = stars;  1 = SF (star-forming);    2 = Q (quiscient);    3 = outliers
#
print('\n"UVC_master_data.py" Section 2: classifying galaxy TYPE as star-forming or quiescent...')
#
#
counting_array3 = np.array([0]*6)
SF_type_uvc = 0
Q_type_uvc = 0
#
for counter in range(len(master_cat_uvc)):
    if master_cat_uvc['USE'][counter] == 1:          # identifies "good" galaxies
        counting_array3[0]+=1
        if master_cat_uvc['vj'][counter] < 0.75:
            if master_cat_uvc['uv'][counter] < 1.3:
                master_cat_uvc['type'][counter] = 1             # identify STAR-FORMING galaxies, type=1
                SF_type_uvc+=1
            else:
                master_cat_uvc['type'][counter] = 2             # identify passive (QUIESCENT) galaxies, type=2
                Q_type_uvc+=1
                #
        elif master_cat_uvc['vj'][counter] >= 0.75:
            if master_cat_uvc['uv'][counter] < ( (0.8 * master_cat_uvc['vj'][counter]) + 0.7 ):
                master_cat_uvc['type'][counter] = 1             # identify STAR-FORMING galaxies, type=1
                SF_type_uvc+=1
            else:
                master_cat_uvc['type'][counter] = 2             # identify passive (QUIESCENT) galaxies, type=2
                Q_type_uvc+=1
            #
    #
#
#
## SECTION (2.1): SUMMARY table
##  summarize data TYPE population as segregated above, and display in a table
#
if summary_flag_2 == 1 or adams_flag == 1:
    ## Summarize initial data stats in table
    type_uvc_names = Column(['Total (USE==1)','SF','Q','SUM (SF+Q)'],name='Property')
    type_uvc0 = Column([counting_array3[0],SF_type_uvc,Q_type_uvc,np.sum([SF_type_uvc,Q_type_uvc])],name='Total')  # total column
    type_uvc_stats = Table([type_uvc_names,type_uvc0])
    #
    print('\nSummary Table 2 - Catalogue by TYPE (PAR): %s'%type_uvc_stats)
    #
#
print('\n"UVC_master_data.py" Section 2 complete.')
#
#
#
#
#
## SECTION (3) : apply MEMBER FILTER: 'member=1' identifies our field sample for the SMF
#
## compute redshift range of galaxies in cluster sample
if z_field_bounds_flag == 1:
    lower_bound = (min(z_cluster) - z_cutoff[1]) / (1 + z_cutoff[1])
    upper_bound = (max(z_cluster) + z_cutoff_field[1]) / (1 - z_cutoff_field[1])
    z_field_bounds = [lower_bound, upper_bound]
#
counting_array4 = np.array([0]*6)    # counting array to track total number of objects in each cluster
#
SF_field_uvc_list = []
Q_field_uvc_list = []
#
## isolate all galaxies (SF & Q) in the redshift range ~0.3 < z < ~0.55, for the field sample of the SMF
for counter in range(len(master_cat_uvc)):
    if master_cat_uvc['USE'][counter] == 1:
        counting_array4[0]+=1
        if master_cat_uvc['z_peak'][counter] < z_field_bounds[1] and master_cat_uvc['z_peak'][counter] > z_field_bounds[0]:
            master_cat_uvc['member'][counter] = 1                                   # member=1: preliminary member of FIELD SAMPLE
            counting_array4[1]+=1     # count field sample
            if master_cat_uvc['lmass'][counter] > limiting_mass_uvc:
                counting_array4[2]+=1      #count above limiting mass
                if master_cat_uvc['type'][counter] == 1:
                    counting_array4[3]+=1     # count SF
                    SF_field_uvc_list.append(master_cat_uvc['lmass'][counter])
                elif master_cat_uvc['type'][counter] == 2:
                    counting_array4[4]+=1     # count Q
                    Q_field_uvc_list.append(master_cat_uvc['lmass'][counter])
            #
            #######
        ####
    #####
#####
#
#
#
#
## SECTION (3.1): SUMMARY table
##  summarize data MEMBER population as segregated above, and display in a table
#
if summary_flag_3 == 1 or adams_flag == 1:
    ## Summarize initial data stats in table
    mem_uvc_names = Column(['Total (USE==1)','Total field sample (all masses)','Total w/ mass > 10^9','SF','Q'],name='Property')
    mem_uvc0 = Column([counting_array4[0],counting_array4[1],counting_array4[2],counting_array4[3],counting_array4[4]],name='Total')  # total column
    mem_uvc_stats = Table([mem_uvc_names,mem_uvc0])
    #
    print('\nSummary Table 3 - Catalogue by MEMBER (UVC):\n%s'%mem_uvc_stats,'\n\nNOTE: name of lists to be made into histograms for Field SMF are "SF/Q_field_uvc_list"')
    #
#
#
print('\n"UVC_master_data.py" Section 3 complete.')
#
#
#
#
#
### TEMPORARY
#
# num_points = int((round((range2[1]-range2[0])/bin_width))+1)       # compute # of data points;  bin_width set in "main_project_file.py"; range2 set in 'data_mass_completeness*.py'
# num_bins = np.linspace(range2[0],range2[1],num_points)#

#
## TIME_FLAG END
#
if time_flag == 1:
    if project_time_flag == 1:
        pass
    else:
        print('Program "UVC_master_data.py" took: %s seconds to run.\n\n' % (time.time() - start_time))
#
#
#
print('\n\n"UVC_master_data.py"  terminated successfully.\n')
#
#
#
###################     PROGRAM END     ###################
#
