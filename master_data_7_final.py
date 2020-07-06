# Created on Wed Jun 24 10:23:11 2020
#
#### UPDATE 06/24/20 ####
## The previous versions of this file all used HFF DR v3.5, a preliminary data release from the summer of 2017. The paper Shipley et al. 2018 (which I reference extensively in my work) corresponds to the finished data release, v3.9. That final, most up-tp-date data release is what is implemented in the following work. 
#
#
### WHAT THIS PROGRAM DOES:
### This script reads in all data for the Hubble Frontier Fields images and prepares data for plotting and analysis. Key information is summarized in the tables: 
###     **"sub_stats","type_stats","SF_spec_stats","Q_spec_stats","SF_phot_stats","Q_phot_stats"**
#
### Data is organized in a single catalogue ("master_cat"), and objects are identified through applying a series of "FILTERS" to designate key populations using a numerical designation for ease of writing code.
#
#
## FILTER 1 - 'cluster':  cluster catalogues are designated as follows: 
## 1: macs0416
## 2: macs1149
## 3: macs0717
## 4: abell370
## 5: abell1063
## 6: abell2744
#
## FILTER 2 - 'sub': :   identifies sub-type of data each object has
## 0: no data (no photometry or spectroscopy)
## 1: spectroscopy & photometry
## 2: photometry only 
## 3: spectroscopy only  (there's just 3 in macs0717) 
## 4: star
#
## FILTER 3 - 'type': :   identifies type of data each object has
## 0: star
## 1: star-forming (SF)
## 2: quiescent (Q)  
## 3: outliers (defined as |del_z\ > 0.15)
#
## FILTER 4 - 'member': :   identifies type of data each object has
## NOTE: this designation is for objects with phot only (sub =2) and spec&phot (sub=1) only, as membership determinaion requires a good 'z_phot' estimate. as such, member=2 & =3 are for spec&phot (sub=1) subsample only, as only they can be classified as false pos/neg
## 0: secure cluster member
## 1: secure field    <-- this comprises the sample of field galaxies to be 
##                        compared with cluster, at similar redshifts
## 2: false positive
## 3: false negative
## 4: field outlier   <-- objects well outside the redshift range of the clusters (e.g. z > 0.55)
## 5: BCGs            <-- identified in the last section of the program, over-writing MEMBER assignment from section 4
#
#
#
### Section summary:
#
### PROGRAM START
#
### (0)    import modules, define functions; in future aggregate FLAGS to Section (0) as well
### (1)    import data into single table, creating KEY TABLE: "master_cat" 
### (1.1)  add filter ("sieves") columns, apply SUB-TYPE FILTER in each 
###        cluster ["nodata", "phot_only","spec_only","both"], 
### (1.2)  add DIAG_FLAG_1: summarize in table "sub_stats"
### (1.3)  convert flux to mag.,
### (2)    calculate various del_z's, 
### (2.1)  identify outliers; add DIAG_FLAG_2: summarize outliers
### (2.2)  compute & DISPLAY OUTLIER FRACTION, SCATTER (i.e. std dev), 
###        and MEAN of |del_z|.
### (3)    distinguish b/w SF/Q: apply TYPE FILTER
### (3.1)  add DIAG_FLAG_3: summarize in table "type_stats"
### (4)    make membership cuts to spec samples (i.e. apply MEMBER FILTER),
###        add DIAG_FLAG_4: apply VARIATIONAL DIAGNOSTIC ANALYSIS to test 
###        different redshift cutoff (OUTPUT FILE: /Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/section_4_false_pos_neg_redshift_cuts_*.txt)
###        *** this section calls: "data_mass_completeness*.py" & 
###        "spec_completeness_binning*.py";
###        NOTE: DIAG_FLAG_4 must equal 0 to proceed beyond this point
### (4.1)  add DIAG_FLAG_5: summarize in table "SF_spec_stats" & "Q_spec_stats"
### (5)    make membership cuts to phot samples
### (5.1)  add DIAG_FLAG_6: summarize in table "SF_phot_stats" & "Q_phot_stats"
### (6)    check for brightest Cluster Galaxies (bCGs) in the catalogue; 
###        add DIAG_FLAG_7: report # of bCGs identified/removed/added
#
#########
#########  **** CONFIRM W/ AM WHETHER THESE NEED TO BE TAKEN OUT OF CATALOAGUE OR INCLUDED    ****
#########  **** Update: no bCGs in catalogue, but they NNED TO BE ADDED BACK. track them down ****
#########  ****                 DELETE THESE COMMENTS ONCE DONE!!!!!                          ****
#########
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
#from spec_membership_selection_file import spec_membership_selection
import pandas as pd
#
#
## TIME_FLAG: START
## superior time_flag which supercedes all others and times the entire program
time_flag = 1     # 0= all timers off;   1=on, time entire program;    2=off, allow individual sections to be timed
#
if time_flag == 1 or project_time_flag == 1:
    time_flag = 1
    start_time = time.time()
#
#
## FLAGS !!!
## MAY NEED TO EDIT
## MASTER diagnostic flag:   0=all flags off;   1=all flags on;    2=individual flags may be turned on
diag_flag_master = 1
#
variational_anaylsis_master_flag = 0         # this allows you to set the above master flag to 1, without entering the variational analysis
#
#
## section summary flags   (all are turned on if "adams_flag"==1)
## these flags print section summary tables for bookkeeping purposes
summary_flag_1 = 1          # S1.2: display diagnostic summary table, describes SUB-type filter
summary_flag_2 = 1          # S2: display outlier fractions & UVJ table
summary_flag_3 = 1          # S3: display TYPE filter summary table (SF/Q)
summary_flag_4 = 1          # S4: MEMBER-filter classification (SPEC subsample), assuming you don't run variational analysis
summary_flag_5 = 1          # S5: MEMBER-filter classification (PHOT subsample)
summary_flag_6 = 1          # S7: Full MEMBERSHIP CUT summary of full Parent sample (SPEC + PHOT subsamples)
summary_flag_7 = 1          # S6: bCGs
## diagnostic flags:
diag_flag_1 = 1             # S2: histograms of del_z spec/phot     
diag_flag_2 = 1 
diag_flag_3 = 1    
diag_flag_4 = 0             
minor_diag_flag_4 = 1
diag_flag_5 = 1             
diag_flag_6 = 1             
diag_flag_7 = 1             
#
if 'project_diagnostic_flag' in locals():
    pass
else:
    project_diagnostic_flag = 2
    project_master_variational_flag = 0
#
if 'adams_flag' in locals():
    pass
else:
    adams_flag = 1
#
# Read in ALL data from WORKING DIRECTORY: NSERC17/HFFtoAdam/working_data
#
##SECTION 1: import all data from HFF team, convert flux to luminosity & gather full 
#
print('"master_data*.py" Section 1: import data beginning...')
## TIME_FLAG_1 START
time_flag_1 = 1     # track & print time to execute current section
#
if time_flag_1 == 1 and time_flag == 2:
    start_time = time.time()
#  
## import catalogues into single table "master_cat"; separate objects for which there is no redshift data (both photo & spec) as "nodata"
#
##create table from EAZY output redshift ".zout" file
z_macs0416 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs0416clu_catalogs/macs0416clu_v3.9.zout',format='ascii')
z_macs1149 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs1149clu_catalogs/macs1149clu_v3.9.zout',format='ascii')
z_macs0717 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs0717clu_catalogs/macs0717clu_v3.9.zout',format='ascii')
z_abell370 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell370clu_catalogs/abell370clu_v3.9.zout',format='ascii')
z_abell1063 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell1063clu_catalogs/abell1063clu_v3.9.zout',format='ascii')   
z_abell2744 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell2744clu_catalogs/abell2744clu_v3.9.zout',format='ascii') 
#   
#create table from FAST ".fout" file (contains mass estimates)
f_macs0416 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs0416clu_catalogs/macs0416clu_v3.9.fout',format='ascii')
f_macs1149 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs1149clu_catalogs/macs1149clu_v3.9.fout',format='ascii')
f_macs0717 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs0717clu_catalogs/macs0717clu_v3.9.fout',format='ascii')
f_abell370 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell370clu_catalogs/abell370clu_v3.9.fout',format='ascii')
f_abell1063 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell1063clu_catalogs/abell1063clu_v3.9.fout',format='ascii')
f_abell2744 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell2744clu_catalogs/abell2744clu_v3.9.fout',format='ascii')
#
# rename columns of the .fout files because for some reason the column names didn't register
col_names_old = ['col1','col2','col3','col4','col5','col6','col7','col8','col9','col10','col11']
col_names_new = ['id','z','ltau','metal','lage','Av','lmass','lsfr','lssfr','la2t','chi2']
for ii in range(len(col_names_new)):
    f_macs0416.rename_column(col_names_old[ii],col_names_new[ii])
    f_macs1149.rename_column(col_names_old[ii],col_names_new[ii])
    f_macs0717.rename_column(col_names_old[ii],col_names_new[ii])
    f_abell370.rename_column(col_names_old[ii],col_names_new[ii])
    f_abell1063.rename_column(col_names_old[ii],col_names_new[ii])
    f_abell2744.rename_column(col_names_old[ii],col_names_new[ii])
#
##read in the whole bloody catalogue
cat_macs0416 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs0416clu_catalogs/hffds_macs0416clu_v3.9.cat',format='ascii')
cat_macs1149 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs1149clu_catalogs/hffds_macs1149clu_v3.9.cat',format='ascii')
cat_macs0717 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs0717clu_catalogs/hffds_macs0717clu_v3.9.cat',format='ascii')
cat_abell370 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell370clu_catalogs/hffds_abell370clu_v3.9.cat',format='ascii')
cat_abell1063 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell1063clu_catalogs/hffds_abell1063clu_v3.9.cat',format='ascii')
cat_abell2744 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell2744clu_catalogs/hffds_abell2744clu_v3.9.cat',format='ascii')
#
##creat table for colours
#macs0416
F0416_u = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs0416clu_catalogs/macs0416clu_v3.9.153.rf',format='ascii')
F0416_v = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs0416clu_catalogs/macs0416clu_v3.9.155.rf',format='ascii')
F0416_j = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs0416clu_catalogs/macs0416clu_v3.9.161.rf',format='ascii')
#macs1149
F1149_u = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs1149clu_catalogs/macs1149clu_v3.9.153.rf',format='ascii')
F1149_v = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs1149clu_catalogs/macs1149clu_v3.9.155.rf',format='ascii')
F1149_j = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs1149clu_catalogs/macs1149clu_v3.9.161.rf',format='ascii')
#macs0717
F0717_u = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs0717clu_catalogs/macs0717clu_v3.9.153.rf',format='ascii')
F0717_v = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs0717clu_catalogs/macs0717clu_v3.9.155.rf',format='ascii')   
F0717_j = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/macs0717clu_catalogs/macs0717clu_v3.9.161.rf',format='ascii')
#abell370
F370_u = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell370clu_catalogs/abell370clu_v3.9.153.rf',format='ascii')
F370_v = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell370clu_catalogs/abell370clu_v3.9.155.rf',format='ascii')
F370_j = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell370clu_catalogs/abell370clu_v3.9.161.rf',format='ascii')
#abell1063
F1063_u = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell1063clu_catalogs/abell1063clu_v3.9.153.rf',format='ascii')
F1063_v = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell1063clu_catalogs/abell1063clu_v3.9.155.rf',format='ascii')
F1063_j = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell1063clu_catalogs/abell1063clu_v3.9.161.rf',format='ascii')
#abell2744
F2744_u = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell2744clu_catalogs/abell2744clu_v3.9.153.rf',format='ascii')
F2744_v = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell2744clu_catalogs/abell2744clu_v3.9.155.rf',format='ascii')
F2744_j = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/abell2744clu_catalogs/abell2744clu_v3.9.161.rf',format='ascii')
##aggregate into a single table
macs0416 = Table([z_macs0416['id'],z_macs0416['z_peak'],z_macs0416['z_spec'],F0416_u['L153'],F0416_v['L155'],F0416_j['L161'],F0416_u['DM'],f_macs0416['lmass'],f_macs0416['lsfr'],f_macs0416['lssfr'],cat_macs0416['flux_radius'],cat_macs0416['star_flag'],cat_macs0416['use_phot'],cat_macs0416['f_F160W'],cat_macs0416['e_F160W'],cat_macs0416['flag_F160W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','use_phot','f_F160W','e_F160W','flag_F160W'))
macs1149 = Table([z_macs1149['id'],z_macs1149['z_peak'],z_macs1149['z_spec'],F1149_u['L153'],F1149_v['L155'],F1149_j['L161'],F1149_u['DM'],f_macs1149['lmass'],f_macs1149['lsfr'],f_macs1149['lssfr'],cat_macs1149['flux_radius'],cat_macs1149['star_flag'],cat_macs1149['use_phot'],cat_macs1149['f_F160W'],cat_macs1149['e_F160W'],cat_macs1149['flag_F160W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','use_phot','f_F160W','e_F160W','flag_F160W'))
macs0717 = Table([z_macs0717['id'],z_macs0717['z_peak'],z_macs0717['z_spec'],F0717_u['L153'],F0717_v['L155'],F0717_j['L161'],F0717_u['DM'],f_macs0717['lmass'],f_macs0717['lsfr'],f_macs0717['lssfr'],cat_macs0717['flux_radius'],cat_macs0717['star_flag'],cat_macs0717['use_phot'],cat_macs0717['f_F160W'],cat_macs0717['e_F160W'],cat_macs0717['flag_F160W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','use_phot','f_F160W','e_F160W','flag_F160W'))
abell370 = Table([z_abell370['id'],z_abell370['z_peak'],z_abell370['z_spec'],F370_u['L153'],F370_v['L155'],F370_j['L161'],F370_u['DM'],f_abell370['lmass'],f_abell370['lsfr'],f_abell370['lssfr'],cat_abell370['flux_radius'],cat_abell370['star_flag'],cat_abell370['use_phot'],cat_abell370['f_F160W'],cat_abell370['e_F160W'],cat_abell370['flag_F160W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','use_phot','f_F160W','e_F160W','flag_F160W'))
abell1063 = Table([z_abell1063['id'],z_abell1063['z_peak'],z_abell1063['z_spec'],F1063_u['L153'],F1063_v['L155'],F1063_j['L161'],F1063_u['DM'],f_abell1063['lmass'],f_abell1063['lsfr'],f_abell1063['lssfr'],cat_abell1063['flux_radius'],cat_abell1063['star_flag'],cat_abell1063['use_phot'],cat_abell1063['f_F160W'],cat_abell1063['e_F160W'],cat_abell1063['flag_F160W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','use_phot','f_F160W','e_F160W','flag_F160W'))
abell2744 = Table([z_abell2744['id'],z_abell2744['z_peak'],z_abell2744['z_spec'],F2744_u['L153'],F2744_v['L155'],F2744_j['L161'],F2744_u['DM'],f_abell2744['lmass'],f_abell2744['lsfr'],f_abell2744['lssfr'],cat_abell2744['flux_radius'],cat_abell2744['star_flag'],cat_abell2744['use_phot'],cat_abell2744['f_F160W'],cat_abell2744['e_F160W'],cat_abell2744['flag_F160W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','use_phot','f_F160W','e_F160W','flag_F160W'))
#
#
## NOTE: cat_*['flag_F160W'] identifies BCGs IN EACH FILTER. ***** CONFIRM w/ AM ***** that I'm using the right filter; BCGs are identified as ['flag_F160W']==4, see Shipley et al 2018. Section 3.9
#
## create columns and append to master_cat to identify sub-sample, type, 
## and fraction for catalogue: 
##   type: 0 = stars,  1 = SF,  2 = Q
##   member:  0 = not in spec sub-sample,  1 = secure in cluster,  2 = secure in field
##            3 = false positive,  4 = false negative
## Note: "member" only applies to the spec sub-sample
#
D1 = Column([1]*len(macs0416),name='cluster')
D2 = Column([2]*len(macs1149),name='cluster')        #cluster designation columns
D3 = Column([3]*len(macs0717),name='cluster')
D4 = Column([4]*len(abell370),name='cluster')
D5 = Column([5]*len(abell1063),name='cluster')
D6 = Column([6]*len(abell2744),name='cluster')
#
macs0416.add_column(D1)       
macs1149.add_column(D2)
macs0717.add_column(D3)
abell370.add_column(D4)
abell1063.add_column(D5)
abell2744.add_column(D6)
#
global master_cat
master_cat = Table(np.concatenate((macs0416,macs1149,macs0717,abell370,abell1063,abell2744), axis=0))  #create a master catalogue of all clusters
cluster_names = ['M0416','M1149','M0717','A370','A1063','A2744']
#
## add "empty" columns (with value set = 99) for the remaining sieves: sub, type, member for [phot or spec],[SF or Q],[cluster member, field, false pos/neg] respectively
E1 = Column([-99]*len(master_cat), name='sub', dtype=np.int8)    # create columns
E2 = Column([-99]*len(master_cat), name='type', dtype=np.int8)
E3 = Column([-99]*len(master_cat), name='member', dtype=np.int8)
master_cat.add_columns([E1,E2,E3],[-1,-1,-1])                   # add columns to the end of table
#
#
## SECTION (1.1) - SUB-TYPE FILTER: this section classifies all objects as either:
##  sub = 0: no data (neither spec. nor phot.)
##  sub = 1: both
##  sub = 2: phot only
##  sub = 3: spec only
##  sub = 4: stars
#
# sift for targets without no data, spec only, phot only, and both
spec_only = np.array([0]*6)    # to keep track by cluster
phot_only = np.array([0]*6)
use_phot = np.array([[0]*6]*2)      # row1=good (use_phot = 1);  row2=bad (use_phot !=1)
error = 0
other = 0
both = np.array([0]*6)
no_data = np.array([0]*6)
stars_sub = np.array([0]*6)
#
for counter in range(len(master_cat)):
    if master_cat['z_spec'][counter] > 0 and master_cat['z_peak'][counter] > 0:  #entries w/ both spectroscopy and photometry  (SPECTROSCOPIC sub-sample)
        if master_cat['use_phot'][counter] ==0:     # use_phot = 0 means bad photometry
            for jj in range(len(spec_only)):
                if master_cat['cluster'][counter] == (jj+1):      # identify # of spec only objects by cluster
                    #use_phot[1][jj]+=1
                    if master_cat['star_flag'][counter] == 1:
                        master_cat['sub'][counter] = 4          # APPLY FILTER: sub=4 for STARS
                        stars_sub[jj]+=1
                    else:
                        spec_only[(jj)]+=1
                        master_cat['sub'][counter] = 3          # APPLY FILTER: sub=3 for objects w/ SPEC ONLY
        elif master_cat['use_phot'][counter] ==1:     # use_phot = 1 means good photometry
            for jj in range(len(both)):
                if master_cat['cluster'][counter] == (jj+1):      # identify # of (spec & phot) objects by cluster
                    #use_phot[0][jj]+=1
                    if master_cat['star_flag'][counter] == 1:
                        master_cat['sub'][counter] = 4          # APPLY FILTER: sub=4 for STARS; THIS SHOULDNT DO ANYTHING GIVEN THE "use_phot==1" DESIGNATION
                        stars_sub[jj]+=1
                    else:
                        both[(jj)]+=1
                        master_cat['sub'][counter] = 1          # APPLY FILTER: sub=1 for objects w/ BOTH SPEC & PHOT
        else: error+=1                              # just to keep track of erroneous objects
    elif master_cat['z_spec'][counter] > 0 and master_cat['z_peak'][counter] < 0:  #entries w/ spectroscopy alone
        for jj in range(len(spec_only)):
            if master_cat['cluster'][counter] == (jj+1):      # identify # of spec only objects by cluster
                if master_cat['star_flag'][counter] == 1:
                    master_cat['sub'][counter] = 4          # APPLY FILTER: sub=4 for STARS
                    stars_sub[jj]+=1
                else:
                    spec_only[(jj)]+=1
                    master_cat['sub'][counter] = 3          # APPLY FILTER: sub=3 for objects w/ SPEC ONLY
    elif master_cat['z_spec'][counter] < 0 and master_cat['z_peak'][counter] > 0:  #entries w/ photometry alone (PHOTOMETRIC sub-sample)
        if master_cat['use_phot'][counter] ==0:
            for jj in range(len(both)):
                if master_cat['cluster'][counter] == (jj+1):      # identify # of NO DATA objects by cluster
                    #use_phot[1][jj]+=1
                    if master_cat['star_flag'][counter] == 1:
                        master_cat['sub'][counter] = 4          # APPLY FILTER: sub=4 for STARS
                        stars_sub[jj]+=1
                    else:
                        no_data[(jj)]+=1
                        master_cat['sub'][counter] = 0          # APPLY FILTER: sub=0 for objects w/ NEITHER SPEC NOR PHOT
        elif master_cat['use_phot'][counter] ==1:
            for jj in range(len(both)):
                if master_cat['cluster'][counter] == (jj+1):      # identify # of phot only objects by cluster
                    #use_phot[0][jj]+=1
                    if master_cat['star_flag'][counter] == 1:
                        master_cat['sub'][counter] = 4          # APPLY FILTER: sub=4 for STARS
                        stars_sub[jj]+=1
                    else:
                        phot_only[(jj)]+=1
                        master_cat['sub'][counter] = 2          # APPLY FILTER: sub=2 for objects w/ PHOT ONLY
    elif master_cat['z_spec'][counter] < 0 and master_cat['z_peak'][counter] < 0:  #entries w/ no z estimates at all
        for jj in range(len(both)):
            if master_cat['cluster'][counter] == (jj+1):      # identify # of NO DATA objects by cluster
                if master_cat['star_flag'][counter] == 1:
                    master_cat['sub'][counter] = 4          # APPLY FILTER: sub=4 for STARS
                    stars_sub[jj]+=1
                else:
                    no_data[(jj)]+=1
                    master_cat['sub'][counter] = 0          # APPLY FILTER: sub=0 for objects w/ NEITHER SPEC NOR PHOT
    if master_cat['use_phot'][counter] == 1:
        for jj in range(len(both)):
            if master_cat['cluster'][counter] == (jj+1):
                use_phot[0][jj]+=1
    else: 
        for jj in range(len(both)):
            if master_cat['cluster'][counter] == (jj+1):
                use_phot[1][jj]+=1
#
#
## SECTION (1.2): SUMMARY table
#
if summary_flag_1 == 1 or adams_flag == 1:
        ## Summarize initial data stats in table
        sub_names = Column(['FULL PHOT (Parent)','Total','spec & phot','only phot','spec only','no data','stars','SUM'],name='Property')
        col_names = cluster_names
        sub0 = Column([np.sum([phot_only,both]),np.sum([np.sum(both),np.sum(phot_only),np.sum(spec_only),np.sum(no_data),np.sum(stars_sub)]),np.sum(both),np.sum(phot_only),np.sum(spec_only),np.sum(no_data),np.sum(stars_sub),np.sum([both,phot_only,spec_only,no_data,stars_sub])],name='Total')  # total column
        sub_stats = Table([sub_names,sub0])
        for ii in range(len(spec_only)):
            sub_col = Column([np.sum([phot_only[ii],both[ii]]),np.sum([both[ii],phot_only[ii],spec_only[ii],no_data[ii],stars_sub[ii]]),both[ii],phot_only[ii],spec_only[ii],no_data[ii],stars_sub[ii],np.sum([both[ii],phot_only[ii],spec_only[ii],no_data[ii],stars_sub[ii]])],name=col_names[ii])  # add columns to table one cluster at a time
            sub_stats.add_column(sub_col)
        #
        print('\nSummary Table 1 - Catalogue by SUB-type:')
        print(sub_stats)
        print('\nOther skipped objects: %s'%other,'\nNOTE: "use_phot==1": %s'%np.sum(use_phot[0]),';  "use_phot==0": %s'%np.sum(use_phot[1]))
        print('NOTE: phot only + (spec+phot) samples are the FULL PHOT (Parent) sample w/: %s'%np.sum(phot_only+both),'\nNOTE: Difference b/w Parent sample & "use_phot==1": %s'%np.abs((np.sum(phot_only+both)-np.sum(use_phot[0]))))
#
#
## SECTION (1.3): convert FLUX TO MAGNITUDE; using well-known mag = -2.5*log_10(flux) + zero_pt. zero_pt = 25
#
# add columns for luminosity calculations
empty_u = Column([99]*len(master_cat), name='L_u', dtype=np.float64)
empty_v = Column([99]*len(master_cat), name='L_v', dtype=np.float64)
empty_j = Column([99]*len(master_cat), name='L_j', dtype=np.float64)
empty_uv = Column([99]*len(master_cat), name='uv', dtype=np.float64)
empty_vj = Column([99]*len(master_cat), name='vj', dtype=np.float64)
master_cat.add_columns([empty_u,empty_v,empty_j,empty_uv,empty_vj])
#
## convert flux to magnitude (erroneously labelled as luminosity, e.g. L_u for magnitude in UV), get color indices U-V, V-J, add to table.
## the number of objects this calculation is performed on should be equal to the number of (phot_only + (phot+spec)) objects identified above
#
objects_99 = np.array([0]*6)
use_phot = np.array([[0]*6]*2)      # row1=good (use_phot = 1);  row2=bad (use_phot !=1)
skipped_UVJ = np.array([0]*6)    # to track stars & objects w/o photometry
good_objects = np.array([0]*6)
#
for counter in range(len(master_cat)):
    if master_cat['use_phot'][counter] == 1:
        for jj in range(len(skipped_UVJ)):
            if master_cat['cluster'][counter] == (jj+1):      # identify # of NO DATA objects by cluster
                use_phot[0][jj]+=1
                if master_cat['sub'][counter] ==0 or master_cat['sub'][counter] ==3 or master_cat['sub'][counter] ==4:      # skip objects w/ "no data" (sub=0); "spec only" (sub=3); "stars" (sub=4)
                    skipped_UVJ[jj]+=1     # objects not in our FULL PHOT sample
                elif master_cat['sub'][counter] == -99:
                    objects_99[jj]+=1      # error objects that were not classified in section 1, if they even exist?
                elif master_cat['sub'][counter] == 1 or master_cat['sub'][counter] == 2:   # this should select all "use_phot=1" objects
                    master_cat['L_u'][counter] = -2.5*np.log10(master_cat['u'][counter]) + 25
                    master_cat['L_v'][counter] = -2.5*np.log10(master_cat['v'][counter]) + 25
                    master_cat['L_j'][counter] = -2.5*np.log10(master_cat['j'][counter]) + 25
                    master_cat['uv'][counter] = master_cat['L_u'][counter] - master_cat['L_v'][counter]
                    master_cat['vj'][counter] = master_cat['L_v'][counter] - master_cat['L_j'][counter]
                    good_objects[jj]+=1

    else:
        for jj in range(len(use_phot[1])):
            if master_cat['cluster'][counter] == (jj+1):     
                use_phot[1][jj]+=1
#################
#
#
if summary_flag_2 == 1 or adams_flag == 1:
        ## Summary table of U V J and excesses calc., total numbers should agree w/ those reported in "sub_stats"
        UVJ_names = Column(['Parent sample','Objects w/ phot','Objects w/o phot','Objects not classified','SUM'],name='Property')
        col_names = cluster_names
        UVJ0 = Column([np.sum([phot_only,both]), np.sum(good_objects),np.sum(skipped_UVJ),np.sum(objects_99),np.sum([good_objects,skipped_UVJ,objects_99])],name='Total')  # total column
        UVJ_stats = Table([UVJ_names,UVJ0])
        for ii in range(len(spec_only)):
            UVJ_col = Column([np.sum([phot_only[ii],both[ii]]),good_objects[ii],skipped_UVJ[ii],objects_99[ii],np.sum([good_objects[ii],skipped_UVJ[ii],objects_99[ii]])],name=col_names[ii])  # add columns to table one cluster at a time
            UVJ_stats.add_column(UVJ_col)
        #####
        print('\nUVJ and excesses calculation (identified w/ "use_phot==1": \n%s'%UVJ_stats,'\n\nNOTE: "objects w/o phot" are those w/ "use_flag" == 1 that have been classified as no good data/bad phot.\nNote: "use_phot" = 1: %s'%np.sum(use_phot[0]),'; "use_phot" = 0: %s'%np.sum(use_phot[1]))
#########
#
print('"master_data*.py" Section 1 complete.\n')
#
## TIME_FLAG_1 END
#
if time_flag_1 == 1 and time_flag == 2:
    print("'master_data*.py' Section 1 took: %s seconds.\n\n" % (time.time() - start_time))
#
#
#
#
#
## SECTION (2): calulate DEL_Z's & separate OUTLIERS
#
print("\n'master_data*.py' Section 2: calculating del_z's and identifying outliers...\n")
#
## TIME_FLAG_2 START
time_flag_2 = 1     # track & print time to execute current section
#
if time_flag_2 == 1 and time_flag == 2:
    start_time = time.time()
#  
#####Note: the photometric redshift used is 'z_peak' column from data
#
#(i) calculate delta_z for targets w/ both z_spec & z_phot
master_cat.add_column(Column([-99]*len(master_cat),name='del_z', dtype=np.float64))           # del_z = (z_phot - z_spec) / (1 + z_spec)
master_cat.add_column(Column([-99]*len(master_cat),name='z_clusterspec', dtype=np.float64))   # del_z = (z_spec - z_cl) / (1 + z_spec)
master_cat.add_column(Column([-99]*len(master_cat),name='z_clusterphot', dtype=np.float64))   # del_z = (z_phot - z_cl) / (1 + z_phot)
#
# store cluster redshifts; obtained from https://archive.stsci.edu/prepds/frontier/
z_cluster = [0.396,0.543,0.545,0.375,0.348,0.308]
# cluster_names = ['M0416','M1149','M0717','A370','A1063','A2744']
#
## calucalte del_z, z_clusterspec, z_clusterphot for outlier cut (defined above); these will be used to make cuts (member, field, false pos/neg) to spec sample, from which we will use relative fractions by mass to correct the photometric sample for completeness.
#
#
spec_subsample = np.array([0]*6)   # to track objects w/ BOTH (spec + phot), by cluster
phot_subsample = np.array([0]*6)   # to track objects w/ ONLY PHOT
skipped_delz = np.array([0]*6)     # to track all other objects
#
for counter in range(len(master_cat)):
    if master_cat['sub'][counter] == 1:   # sub=1 identifies spec&phot subsample
        master_cat['del_z'][counter] = ((master_cat['z_peak'][counter] - master_cat['z_spec'][counter]) / (1 + master_cat['z_spec'][counter]))
        for ii in range(len(z_cluster)):
            if master_cat['cluster'][counter] == (ii+1):
                spec_subsample[ii]+=1
                master_cat['z_clusterspec'][counter] = ((master_cat['z_spec'][counter] - z_cluster[ii]) / (1 + master_cat['z_spec'][counter]))
                master_cat['z_clusterphot'][counter] = ((master_cat['z_peak'][counter] - z_cluster[ii]) / (1 + master_cat['z_peak'][counter]))
    elif master_cat['sub'][counter] == 2:   # sub=2 identifies phot-only subsample
        for ii in range(len(z_cluster)):
            if master_cat['cluster'][counter] == (ii+1):
                phot_subsample[ii]+=1
                master_cat['z_clusterphot'][counter] = ((master_cat['z_peak'][counter] - z_cluster[ii]) / (1 + master_cat['z_peak'][counter]))
    else:
        for ii in range(len(z_cluster)):
            if master_cat['cluster'][counter] == (ii+1):
                skipped_delz[ii]+=1
                #
            #####
        #####
#########
#
## OPTIONAL DIAGNOSTIC: make a histogram of what was just calculated
if (diag_flag_1 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    exec(open('delz_hist_diagnostic.py').read())
#
#
## SECTION (2.1): separate OUTLIERS from both phot & spec sub-sample, defined as |del_z| < 0.15. apply FILTER TYPE = 3 for outliers;   
#
outliers = np.array([0]*6)      # initialize array to track outliers by cluster, for computing outlier fraction later
sum_delz = []                   # for computing mean |del_z|
count_stars = 0
#
for counter in range(len(master_cat)):
    if master_cat['sub'][counter] == 1:                        # sub=1 for objects with both spec & phot; total # of such objects tracked by cluster above in the array called "both"
        if np.abs(master_cat['del_z'][counter]) > 0.15:        # |del_z| > 0.15 for outliers, threshold chosen for historical reasons to facilitate comparison with other studies
            master_cat['type'][counter] = 3                  # type=3 identifies outliers
            for ii in range(len(outliers)):
                if master_cat['cluster'][counter] == (ii+1):   # keep track of outliers by cluster
                    outliers[ii]+=1
            sum_delz.append(np.abs(master_cat['del_z'][counter]))    # keep track of all |del_z| measurements for stats computation
        else:
            sum_delz.append(np.abs(master_cat['del_z'][counter]))
        if master_cat['type'][counter] == 3 and master_cat['star_flag'][counter] == 1:                  # overwrite designation for stars
            master_cat['type'][counter] = 0                        # type=0 for stars
            outliers-=1
            count_stars+=1
#
#
## SECTION (2.2): compute & DISPLAY OUTLIER FRACTION, SCATTER (i.e. std dev), and MEAN of |del_z|.
#
if summary_flag_2 == 1 or adams_flag == 1:
    ## Summary table of U V J and excesses calc., total numbers should agree w/ those reported in "sub_stats"
    delz_names = Column(['Parent sample','Spec subsample','Phot subsample','SUM'],name='Property')
    col_names = cluster_names
    delz0 = Column([np.sum([phot_only,both]), np.sum(spec_subsample),np.sum(phot_subsample),np.sum([spec_subsample,phot_subsample])],name='Total')  # total column
    delz_stats = Table([delz_names,delz0])
    for ii in range(len(spec_only)):
        delz_col = Column([np.sum([phot_only[ii],both[ii]]),spec_subsample[ii],phot_subsample[ii],np.sum([spec_subsample[ii],phot_subsample[ii]])],name=col_names[ii])  # add columns to table one cluster at a time
        delz_stats.add_column(delz_col)
    #####
    print('\nSummary Table 2: delz calculation: \n%s'%delz_stats)
    print('\nNOTE: SUM + skipped objects: %s'%np.sum([spec_subsample,phot_subsample]),' + %s'%np.sum(skipped_delz),' = %s'%np.sum([spec_subsample,phot_subsample,skipped_delz]))
    ###
    delz_mean = np.mean(sum_delz)
    delz_scatter = np.std(sum_delz)
    print('\nOUTLIERS total (in spec subsample): %s' % np.sum(outliers),'\nOutlier fraction: %s' % (np.sum(outliers)/np.sum(both)))
    print('|del_z| mean: %s'%delz_mean)
    print('|del_z| scatter: %s\n'%delz_scatter)
    print('NOTE: the new "use-able" number for the spec subsample is: %s'%(np.sum(both)-np.sum(outliers)))
#
#
## TIME_FLAG_2 END
#
if time_flag_2 == 1 and time_flag == 2:
    print('"master_data*.py"Section 2 took: %s seconds.\n\n' % (time.time() - start_time))
else:
    print('"master_data*.py" Section 2 complete.')
#  
#
#
#
#
## SECTION (3): add FILTER TYPE: separate SF/Q for both subsamples based on van der Burg (2013) 
## colour criteria;    filter name: 'type'  IDs below
##   0 = stars;  1 = SF (star-forming);    2 = Q (quiscient);    3 = outliers
#
print('\n"master_data*.py" Section 3: classifying galaxy TYPE as star-forming or quiescent...')
#
## TIME_FLAG_3 START
time_flag_3 = 1     # track & print time to execute current section
#
if time_flag_3 == 1 and time_flag == 2:
    start_time = time.time()
#  
SF_type = np.array([[0]*6]*2)                         # initialize arrays;  row1=spec;  row2=phot
Q_type = np.array([[0]*6]*2)
stars_type = np.array([0]*6)                          # count # of stars by cluster
lost_type = np.array([0]*6)                           # objects lost due no data (sub=0), spec only (sub=3), 
#
for counter in range(len(master_cat)):
    if master_cat['star_flag'][counter] == 1:          # identify STARS, filter type=0
        master_cat['type'][counter] = 0              
        for ii in range(len(stars_type)):
            if master_cat['cluster'][counter] == (ii+1):   # keep stars of outliers by cluster
                stars_type[ii]+=1
    elif master_cat['sub'][counter]==1 or master_cat['sub'][counter]==2:    # sub=1 for (spec & phot) subsample, sub=2 for phot subsample; i.e. look at all objects with photometry   <-- THIS IS THE LINE THAT MATTERS!
        if master_cat['uv'][counter] > 1.3 and master_cat['vj'][counter] < 1.6 and master_cat['uv'][counter] > ((0.88*master_cat[counter]['vj']) + 0.6): 
            if master_cat['type'][counter] !=3:             # skip outliers
                master_cat['type'][counter] = 2             # identify passive (QUIESCENT) galaxies, type=2
                for ii in range(len(Q_type[0])):
                    if master_cat['cluster'][counter] == (ii+1):   # keep stars of outliers by cluster
                        if master_cat['sub'][counter] == 1:
                            Q_type[0][ii]+=1
                        else: Q_type[1][ii]+=1
        else:
            if master_cat['type'][counter] !=3:             # skip outliers    
                master_cat['type'][counter] = 1             # identify STAR-FORMING galaxies, type=1
                for ii in range(len(SF_type[0])):
                    if master_cat['cluster'][counter] == (ii+1):   # keep 
                        if master_cat['sub'][counter] == 1:
                            SF_type[0][ii]+=1
                        else:
                            SF_type[1][ii]+=1    
    else:
        for ii in range(len(lost_type)):
            if master_cat['cluster'][counter] == (ii+1):   # keep stars of outliers by cluster
                lost_type[ii]+=1      # objects lost due to spec only (sub=3),no data (sub=0), and possibly outliers (type=3) 
#        
#
## SECTION (3.1): SUMMARY table
##  summarize data TYPE population as segregated above, and display in a table
#
if summary_flag_3 == 1 or adams_flag == 1:
    ## Summarize initial data stats in table
    type_names = Column(['Parent sample','SF - total','SF - phot','SF - spec','Q - total','Q - phot','Q - spec','Outliers','SUM (less totals)'],name='Property')
    col_names = cluster_names
    type0 = Column([np.sum([phot_only,both]),np.sum(SF_type),np.sum(SF_type[1]),np.sum(SF_type[0]),np.sum(Q_type),np.sum(Q_type[1]),np.sum(Q_type[0]),np.sum(outliers),np.sum([np.sum(SF_type),np.sum(Q_type),np.sum(outliers)])],name='Total')  # total column
    type_stats = Table([type_names,type0])
    for ii in range(len(spec_only)):
        type_col = Column([np.sum([phot_only[ii],both[ii]]),(SF_type[1][ii]+SF_type[0][ii]),SF_type[1][ii],SF_type[0][ii],(Q_type[1][ii]+Q_type[0][ii]),Q_type[1][ii],Q_type[0][ii],outliers[ii],np.sum([SF_type[1][ii],SF_type[0][ii],Q_type[1][ii],Q_type[0][ii],outliers[ii]])],name=col_names[ii])
        type_stats.add_column(type_col)  # add columns to table one cluster at a time
    #
    print('Summary Table 3 - Catalogue by TYPE:')
    print(type_stats)
    print('\nNOTE: "Other" is comprised of objects without photometry (i.e. stars, and objects with either bad photometry or both bad photometry and bad spectroscopy).\n')
#
## TIME_FLAG_3 END
#
if time_flag_3 == 1 and time_flag == 2:
    print('"master_data*.py" Section 3 took: %s seconds.\n\n' % (time.time() - start_time))
else:
    print('"master_data*.py" Section 3 complete.')
#
#
#
#
#
## SECTION (4) : apply MEMBER FILTER based on cuts discussed with AM in 2017 (see comments below), based on cuts made in VDB2013. MEMBER isolates: 0=cluster member (secure); 1=field (secure); 2=false pos; 3=false neg; UPDATE: write an iterative loop which tests varying values for redshift membership cutoffs, bins false pos/neg, determines spec. completeness correction factors, measure them against each other, rank them from "best to worst", and print result to an output document. Then INSPECT BY EYE and pick the one you want, and hard code the cuts after the diag_flag_* loop
#
print('\n"master_data*.py" Section 4: spectroscopic membership cuts and classifying cluster members, field, false pos/neg...')
## CRITERIA:
#
## SF: cluster = 0: abs(z_clusterspec) < 0.01 & abs(z_cluster phot) < 0.03; 
##      field = 1: abs(z_clusterspec) > 0.02 & abs(z_cluster phot) > 0.1; 
## false pos = 2: abs(z_clusterspec) > 0.01 & abs(z_cluster phot) < 0.03;
## false neg = 3: abs(z_clusterspec) < 0.01 & abs(z_cluster phot) > 0.03;
## Q: same cutoff for z_clusterspec, cutoff for z_clusterphot > 0.07
#
## Note: this is only done for spec sample; results used to correct photo sample for completeness at end of analysis in file "master_smf*.py"
#
### SF sub-sample: SF_spec into secure member, field, false pos/neg
#
#
## TIME_FLAG_4 START
time_flag_4 = 1     # 0=off;   1=on;   track & print time to execute current section
#
if time_flag_4 == 1 and time_flag == 2:
    start_time = time.time()
#
#
##  test different redshift cutoffs and print result to output document
#
#
space = ' '
if variational_anaylsis_master_flag == 1 or project_master_variational_flag == 1:
    #
    print('\n\nDIAG_4: ENTERING redshift cutoff VARIATIONAL DIAGNOSTIC...\n\n')
    ## the following code is part of a diagnostic to test the # of false pos/neg produced by altering the photometric/spectroscopic redshift cutoff. upon completion it will be commented out permenantly.
    #
    ## MAY NEED TO EDIT: "increment" & "*cutoff_range", for both SPEC (immediately below) and PHOT (after the first 'for' loop)
    ## SPEC cuts
    increment = np.array([0.01,0.01])   # [spec,phot]
    z_cutoff_spec_range = np.array([0.010,(0.050+increment[0])])
    # define cut-offs for SF & Q
    z_cutoff_spec = np.round(np.arange(z_cutoff_spec_range[0],z_cutoff_spec_range[1],increment[0]),decimals=3)       # create array from [0.01,0.05] in steps of "increment"
    #
    #
    ## these lists don't have a purpose but might come in handy later
    ## cluster mass lists to see false pos/neg by mass bin by cluster
    #SF_pos = [[],[],[],[],[],[]]
    #SF_neg = [[],[],[],[],[],[]]
    #Q_pos = [[],[],[],[],[],[]]
    #Q_neg = [[],[],[],[],[],[]]    
    #
    #
    ## DIAGNOSTIC loop to test different values of z_cutoff STARTS here
        # open a file to print to
    #f = open('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/diagnostic_outputs/master_data_section_4_false_pos_neg_redshift_cuts_VARIATIONAL.txt','w+')
    #
    ## explain what the file does
    #header1 = '\n\n### The following output lists cluster membership, false\npos/neg, & overall false pos/neg ratio (see column\nnames below). The program varies the spectroscopic\nmembership cutoff in the range [0.01,0.05] in increments\nof 0.005, testing all photometric cutoffs in the\nrange [0.01,0.09] in increments of 0.01.'
    #header2 = '\n#\n#\n### Columns:  del_z cutoff  Type   Member   False pos.   False neg.   % acceptance   #   False pos/neg ratio\n#\n#\n'
    #writer = asterisks+'%s\n'%header1+asterisks+header2+'\n'
    #f.write(writer)
    #
    ## create a list to aggregate all the different dataframes (created from importing results of files created below)
    df_list = []
    ## create file looping through different values for spec/phot membeship cutoff
    for cutoff_spec in range(len(z_cutoff_spec)):
        #
        if cutoff_spec == 10:
            break
        ## print diagnostic progress
        if ii%5 == 0:
            progress = (cutoff_spec / len(z_cutoff_spec)*100)    # calculate progress through z_cutoff_spec_range
            print('%s'%ii,'th cut of %s'%len(z_cutoff_spec),' being investigated.\nz_spec cut: %.3f'%z_cutoff_spec[cutoff_spec],'; in range %.3f'%z_cutoff_spec[0],' to %.3f'%z_cutoff_spec[-1],'\nCurrent progress: %.3f'%progress,'%')
        ## PHOT cuts
        z_cutoff_phot_range = [max(z_cutoff_spec[cutoff_spec],0.03),(0.09+(increment[1]))]  ## MAY NEED TO EDIT phot range
        z_cutoff_phot = np.round(np.arange(z_cutoff_phot_range[0],z_cutoff_phot_range[1],increment[1]),decimals=3) # create array from [0.01,0.05] in steps of increment_phot; replace in loop below with z_cutoff once cutoffs are determined
        for cutoff_phot in range(len(z_cutoff_phot)):
            #
            ## define spec & phot cutoffs
            z_cutoff = np.array([z_cutoff_spec[cutoff_spec],z_cutoff_phot[cutoff_phot]])
            #
            ## check it directories to store outputs exist. if not, create them
            output_dir = '/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/diagnostic_outputs/spec_binning/z_spec_%.3f'%z_cutoff[0]
            check_folder = os.path.isdir(output_dir)
            #
            ## If folder doesn't exist, then create it.
            if not check_folder:
                os.makedirs(output_dir)
                if project_diagnostic_flag == 1 or diag_flag_master == 2:
                    print("\nCreated folder : "+output_dir)
            else:
                pass#print(output_dir, "\nfolder already exists.")
            #
            #
            #
            ########################
            ########################
            # WE NOW CALL SEPARATE PROGRAMS TO ACHIEVE SEPARATE TASKS
            ########################
            ########################
            #
            #
            #
            ## Call "spec_membership_selection.py" to determine the preliminary spectroscopic sample
            #
            #
            exec(open('spec_membership_selection_file.py').read()) 
            #
            #
            #
            ## Call "phot_membership_selection.py" to determine the preliminary photometric sample
            #
            #
            exec(open('phot_membership_selection_file.py').read()) 
            #
            #
            #
            ## Call "data_mass_completeness*.py", to determine the limiting mass of each cluster for the redshift cuts you just adopted     #
            #
            exec(open('data_mass_completeness_5.py').read()) 
            #
            #
            ## Now call the "binning" program, "spec_completeness_binning.py", which takes z_cutoff[spec,phot] & limiting_mass as arguments. this file will evaluate different binning options for false pos/neg (above limiting mass of their respective cluster) in order to find the combination of z_cutoff[spec,phot] which yields correction factors clostest to 1 overall
            #
            #
            exec(open('spec_completeness_binning.py').read())      #opens and executes the script 'spec_completeness_binning.py' 
            #
            #
            ## now open the file just created, aggregate all such files into a single table, and sorts them according to metric of merit
            #
            filename = '/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/diagnostic_outputs/spec_binning/z_spec_%.3f'%z_cutoff[0]+'/binning_%.3f'%z_cutoff[0]+'_spec_cutoff_%.3f'%z_cutoff[1]+'_phot_cutoff.txt'
            df = pd.read_csv(filename, delimiter=',',index_col=None, header=0)
            df_list.append(df)
    ## save the list when ordered by z_cutoff
    df_result = pd.concat(df_list, axis=0, ignore_index=True)
    df_result = df_result.sort_values(by=['z_spec_cutoff','z_phot_cutoff','M.o.M.'])       
    #
    ## print the results to a file 
    filename ='/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/diagnostic_outputs/spec_binning/variational_analysis_results_%s'%increment[0]+'_sorted_by_ZCUT.txt'
    df_result.to_csv(filename,sep=' ',na_rep='NaN',index=False,header=True)  
    #
    ## save the list when ordered by z_cutoff
    df_result = pd.concat(df_list, axis=0, ignore_index=True)
    df_result = df_result.sort_values(by=['M.o.M.'])       
    #
    ## print the results to a file 
    filename ='/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/diagnostic_outputs/spec_binning/variational_analysis_results_%s'%increment[0]+'_sorted_by_VAR.txt'
    df_result.to_csv(filename,sep=' ',na_rep='NaN',index=False,header=True)
    #
    #
    ## TIME_FLAG END
    #
    if time_flag_4 == 1 or time_flag == 2:
        print('\n"master_data*.py" Section 4 diagnostic complete.\n\nSection 4 took: %s seconds to run variational analysis on spec/phot membership cutoffs.\n\n' % (time.time() - start_time))
    #
    print('PROGRAM SHOULD EXIT AFTER PRINTING THIS STATEMENT')
    sys.exit()
    print('PROGRAM SHOULD HAVE EXITED BEFORE PRINTING THIS STATEMENT')           
#
#
#
#if adams_flag == 0:
#    pass
#
#
## END of DIAGNOSTIC loop
#
else:
    # define cut-offs for SF & Q
    #z_cutoff = [0.02,0.06]#[0.05,0.05]         # NOTE: THIS IS NOW SET IN "main_project_file.py"
    #
    ## Call "spec_membership_selection.py" to determine spec. membership based on hard-coded cuts you just made
    #master_cat, mem_spec, field_spec, pos_spec, neg_spec, lost_due_to_buffer_spec = spec_membership_selection(master_cat,z_cutoff)
    #
    exec(open('spec_membership_selection_file.py').read()) 
    #
    #
    #
#    
#
## calculate membership acceptance fraction
mem_fraction_spec = np.array([0]*2,dtype='float32')
#
mem_fraction_spec[0] = np.round(np.sum(mem_spec[0])/np.sum([mem_spec[0],field_spec[0],pos_spec[0],neg_spec[0]]),decimals=3)
mem_fraction_spec[1] = np.round(np.sum(mem_spec[1])/np.sum([mem_spec[1],field_spec[1],pos_spec[1],neg_spec[1]]),decimals=3)
#
#    
## SECTION (4.1): SUMMARY table
##  summarize data population as segregated above, and display in a table
#
if summary_flag_4 == 1 or adams_flag == 1:
    ## Summarize initial data stats in table
    spec_member_names = Column(['SPEC subsample total','SF - member','SF - field','SF - fasle pos','SF - false neg','SF - LDTB','Q - member','Q - field','Q - false pos','Q - false neg','Q - LDTB','SUM'],name='Property')
    col_names = cluster_names
    # SF table
    spec_member0 = Column([(np.sum(both)-np.sum(outliers)),np.sum(mem_spec[0]),np.sum(field_spec[0]),np.sum(pos_spec[0]),np.sum(neg_spec[0]),np.sum(lost_due_to_buffer_spec[0]),np.sum(mem_spec[1]),np.sum(field_spec[1]),np.sum(pos_spec[1]),np.sum(neg_spec[1]),np.sum(lost_due_to_buffer_spec[1]),np.sum([np.sum(mem_spec),np.sum(field_spec),np.sum(pos_spec),np.sum(neg_spec),np.sum(lost_due_to_buffer_spec)])],name='Total')  # total column
    spec_member_stats = Table([spec_member_names,spec_member0])
    for ii in range(len(mem_spec[0])):
        col = Column([(both[ii]-outliers[ii]),mem_spec[0][ii],field_spec[0][ii],pos_spec[0][ii],neg_spec[0][ii],lost_due_to_buffer_spec[0][ii],mem_spec[1][ii],field_spec[1][ii],pos_spec[1][ii],neg_spec[1][ii],lost_due_to_buffer_spec[1][ii],np.sum([mem_spec[0][ii],field_spec[0][ii],pos_spec[0][ii],neg_spec[0][ii],lost_due_to_buffer_spec[0][ii],mem_spec[1][ii],field_spec[1][ii],pos_spec[1][ii],neg_spec[1][ii],lost_due_to_buffer_spec[1][ii]])],name=col_names[ii])
        spec_member_stats.add_column(col)  # add columns to table one cluster at a time
    #
    print('\nSummary Table 4 - SPEC Subsample\nCatalogue by MEMBER:')
    print(spec_member_stats)
    print("\nNOTE: LDTB = Lost Due To Buffer b/w member & field (recall: the def'n of 'field' doesn't start where the def'n of 'member' stops. There's a buffer in between them.\n")
    print('SF total: %s'%np.sum([mem_spec[0],field_spec[0],pos_spec[0],neg_spec[0],lost_due_to_buffer_spec[0]]),'\nQ total: %s'%np.sum([mem_spec[1],field_spec[1],pos_spec[1],neg_spec[1],lost_due_to_buffer_spec[1]]),'\nSF total + Q total + Outliers = %s'%np.sum([mem_spec[0],field_spec[0],pos_spec[0],neg_spec[0],lost_due_to_buffer_spec[0]]),' + %s'%np.sum([mem_spec[1],field_spec[1],pos_spec[1],neg_spec[1],lost_due_to_buffer_spec[1]]),' + %s'%np.sum(outliers),' = %s'%np.sum([np.sum(mem_spec),np.sum(field_spec),np.sum(pos_spec),np.sum(neg_spec),np.sum(lost_due_to_buffer_spec),np.sum(outliers)]))
    #####
#
#                       
#
## TIME_FLAG_4 END
#
if time_flag_4 == 1 and time_flag == 2:
    print('"master_data*.py" Section 4 took: %s seconds.\n\n' % (time.time() - start_time))
else:
    print('"master_data*.py" Section 4 complete.')
#
#
#
#
#
## SETION (5): make cuts to PHOTOMETRIC SUBSAMPLE based on defintion for del_z:  del_z = (z_phot - z_cl)/(1 + z_phot) < **some number** (to match cut made above in z_cutoff[1]; apply MEMBER FILTER for PHOTOMETRIC SUBSAMPLE MEMBERS); same photometric cut made to specroscopic sub-sample. this is a preliminary measure for determining the photometric sample, final corrections will be made by mass bins to match false pos/neg fractions in spec. sample per van der Burg (2013)
#
## apply cut at |z_clusterphot| < z_cutoff[1] to separate cluster members from field for targets with photometry only. store in MEMBER FILTER: 
## 0 = cluster member; 1 = field; 2 = false pos; 3 = false neg; 4 = field outlier
#
## recall: z_clusterphot defined as (z_peak - z_cluster / 1 + z_peak)
#
#
## TIME_FLAG_5 START
time_flag_5 = 1     # track & print time to execute current section
#
if time_flag_5 == 1 and time_flag == 2:
    start_time = time.time()
#
print('\n"master_data*.py" Section 5: photometric cluster membership selection beginning...')          
#
## Call the program that selects phot members based on criteria in Section (4)
exec(open('phot_membership_selection_file.py').read())            
#
#
## SECTION (5.1): SUMMARY table
##  summarize data population as segregated above, and display in a table
#
if summary_flag_5 == 1 or adams_flag == 1:
    ## Summarize initial data stats in table
    phot_member_names = Column(['PHOT subsample total','SF - Total','SF - member','SF - field','SF - outlier','SF - LDTB','Q - Total','Q - member','Q - field','Q - outlier','Q - LDTB','SUM (less totals)'],name='Property')
    col_names = cluster_names
    # SF table
    phot_member0 = Column([np.sum(phot_only),np.sum([np.sum(mem_phot[0]),np.sum(field_phot[0]),np.sum(lost_due_to_buffer_phot[0]),np.sum(field_outliers[0])]),np.sum(mem_phot[0]),np.sum(field_phot[0]),np.sum(field_outliers[0]),np.sum(lost_due_to_buffer_phot[0]),np.sum([np.sum(mem_phot[1]),np.sum(field_phot[1]),np.sum(lost_due_to_buffer_phot[1]),np.sum(field_outliers[1])]),np.sum(mem_phot[1]),np.sum(field_phot[1]),np.sum(field_outliers[1]),np.sum(lost_due_to_buffer_phot[1]),np.sum([np.sum(mem_phot),np.sum(field_phot),np.sum(field_outliers),np.sum(lost_due_to_buffer_phot)])],name='Total')  # total column
    phot_member_stats = Table([phot_member_names,phot_member0])
    for ii in range(len(mem_phot[0])):
        col = Column([phot_only[ii],np.sum([mem_phot[0][ii],field_phot[0][ii],field_outliers[0][ii],lost_due_to_buffer_phot[0][ii]]),mem_phot[0][ii],field_phot[0][ii],field_outliers[0][ii],lost_due_to_buffer_phot[0][ii],np.sum([mem_phot[1][ii],field_phot[1][ii],field_outliers[1][ii],lost_due_to_buffer_phot[1][ii]]),mem_phot[1][ii],field_phot[1][ii],field_outliers[1][ii],lost_due_to_buffer_phot[1][ii],np.sum([mem_phot[0][ii],mem_phot[1][ii],field_phot[0][ii],field_phot[1][ii],field_outliers[0][ii],field_outliers[1][ii],lost_due_to_buffer_phot[0][ii],lost_due_to_buffer_phot[1][ii]])],name=col_names[ii])
        phot_member_stats.add_column(col)  # add columns to table one cluster at a time
    #
    #
    print('\nSummary table 5 - PHOT-ONLY Subsample\nCatalogue by MEMBER - Star-forming:')
    print(phot_member_stats)
    print('\nNOTE: LDTB = Lost Due To Buffer b/w member & field.\nTotal Field Outliers (z>0.6 or z<0.25): %s' % np.sum(field_outliers))
#####
#
    mem_phot_fraction = np.array([0]*2,dtype='float32')     # to keep track of membership acceptance fraction
    mem_phot_fraction[0] = (np.sum(mem_phot[0]) / n_SF)
    mem_phot_fraction[1] = (np.sum(mem_phot[1]) / n_Q)
    print('\nOverall membership fraction: \nSF: %s'%mem_phot_fraction[0],' & Q: %s'%mem_phot_fraction[1],'   for cutoff: %s'%str(z_cutoff))
    #
    print('\nTotal catalogue length: %s'%len(master_cat))
    print('PHOT ONLY sub-sample: %s' %n_phot_only)
    print('SF (members + field): %s' % np.sum([mem_phot[0],field_phot[0]]))
    print('Q (members + field): %s' % np.sum([mem_phot[1],field_phot[1]]))
    print('\nField outliers (z>0.6 or z<0.3): %s' % np.sum(field_outliers))
    print('Lost due to buffer b/w definition of cluster member/field: %s'%np.sum(lost_due_to_buffer_phot))
    print('Stars & outliers: %s' % stars_outliers)
    print('Sum of the above: %s'%np.sum([np.sum(lost_due_to_buffer_phot),np.sum(field_outliers),np.sum([mem_phot[1],field_phot[1]]),np.sum([mem_phot[0],field_phot[0]])]))
    print('Difference between # in PHOT-ONLY subsample & sum above: %s'%(n_phot_only - np.sum([np.sum(lost_due_to_buffer_phot),np.sum(field_outliers),np.sum([mem_phot[1],field_phot[1]]),np.sum([mem_phot[0],field_phot[0]])])))
    print('Other (not in phot only subsample): %s'%other_phot)
#####                      
#
## TIME_FLAG_5 END
#
if time_flag_5 == 1 and time_flag == 2:
    print('"master_data*.py" Section 5 took: %s seconds.\n\n' % (time.time() - start_time))
else:
    print('"master_data*.py" Section 5 complete.')
#
## SECTION (5.2): Overall MEMBERSHIP summary table (spec + phot subsamples)
#
if summary_flag_6 == 1 or adams_flag == 1:
    ## Summarize initial data stats in table
    member_names = Column(['Parent total','PHOT member','PHOT field (ALL)','PHOT TOTAL','SPEC member','SPEC field','SPEC false pos','SPEC false neg','SPEC TOTAL','SUM (less totals)'],name='Property')
    col_names = cluster_names
    # SF table
    member0 = Column([np.sum([phot_only,both]),np.sum(mem_phot),np.sum([field_phot,field_outliers]),np.sum([mem_phot,field_phot,field_outliers]),np.sum([mem_spec[0],mem_spec[1]]),np.sum(field_spec),np.sum(pos_spec),np.sum(neg_spec),np.sum([mem_spec,field_spec,pos_spec,neg_spec]),np.sum([mem_phot,mem_spec,field_phot,field_outliers,field_spec,pos_spec,neg_spec])],name='Total')  # total column
    #member0 = Column([np.sum([phot_only,both]),np.sum(mem_phot),np.sum(mem_phot[0]),np.sum(mem_phot[1]),np.sum([field_phot,field_outliers]),np.sum([mem_phot,field_phot,field_outliers]),np.sum([mem_spec[0],mem_spec[1]]),np.sum(mem_spec[0]),np.sum(mem_spec[1]),np.sum(field_spec),np.sum(pos_spec),np.sum(neg_spec),np.sum([mem_spec,field_spec,pos_spec,neg_spec]),np.sum([mem_phot,mem_spec,field_phot,field_outliers,field_spec,pos_spec,neg_spec])],name='Total')  # total column
    member_stats = Table([member_names,member0])
    for ii in range(len(mem_phot[0])):
        col = Column([np.sum([phot_only[ii],both[ii]]),np.sum([mem_phot[0][ii],mem_phot[1][ii]]),np.sum([field_phot[0][ii],field_phot[1][ii],field_outliers[0][ii],field_outliers[1][ii]]),np.sum([mem_phot[0][ii],mem_phot[1][ii],field_phot[0][ii],field_phot[1][ii],field_outliers[0][ii],field_outliers[1][ii]]),np.sum([mem_spec[0][ii],mem_spec[1][ii]]),np.sum([field_spec[0][ii],field_spec[1][ii]]),np.sum([pos_spec[0][ii],pos_spec[1][ii]]),np.sum([neg_spec[0][ii],neg_spec[1][ii]]),np.sum([mem_spec[0][ii],mem_spec[1][ii],field_spec[0][ii],field_spec[1][ii],pos_spec[0][ii],pos_spec[1][ii],neg_spec[0][ii],neg_spec[1][ii]]),np.sum([mem_phot[0][ii],mem_phot[1][ii],field_phot[0][ii],field_phot[1][ii],field_outliers[0][ii],field_outliers[1][ii],mem_spec[0][ii],mem_spec[1][ii],field_spec[0][ii],field_spec[1][ii],pos_spec[0][ii],pos_spec[1][ii],neg_spec[0][ii],neg_spec[1][ii]])],name=col_names[ii])
        #col = Column([np.sum([phot_only[ii],both[ii]]),np.sum([mem_phot[0][ii],mem_phot[1][ii]]),mem_phot[0][ii],mem_phot[1][ii],np.sum([field_phot[0][ii],field_phot[1][ii],field_outliers[0][ii],field_outliers[1][ii]]),np.sum([mem_phot[0][ii],mem_phot[1][ii],field_phot[0][ii],field_phot[1][ii],field_outliers[0][ii],field_outliers[1][ii]]),np.sum([mem_spec[0][ii],mem_spec[1][ii]]),mem_spec[0][ii],mem_spec[1][ii],np.sum([field_spec[0][ii],field_spec[1][ii]]),np.sum([pos_spec[0][ii],pos_spec[1][ii]]),np.sum([neg_spec[0][ii],neg_spec[1][ii]]),np.sum([mem_spec[0][ii],mem_spec[1][ii],field_spec[0][ii],field_spec[1][ii],pos_spec[0][ii],pos_spec[1][ii],neg_spec[0][ii],neg_spec[1][ii]]),np.sum([mem_phot[0][ii],mem_phot[1][ii],field_phot[0][ii],field_phot[1][ii],field_outliers[0][ii],field_outliers[1][ii],mem_spec[0][ii],mem_spec[1][ii],field_spec[0][ii],field_spec[1][ii],pos_spec[0][ii],pos_spec[1][ii],neg_spec[0][ii],neg_spec[1][ii]])],name=col_names[ii])
        member_stats.add_column(col)  # add columns to table one cluster at a time
    #
    LDTB_total = np.sum([lost_due_to_buffer_spec,lost_due_to_buffer_phot])
    #
    print('\nSummary Table 6 - FULL Parent Sample\nCatalogue by MEMBER:')
    print(member_stats)
    print('\nNOTE: Lost Due To Buffer b/w member & field: %s'%LDTB_total,'.\nTotal galaxies considered: Members + Buffer + Outliers(spec) = %s'%np.sum([mem_phot,mem_spec,field_phot,field_outliers,field_spec,pos_spec,neg_spec]),' + %s'%LDTB_total,' + %s'%np.sum(outliers),' = %s'%np.sum([np.sum(mem_phot),np.sum(mem_spec),np.sum(field_phot),np.sum(field_outliers),np.sum(field_spec),np.sum(pos_spec),np.sum(neg_spec),np.sum(LDTB_total),np.sum(outliers)]),'\nNOTE: Total (phot) Field Outliers (z>0.6 or z<0.3): %s'%np.sum(field_outliers))
    print('\nTOTAL MEMBERS SELECTED: phot + spec = %s'%(np.sum([mem_phot[0],mem_phot[1]])),' + %s'%(np.sum([mem_spec[0],mem_spec[1]])),' = %s'%np.sum([mem_phot,mem_spec]))
    #
#
#
#                        
## SECTION (6): BCGs. Brightest Cluster Galaxies (BCGs) need to be taken out of the cluster sample as they are unique to overly dense environments and so lack a counterpart in the field against which to make a fair comparison. as such, we remove them from our sample before making the SMF
#
print('\n"master_data*.py" Section 6: removing BCGs beginning...')
#
#
## TIME_FLAG_6 START
time_flag_6 = 0     # track & print time to execute current section
#
if time_flag_6 == 1 and time_flag == 2:
    start_time = time.time()
#
## BCGs are identified in the "flag_{band}" column, flag_{ban} = 4 for BCGs. I will use band=F160W to identify BCGs
#

exec(open('bcg_hist_diagnostic.py').read())            
#
BCG_delz_means = np.array([0]*6,dtype='float32')
small_BCG_delz_means = np.array([0]*6,dtype='float32')
#
for BCG in range(len(num_BCG)):
    BCG_delz_means[BCG] = np.round(np.median(BCG_delz[BCG]),decimals=3)
    small_BCG_delz_means[BCG] = np.round(np.median(small_BCG_delz[BCG]),decimals=3)             
#
if summary_flag_7 == 1 or adams_flag == 1:
    ## Summarize bCG stats in table
    BCG_names = Column(['Total (>12)','Avg. delz_phot','SF_total','SF_spec','SF_phot','SF_sum','Q_total','Q_spec','Q_phot','Q_sum','SF+Q_sum'],name='Property')
    col_names = ['macs0416','macs1149','macs0717','abell370','abell1063','abell2744']
    BCG0 = Column([np.sum(num_BCG),np.round(np.mean(BCG_delz_means),decimals=3),np.sum(BCG_SF),np.sum(BCG_spec[0]),np.sum(BCG_phot[0]),np.sum([BCG_spec[0],BCG_phot[0]]),np.sum(BCG_Q),np.sum(BCG_spec[1]),np.sum(BCG_phot[1]),np.sum([BCG_spec[1],BCG_phot[1]]),np.sum([BCG_spec[0],BCG_phot[0],BCG_spec[1],BCG_phot[1]])],name='Total')  # total column
    BCG_stats = Table([BCG_names,BCG0])
    for ii in range(len(num_BCG)):
        BCG_col = Column([num_BCG[ii],BCG_delz_means[ii],BCG_SF[ii],BCG_spec[0][ii],BCG_phot[0][ii], np.sum([BCG_spec[0][ii],BCG_phot[0][ii]]),BCG_Q[ii],BCG_spec[1][ii],BCG_phot[1][ii], np.sum([BCG_spec[1][ii],BCG_phot[1][ii]]),np.sum([BCG_spec[0][ii],BCG_phot[0][ii],BCG_spec[1][ii],BCG_phot[1][ii]])],name=col_names[ii])               # cluster columns
        BCG_stats.add_column(BCG_col) 
    #
    #
    #
    print('\nSummary Table 7: bCG stats\n%s'%BCG_stats)
    #
    ## print summary of operation
    print('Total potential (M > 11) bCGs identified from catalogue: %s'%np.sum(num_BCG),'.\nBCGs identified by cluster: %s'%num_BCG)
#
#
## TIME_FLAG_6 END
time_flag_6 = 1     # track & print time to execute current section
#
if time_flag_6 == 1 and time_flag == 2:
    print('"master_data*.py" Section 6 took: %s seconds.\n\n' % (time.time() - start_time))
else:
    print('"master_data*.py" Section 6 complete.')
#
#
#
#
#
## TIME_FLAG END
#
if time_flag == 1 or project_time_flag == 1:
    print('Program "master_data*.py" took: %s seconds to run.\n\n' % (time.time() - start_time))
#
#
#
print('\n\n"master_data*.py"  terminated successfully.\n')
#
#
#                        
###################     PROGRAM END     ###################