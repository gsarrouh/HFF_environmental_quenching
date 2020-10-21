# Created on Wed Jun 24 10:23:11 2020
#
#### UPDATE 07/09/20 ####
## The previous versions of this file all used HFF DR v3.5, a preliminary data release from the summer of 2017. The paper Shipley et al. 2018 (which I reference extensively in my work) corresponds to the finished data release, v3.9. That final, most up-tp-date data release is what is implemented in the following work.
#
#
### WHAT THIS PROGRAM DOES:
### This script reads in all data for the Hubble Frontier Fields Parallel Fields images and prepares data for plotting and analysis. Much of the code and commenting is borrowed from "master_data*.py". This file assembles the field sample for the Stellar Mass Function (SMF), isolating all galaxies in the parallel fields in the redshift range 0.3 < z < 0.55 (i.e. the redshift range of the clusters in the HFF).
#
### Data is organized in a single catalogue ("master_cat_par"), and objects are identified through applying a series of "FILTERS" to designate key populations using a numerical designation for ease of writing code.
#
#
## FILTER 1 - 'cluster':  parallel field catalogues are designated as follows:
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
## 3: spectroscopy only
## 4: star
##
## FILTER 3 - 'type': :   identifies type of data each object has
## 0: star
## 1: star-forming (SF)
## 2: quiescent (Q)
## 3: outliers (defined as |del_z\ > 0.15)
#
## FILTER 4 - 'member': :   identifies type of data each object has
## NOTE: this designation is for objects with phot only (sub =2) and spec&phot (sub=1) only, as membership determinaion requires a good 'z_phot' estimate; "member=0" is reserved for galaxies with 0.3 < z < 0.55; "member=1" for all others
## 0: available
## 1: secure field sample member         <-- this comprises the sample of field of the SMF
## 2: redshift outside field sample range
## 3: redshift within range, but galaxy below limiting mass of parallel field
#
#
### Section summary:
#
### PROGRAM START
#
### (0)    import modules, define functions;
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
### (4)    classify FIELD SAMPLE members ( 0.3 < z < 0.55 )
### (4.1)  add DIAG_FLAG_5: summarize in table "SF_spec_stats" & "Q_spec_stats"
### (5)    determine LIMITING MASS for each parallel field; calls the program
###        "data_mass_completeness*_par.py"
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
## FLAGS !!!
## MAY NEED TO EDIT
## MASTER diagnostic flag:   0=all flags off;   1=all flags on;    2=individual flags may be turned on
diag_flag_master = 1
#
## section summary flags   (all are turned on if "adams_flag"==1)
## these flags print section summary tables for bookkeeping purposes
summary_flag_1 = 1          # S1.2: display diagnostic summary table, describes SUB-type filter
summary_flag_2 = 1          # S2: display outlier fractions & UVJ table
summary_flag_3 = 1          # S3: display TYPE filter summary table (SF/Q)
summary_flag_4 = 1          # S4: MEMBER-filter classification
summary_flag_5 = 1          # S4: MEMBER-filter classification, post-limiting mass calculation
#
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
#
if 'project_diagnostic_flag' in locals():
    pass
else:
    project_plot_flag = 2
    project_diagnostic_flag = 2
    project_master_variational_flag = 0
    z_cutoff = [0.012,0.055]     # [spec,phot] cutoffs for cluster membership
    z_cutoff_field = [0.08,0.15]
    limiting_mass_flag = 1
    bin_width = 0.4
    diagnostic_round_flag = 2
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
print('"master_parallel*.py" Section 1: import data beginning...')
#
#
## import catalogues into single table "master_cat"; separate objects for which there is no redshift data (both photo & spec) as "nodata"
#
##create table from EAZY output redshift ".zout" file
z_macs0416 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs0416par_catalogs/macs0416par_v3.9.zout',format='ascii')
z_macs1149 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs1149par_catalogs/macs1149par_v3.9.zout',format='ascii')
z_macs0717 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs0717par_catalogs/macs0717par_v3.9.zout',format='ascii')
z_abell370 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell370par_catalogs/abell370par_v3.9.zout',format='ascii')
z_abell1063 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell1063par_catalogs/abell1063par_v3.9.zout',format='ascii')
z_abell2744 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell2744par_catalogs/abell2744par_v3.9.zout',format='ascii')
#
#create table from FAST ".fout" file (contains mass estimates)
f_macs0416 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs0416par_catalogs/macs0416par_v3.9.fout',format='ascii')
f_macs1149 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs1149par_catalogs/macs1149par_v3.9.fout',format='ascii')
f_macs0717 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs0717par_catalogs/macs0717par_v3.9.fout',format='ascii')
f_abell370 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell370par_catalogs/abell370par_v3.9.fout',format='ascii')
f_abell1063 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell1063par_catalogs/abell1063par_v3.9.fout',format='ascii')
f_abell2744 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell2744par_catalogs/abell2744par_v3.9.fout',format='ascii')
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
cat_macs0416 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs0416par_catalogs/hffds_macs0416par_v3.9.cat',format='ascii')
cat_macs1149 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs1149par_catalogs/hffds_macs1149par_v3.9.cat',format='ascii')
cat_macs0717 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs0717par_catalogs/hffds_macs0717par_v3.9.cat',format='ascii')
cat_abell370 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell370par_catalogs/hffds_abell370par_v3.9.cat',format='ascii')
cat_abell1063 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell1063par_catalogs/hffds_abell1063par_v3.9.cat',format='ascii')
cat_abell2744 = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell2744par_catalogs/hffds_abell2744par_v3.9.cat',format='ascii')
#
##creat table for colours
#macs0416
F0416_u = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs0416par_catalogs/macs0416par_v3.9.153.rf',format='ascii')
F0416_v = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs0416par_catalogs/macs0416par_v3.9.155.rf',format='ascii')
F0416_j = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs0416par_catalogs/macs0416par_v3.9.161.rf',format='ascii')
#macs1149
F1149_u = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs1149par_catalogs/macs1149par_v3.9.153.rf',format='ascii')
F1149_v = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs1149par_catalogs/macs1149par_v3.9.155.rf',format='ascii')
F1149_j = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs1149par_catalogs/macs1149par_v3.9.161.rf',format='ascii')
#macs0717
F0717_u = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs0717par_catalogs/macs0717par_v3.9.153.rf',format='ascii')
F0717_v = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs0717par_catalogs/macs0717par_v3.9.155.rf',format='ascii')
F0717_j = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/macs0717par_catalogs/macs0717par_v3.9.161.rf',format='ascii')
#abell370
F370_u = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell370par_catalogs/abell370par_v3.9.153.rf',format='ascii')
F370_v = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell370par_catalogs/abell370par_v3.9.155.rf',format='ascii')
F370_j = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell370par_catalogs/abell370par_v3.9.161.rf',format='ascii')
#abell1063
F1063_u = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell1063par_catalogs/abell1063par_v3.9.153.rf',format='ascii')
F1063_v = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell1063par_catalogs/abell1063par_v3.9.155.rf',format='ascii')
F1063_j = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell1063par_catalogs/abell1063par_v3.9.161.rf',format='ascii')
#abell2744
F2744_u = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell2744par_catalogs/abell2744par_v3.9.153.rf',format='ascii')
F2744_v = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell2744par_catalogs/abell2744par_v3.9.155.rf',format='ascii')
F2744_j = Table.read('/Users/gsarrouh/Research/NSERC_2017_HFF/nserc17/working_data/HFF_catalgoues/abell2744par_catalogs/abell2744par_v3.9.161.rf',format='ascii')
##aggregate into a single table
macs0416 = Table([z_macs0416['id'],z_macs0416['z_peak'],z_macs0416['z_spec'],F0416_u['L153'],F0416_v['L155'],F0416_j['L161'],F0416_u['DM'],f_macs0416['lmass'],f_macs0416['lsfr'],f_macs0416['lssfr'],cat_macs0416['flux_radius'],cat_macs0416['star_flag'],cat_macs0416['use_phot'],cat_macs0416['f_F160W'],cat_macs0416['e_F160W'],cat_macs0416['flag_F160W'],cat_macs0416['f_F814W'],cat_macs0416['flag_F814W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','use_phot','f_F160W','e_F160W','flag_F160W','f_F814W','flag_F814W'))
macs1149 = Table([z_macs1149['id'],z_macs1149['z_peak'],z_macs1149['z_spec'],F1149_u['L153'],F1149_v['L155'],F1149_j['L161'],F1149_u['DM'],f_macs1149['lmass'],f_macs1149['lsfr'],f_macs1149['lssfr'],cat_macs1149['flux_radius'],cat_macs1149['star_flag'],cat_macs1149['use_phot'],cat_macs1149['f_F160W'],cat_macs1149['e_F160W'],cat_macs1149['flag_F160W'],cat_macs1149['f_F814W'],cat_macs1149['flag_F814W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','use_phot','f_F160W','e_F160W','flag_F160W','f_F814W','flag_F814W'))
macs0717 = Table([z_macs0717['id'],z_macs0717['z_peak'],z_macs0717['z_spec'],F0717_u['L153'],F0717_v['L155'],F0717_j['L161'],F0717_u['DM'],f_macs0717['lmass'],f_macs0717['lsfr'],f_macs0717['lssfr'],cat_macs0717['flux_radius'],cat_macs0717['star_flag'],cat_macs0717['use_phot'],cat_macs0717['f_F160W'],cat_macs0717['e_F160W'],cat_macs0717['flag_F160W'],cat_macs0717['f_F814W'],cat_macs0717['flag_F814W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','use_phot','f_F160W','e_F160W','flag_F160W','f_F814W','flag_F814W'))
abell370 = Table([z_abell370['id'],z_abell370['z_peak'],z_abell370['z_spec'],F370_u['L153'],F370_v['L155'],F370_j['L161'],F370_u['DM'],f_abell370['lmass'],f_abell370['lsfr'],f_abell370['lssfr'],cat_abell370['flux_radius'],cat_abell370['star_flag'],cat_abell370['use_phot'],cat_abell370['f_F160W'],cat_abell370['e_F160W'],cat_abell370['flag_F160W'],cat_abell370['f_F814W'],cat_abell370['flag_F814W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','use_phot','f_F160W','e_F160W','flag_F160W','f_F814W','flag_F814W'))
abell1063 = Table([z_abell1063['id'],z_abell1063['z_peak'],z_abell1063['z_spec'],F1063_u['L153'],F1063_v['L155'],F1063_j['L161'],F1063_u['DM'],f_abell1063['lmass'],f_abell1063['lsfr'],f_abell1063['lssfr'],cat_abell1063['flux_radius'],cat_abell1063['star_flag'],cat_abell1063['use_phot'],cat_abell1063['f_F160W'],cat_abell1063['e_F160W'],cat_abell1063['flag_F160W'],cat_abell1063['f_F814W'],cat_abell1063['flag_F814W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','use_phot','f_F160W','e_F160W','flag_F160W','f_F814W','flag_F814W'))
abell2744 = Table([z_abell2744['id'],z_abell2744['z_peak'],z_abell2744['z_spec'],F2744_u['L153'],F2744_v['L155'],F2744_j['L161'],F2744_u['DM'],f_abell2744['lmass'],f_abell2744['lsfr'],f_abell2744['lssfr'],cat_abell2744['flux_radius'],cat_abell2744['star_flag'],cat_abell2744['use_phot'],cat_abell2744['f_F160W'],cat_abell2744['e_F160W'],cat_abell2744['flag_F160W'],cat_abell2744['f_F814W'],cat_abell2744['flag_F814W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','use_phot','f_F160W','e_F160W','flag_F160W','f_F814W','flag_F814W'))
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
master_cat_par = Table(np.concatenate((macs0416,macs1149,macs0717,abell370,abell1063,abell2744), axis=0))  #create a master catalogue of all clusters
cluster_names_par = ['M0416','M1149','M0717','A370','A1063','A2744']
#
## add "empty" columns (with value set = 99) for the remaining sieves: sub, type, member for [phot or spec],[SF or Q],[cluster member, field, false pos/neg] respectively
E1 = Column([-99]*len(master_cat_par), name='sub', dtype=np.int8)    # create columns
E2 = Column([-99]*len(master_cat_par), name='type', dtype=np.int8)
E3 = Column([-99]*len(master_cat_par), name='member', dtype=np.int8)
master_cat_par.add_columns([E1,E2,E3],[-1,-1,-1])                   # add columns to the end of table
#
#
##
## SECTION (1.1) - SUB-TYPE FILTER: this section classifies all objects as either:
##  sub = 0: no data (neither spec. nor phot.)
##  sub = 1: both
##  sub = 2: phot only
##  sub = 3: spec only
##  sub = 4: stars
#
# sift for targets without no data, spec only, phot only, and both
spec_only_par = np.array([0]*6)    # to keep track by cluster
phot_only_par = np.array([0]*6)
use_phot_par = np.array([[0]*6]*2)      # row1=good (use_phot = 1);  row2=bad (use_phot !=1)
error_par = 0
other_par = 0
both_par = np.array([0]*6)
no_data_par = np.array([0]*6)
stars_sub_par = np.array([0]*6)
#
for counter in range(len(master_cat_par)):
    if master_cat_par['z_spec'][counter] > 0 and master_cat_par['z_peak'][counter] > 0:  #entries w/ both spectroscopy and photometry  (SPECTROSCOPIC sub-sample)
        if master_cat_par['use_phot'][counter] ==0:     # use_phot = 0 means bad photometry
            for jj in range(len(spec_only_par)):
                if master_cat_par['cluster'][counter] == (jj+1):      # identify # of spec only objects by cluster
                    #use_phot[1][jj]+=1
                    if master_cat_par['star_flag'][counter] == 1:
                        master_cat_par['sub'][counter] = 4          # APPLY FILTER: sub=4 for STARS
                        stars_sub_par[jj]+=1
                    else:
                        spec_only_par[(jj)]+=1
                        master_cat_par['sub'][counter] = 3          # APPLY FILTER: sub=3 for objects w/ SPEC ONLY
        elif master_cat_par['use_phot'][counter] ==1:     # use_phot = 1 means good photometry
            for jj in range(len(both_par)):
                if master_cat_par['cluster'][counter] == (jj+1):      # identify # of (spec & phot) objects by cluster
                    #use_phot[0][jj]+=1
                    if master_cat_par['star_flag'][counter] == 1:
                        master_cat_par['sub'][counter] = 4          # APPLY FILTER: sub=4 for STARS; THIS SHOULDNT DO ANYTHING GIVEN THE "use_phot==1" DESIGNATION
                        stars_sub_par[jj]+=1
                    else:
                        both_par[(jj)]+=1
                        master_cat_par['sub'][counter] = 1          # APPLY FILTER: sub=1 for objects w/ BOTH SPEC & PHOT
        else: error_par+=1                              # just to keep track of erroneous objects
    elif master_cat_par['z_spec'][counter] > 0 and master_cat_par['z_peak'][counter] < 0:  #entries w/ spectroscopy alone
        for jj in range(len(spec_only_par)):
            if master_cat_par['cluster'][counter] == (jj+1):      # identify # of spec only objects by cluster
                if master_cat_par['star_flag'][counter] == 1:
                    master_cat_par['sub'][counter] = 4          # APPLY FILTER: sub=4 for STARS
                    stars_sub_par[jj]+=1
                else:
                    spec_only_par[(jj)]+=1
                    master_cat_par['sub'][counter] = 3          # APPLY FILTER: sub=3 for objects w/ SPEC ONLY
    elif master_cat_par['z_spec'][counter] < 0 and master_cat_par['z_peak'][counter] > 0:  #entries w/ photometry alone (PHOTOMETRIC sub-sample)
        if master_cat_par['use_phot'][counter] ==0:
            for jj in range(len(both_par)):
                if master_cat_par['cluster'][counter] == (jj+1):      # identify # of NO DATA objects by cluster
                    #use_phot[1][jj]+=1
                    if master_cat_par['star_flag'][counter] == 1:
                        master_cat_par['sub'][counter] = 4          # APPLY FILTER: sub=4 for STARS
                        stars_sub_par[jj]+=1
                    else:
                        no_data_par[(jj)]+=1
                        master_cat_par['sub'][counter] = 0          # APPLY FILTER: sub=0 for objects w/ NEITHER SPEC NOR PHOT
        elif master_cat_par['use_phot'][counter] ==1:
            for jj in range(len(both_par)):
                if master_cat_par['cluster'][counter] == (jj+1):      # identify # of phot only objects by cluster
                    #use_phot[0][jj]+=1
                    if master_cat_par['star_flag'][counter] == 1:
                        master_cat_par['sub'][counter] = 4          # APPLY FILTER: sub=4 for STARS
                        stars_sub_par[jj]+=1
                    else:
                        phot_only_par[(jj)]+=1
                        master_cat_par['sub'][counter] = 2          # APPLY FILTER: sub=2 for objects w/ PHOT ONLY
    elif master_cat_par['z_spec'][counter] < 0 and master_cat_par['z_peak'][counter] < 0:  #entries w/ no z estimates at all
        for jj in range(len(both_par)):
            if master_cat_par['cluster'][counter] == (jj+1):      # identify # of NO DATA objects by cluster
                if master_cat_par['star_flag'][counter] == 1:
                    master_cat_par['sub'][counter] = 4          # APPLY FILTER: sub=4 for STARS
                    stars_sub_par[jj]+=1
                else:
                    no_data_par[(jj)]+=1
                    master_cat_par['sub'][counter] = 0          # APPLY FILTER: sub=0 for objects w/ NEITHER SPEC NOR PHOT
    if master_cat_par['use_phot'][counter] == 1:
        for jj in range(len(both_par)):
            if master_cat_par['cluster'][counter] == (jj+1):
                use_phot_par[0][jj]+=1
    else:
        for jj in range(len(both_par)):
            if master_cat_par['cluster'][counter] == (jj+1):
                use_phot_par[1][jj]+=1
#
#
## SECTION (1.2): SUMMARY table
#
if summary_flag_1 == 1 or adams_flag == 1:
        ## Summarize initial data stats in table
        sub_par_names = Column(['FULL PHOT (Parent)','Total','spec & phot','only phot','spec only','no data','stars','SUM'],name='Property')
        col_names = cluster_names_par
        sub_par0 = Column([np.sum([phot_only_par,both_par]),np.sum([np.sum(both_par),np.sum(phot_only_par),np.sum(spec_only_par),np.sum(no_data_par),np.sum(stars_sub_par)]),np.sum(both_par),np.sum(phot_only_par),np.sum(spec_only_par),np.sum(no_data_par),np.sum(stars_sub_par),np.sum([both_par,phot_only_par,spec_only_par,no_data_par,stars_sub_par])],name='Total')  # total column
        sub_par_stats = Table([sub_par_names,sub_par0])
        for ii in range(len(spec_only_par)):
            sub_par_col = Column([np.sum([phot_only_par[ii],both_par[ii]]),np.sum([both_par[ii],phot_only_par[ii],spec_only_par[ii],no_data_par[ii],stars_sub_par[ii]]),both_par[ii],phot_only_par[ii],spec_only_par[ii],no_data_par[ii],stars_sub_par[ii],np.sum([both_par[ii],phot_only_par[ii],spec_only_par[ii],no_data_par[ii],stars_sub_par[ii]])],name=col_names[ii])  # add columns to table one cluster at a time
            sub_par_stats.add_column(sub_par_col)
        #
        print('\nSummary Table 1 - Catalogue by SUB-type (PAR):')
        print(sub_par_stats)
        print('\nOther skipped objects: %s'%other_par,'\nNOTE: "use_phot==1": %s'%np.sum(use_phot_par[0]),';  "use_phot==0": %s'%np.sum(use_phot_par[1]))
        print('NOTE: phot only + (spec+phot) samples are the FULL PHOT (Parent) sample w/: %s'%np.sum(phot_only_par+both_par),'\nNOTE: Difference b/w Parent sample & "use_phot==1": %s'%np.abs((np.sum(phot_only_par+both_par)-np.sum(use_phot_par[0]))))
#
#
## SECTION (1.3): convert FLUX TO MAGNITUDE; using well-known mag = -2.5*log_10(flux) + zero_pt. zero_pt = 25
#
# add columns for luminosity calculations
empty_u = Column([99]*len(master_cat_par), name='L_u', dtype=np.float64)
empty_v = Column([99]*len(master_cat_par), name='L_v', dtype=np.float64)
empty_j = Column([99]*len(master_cat_par), name='L_j', dtype=np.float64)
empty_uv = Column([99]*len(master_cat_par), name='uv', dtype=np.float64)
empty_vj = Column([99]*len(master_cat_par), name='vj', dtype=np.float64)
master_cat_par.add_columns([empty_u,empty_v,empty_j,empty_uv,empty_vj])
#
## convert flux to magnitude (erroneously labelled as luminosity, e.g. L_u for magnitude in UV), get color indices U-V, V-J, add to table.
## the number of objects this calculation is performed on should be equal to the number of (phot_only + (phot+spec)) objects identified above
#
objects_99_par = np.array([0]*6)
use_phot_par = np.array([[0]*6]*2)      # row1=good (use_phot = 1);  row2=bad (use_phot !=1)
skipped_UVJ_par = np.array([0]*6)    # to track stars & objects w/o photometry
good_objects_par = np.array([0]*6)
#
for counter in range(len(master_cat_par)):
    if master_cat_par['use_phot'][counter] == 1:
        for jj in range(len(skipped_UVJ_par)):
            if master_cat_par['cluster'][counter] == (jj+1):      # identify # of NO DATA objects by cluster
                use_phot_par[0][jj]+=1
                if master_cat_par['sub'][counter] ==0 or master_cat_par['sub'][counter] ==3 or master_cat_par['sub'][counter] ==4:      # skip objects w/ "no data" (sub=0); "spec only" (sub=3); "stars" (sub=4)
                    skipped_UVJ_par[jj]+=1     # objects not in our FULL PHOT sample
                elif master_cat_par['sub'][counter] == -99:
                    objects_99_par[jj]+=1      # error objects that were not classified in section 1, if they even exist?
                elif master_cat_par['sub'][counter] == 1 or master_cat_par['sub'][counter] == 2:   # this should select all "use_phot=1" objects
                    master_cat_par['L_u'][counter] = -2.5*np.log10(master_cat_par['u'][counter]) + 25
                    master_cat_par['L_v'][counter] = -2.5*np.log10(master_cat_par['v'][counter]) + 25
                    master_cat_par['L_j'][counter] = -2.5*np.log10(master_cat_par['j'][counter]) + 25
                    master_cat_par['uv'][counter] = master_cat_par['L_u'][counter] - master_cat_par['L_v'][counter]
                    master_cat_par['vj'][counter] = master_cat_par['L_v'][counter] - master_cat_par['L_j'][counter]
                    good_objects_par[jj]+=1

    else:
        for jj in range(len(use_phot_par[1])):
            if master_cat_par['cluster'][counter] == (jj+1):
                use_phot_par[1][jj]+=1
#################
#
#
if summary_flag_2 == 1 or adams_flag == 1:
        ## Summary table of U V J and excesses calc., total numbers should agree w/ those reported in "sub_stats"
        UVJ_par_names = Column(['Parent sample','Objects w/ phot','Objects w/o phot','Objects not classified','SUM'],name='Property')
        col_names = cluster_names_par
        UVJ_par0 = Column([np.sum([phot_only_par,both_par]), np.sum(good_objects_par),np.sum(skipped_UVJ_par),np.sum(objects_99_par),np.sum([good_objects_par,skipped_UVJ_par,objects_99_par])],name='Total')  # total column
        UVJ_par_stats = Table([UVJ_par_names,UVJ_par0])
        for ii in range(len(spec_only_par)):
            UVJ_par_col = Column([np.sum([phot_only_par[ii],both_par[ii]]),good_objects_par[ii],skipped_UVJ_par[ii],objects_99_par[ii],np.sum([good_objects_par[ii],skipped_UVJ_par[ii],objects_99_par[ii]])],name=col_names[ii])  # add columns to table one cluster at a time
            UVJ_par_stats.add_column(UVJ_par_col)
        #####
        print('\nUVJ and excesses calculation (identified w/ "use_phot==1"): \n%s'%UVJ_par_stats,'\n\nNOTE: "objects w/o phot" are those w/ "use_flag" == 1 that have been classified as no good data/bad phot.\nNote: "use_phot" = 1: %s'%np.sum(use_phot_par[0]),'; "use_phot" = 0: %s'%np.sum(use_phot_par[1]))
#########
#
print('"master_parallel*.py" Section 1 complete.\n')
#
#
#
#
#
## SECTION (2): calulate DEL_Z's & separate OUTLIERS
#
print("\n'master_parallel*.py' Section 2: calculating del_z's and identifying outliers...\n")
#
#
#####Note: the photometric redshift used is 'z_peak' column from data
#
#(i) calculate delta_z for targets w/ both z_spec & z_phot
master_cat_par.add_column(Column([-99]*len(master_cat_par),name='del_z', dtype=np.float64))           # del_z = (z_phot - z_spec) / (1 + z_spec)
master_cat_par.add_column(Column([-99]*len(master_cat_par),name='z_clusterspec', dtype=np.float64))   # del_z = (z_spec - z_cl) / (1 + z_spec)
master_cat_par.add_column(Column([-99]*len(master_cat_par),name='z_clusterphot', dtype=np.float64))   # del_z = (z_phot - z_cl) / (1 + z_phot)
#
# store cluster redshifts; obtained from https://archive.stsci.edu/prepds/frontier/
z_cluster = [0.396,0.543,0.545,0.375,0.348,0.308]
# cluster_names = ['M0416','M1149','M0717','A370','A1063','A2744']
#
## calucalte del_z, z_clusterspec, z_clusterphot for outlier cut (defined above); these will be used to make cuts (member, field, false pos/neg) to spec sample, from which we will use relative fractions by mass to correct the photometric sample for completeness.
#
#
spec_subsample_par = np.array([0]*6)   # to track objects w/ BOTH (spec + phot), by cluster
phot_subsample_par = np.array([0]*6)   # to track objects w/ ONLY PHOT
skipped_delz_par = np.array([0]*6)     # to track all other objects
#
for counter in range(len(master_cat_par)):
    if master_cat_par['sub'][counter] == 1:   # sub=1 identifies spec&phot subsample
        master_cat_par['del_z'][counter] = ((master_cat_par['z_peak'][counter] - master_cat_par['z_spec'][counter]) / (1 + master_cat_par['z_spec'][counter]))
        #for ii in range(len(z_cluster)):
        #    if master_cat['cluster'][counter] == (ii+1):
        #        spec_subsample[ii]+=1
        #        master_cat['z_clusterspec'][counter] = ((master_cat['z_spec'][counter] - z_cluster[ii]) / (1 + master_cat['z_spec'][counter]))
        #        master_cat['z_clusterphot'][counter] = ((master_cat['z_peak'][counter] - z_cluster[ii]) / (1 + master_cat['z_peak'][counter]))
    #elif master_cat['sub'][counter] == 2:   # sub=2 identifies phot-only subsample
    #    for ii in range(len(z_cluster)):
    #        if master_cat['cluster'][counter] == (ii+1):
    #            phot_subsample[ii]+=1
    #            master_cat['z_clusterphot'][counter] = ((master_cat['z_peak'][counter] - z_cluster[ii]) / (1 + master_cat['z_peak'][counter]))
    else:
        for ii in range(len(z_cluster)):
            if master_cat_par['cluster'][counter] == (ii+1):
                skipped_delz_par[ii]+=1
                #
            #####
        #####
#########
#
#
## SECTION (2.1): separate OUTLIERS from both phot & spec sub-sample, defined as |del_z| < 0.15. apply FILTER TYPE = 3 for outliers;
#
outliers_par = np.array([0]*6)      # initialize array to track outliers by cluster, for computing outlier fraction later
sum_delz_par = []                   # for computing mean |del_z|
count_stars_par = 0
#
for counter in range(len(master_cat_par)):
    if master_cat_par['sub'][counter] == 1:                        # sub=1 for objects with both spec & phot; total # of such objects tracked by cluster above in the array called "both"
        if np.abs(master_cat_par['del_z'][counter]) > 0.15:        # |del_z| > 0.15 for outliers, threshold chosen for historical reasons to facilitate comparison with other studies
            master_cat_par['type'][counter] = 3                  # type=3 identifies outliers
            for ii in range(len(outliers_par)):
                if master_cat_par['cluster'][counter] == (ii+1):   # keep track of outliers by cluster
                    outliers_par[ii]+=1
            #sum_delz.append(np.abs(master_cat['del_z'][counter]))    # keep track of all |del_z| measurements for stats computation
        else:
            sum_delz_par.append(np.abs(master_cat_par['del_z'][counter]))
        if master_cat_par['type'][counter] == 3 and master_cat_par['star_flag'][counter] == 1:                  # overwrite designation for stars
            master_cat_par['type'][counter] = 0                        # type=0 for stars
            outliers_par-=1
            count_stars_par+=1
#
#
## SECTION (2.2): compute & DISPLAY OUTLIER FRACTION, SCATTER (i.e. std dev), and MEAN of |del_z|.
#
if summary_flag_2 == 1 or adams_flag == 1:
    #
    #####
    print('\nSummary Table 2: delz calculation (PAR): \n')
    print('\nNOTE: Total good objects: %s'%len(sum_delz_par))
    ###
    delz_median_par = np.median(sum_delz_par)
    delz_scatter_par = np.std(sum_delz_par)
    print('\nOUTLIERS total (in spec subsample): %s' % np.sum(outliers_par),'\nOutlier fraction: %s' % (np.sum(outliers_par)/np.sum(both_par)))
    print('|del_z| median: %s'%delz_median_par)
    print('|del_z| scatter: %s\n'%delz_scatter_par)
    print('NOTE: the new "use-able" number for the spec subsample is: %s'%(np.sum(both_par)-np.sum(outliers_par)))
#
#
print('"master_parallel*.py" Section 2 complete.\n')
#
#
#
#
#
## SECTION (3): add FILTER TYPE: separate SF/Q for both subsamples;    filter name: 'type'  IDs below; SELECTION CRITERIA from Shipley et al. 2018, Section (5.3)
##   0 = stars;  1 = SF (star-forming);    2 = Q (quiscient);    3 = outliers
#
print('\n"master_parallel*.py" Section 3: classifying galaxy TYPE as star-forming or quiescent...')
#
#
SF_type_par = np.array([[0]*6]*2)                         # initialize arrays;  row1=spec;  row2=phot
Q_type_par = np.array([[0]*6]*2)
stars_type_par = np.array([0]*6)                          # count # of stars by cluster
lost_type_par = np.array([0]*6)                           # objects lost due no data (sub=0), spec only (sub=3),
stars_outliers_par = np.array([0]*6)
#
for counter in range(len(master_cat_par)):
    if master_cat_par['star_flag'][counter] == 1:          # identify STARS, filter type=0
        master_cat_par['type'][counter] = 0
        for ii in range(len(stars_type_par)):
            if master_cat_par['cluster'][counter] == (ii+1):   # keep stars of outliers by cluster
                stars_type_par[ii]+=1
    elif master_cat_par['sub'][counter]==1 or master_cat_par['sub'][counter]==2:    # sub=1 for (spec & phot) subsample, sub=2 for phot subsample; i.e. look at all objects with photometry   <-- THIS IS THE LINE THAT MATTERS!
        if master_cat_par['vj'][counter] < 0.75:
            if master_cat_par['uv'][counter] < 1.3:
                if master_cat_par['type'][counter] ==3 or master_cat_par['type'][counter] ==0:        # skip outliers & stars
                    for ii in range(len(stars_outliers_par)):
                        if master_cat_par['cluster'][counter] == (ii+1):   # keep stars of outliers by cluster
                            if master_cat_par['sub'][counter] == 1:
                                stars_outliers_par[ii]+=1
                else:
                    master_cat_par['type'][counter] = 1             # identify STAR-FORMING galaxies, type=1
                    for ii in range(len(SF_type_par[0])):
                        if master_cat_par['cluster'][counter] == (ii+1):   # keep stars of outliers by cluster
                            if master_cat_par['sub'][counter] == 1:
                                SF_type_par[0][ii]+=1               # track spec
                            else: SF_type_par[1][ii]+=1             # track phot
            else:
                if master_cat_par['type'][counter] ==3 or master_cat_par['type'][counter] ==0:        # skip outliers & stars
                    for ii in range(len(stars_outliers_par)):
                        if master_cat_par['cluster'][counter] == (ii+1):   # keep stars of outliers by cluster
                            if master_cat_par['sub'][counter] == 1:
                                stars_outliers_par[ii]+=1
                else:
                    master_cat_par['type'][counter] = 2             # identify passive (QUIESCENT) galaxies, type=2
                    for ii in range(len(Q_type_par[0])):
                        if master_cat_par['cluster'][counter] == (ii+1):   # keep
                            if master_cat_par['sub'][counter] == 1:
                                Q_type_par[0][ii]+=1               # track spec
                            else:
                                Q_type_par[1][ii]+=1               # track phot
        elif master_cat_par['vj'][counter] >= 0.75:
            if master_cat_par['uv'][counter] < ( (0.8 * master_cat_par['vj'][counter]) + 0.7 ):
                if master_cat_par['type'][counter] ==3 or master_cat_par['type'][counter] ==0:        # skip outliers & stars
                    for ii in range(len(stars_outliers_par)):
                        if master_cat_par['cluster'][counter] == (ii+1):   # keep stars of outliers by cluster
                            if master_cat_par['sub'][counter] == 1:
                                stars_outliers_par[ii]+=1
                else:
                    master_cat_par['type'][counter] = 1             # identify STAR-FORMING galaxies, type=1
                    for ii in range(len(Q_type_par[0])):
                        if master_cat_par['cluster'][counter] == (ii+1):   # keep stars of outliers by cluster
                            if master_cat_par['sub'][counter] == 1:
                                SF_type_par[0][ii]+=1
                            else: SF_type_par[1][ii]+=1
            else:
                if master_cat_par['type'][counter] ==3 or master_cat_par['type'][counter] ==0:        # skip outliers & stars
                    for ii in range(len(stars_outliers_par)):
                        if master_cat_par['cluster'][counter] == (ii+1):   # keep stars of outliers by cluster
                            if master_cat_par['sub'][counter] == 1:
                                stars_outliers_par[ii]+=1
                else:
                    master_cat_par['type'][counter] = 2             # identify passive (QUIESCENT) galaxies, type=2
                    for ii in range(len(SF_type_par[0])):
                        if master_cat_par['cluster'][counter] == (ii+1):   # keep
                            if master_cat_par['sub'][counter] == 1:
                                Q_type_par[0][ii]+=1               # track spec vs phot
                            else:
                                Q_type_par[1][ii]+=1
    else:
        for ii in range(len(lost_type_par)):
            if master_cat_par['cluster'][counter] == (ii+1):   # keep stars of outliers by cluster
                lost_type_par[ii]+=1      # objects lost due to spec only (sub=3),no data (sub=0), and possibly outliers (type=3)
#
#
## SECTION (3.1): SUMMARY table
##  summarize data TYPE population as segregated above, and display in a table
#
if summary_flag_3 == 1 or adams_flag == 1:
    ## Summarize initial data stats in table
    type_par_names = Column(['Parent sample','SF - total','SF - phot','SF - spec','Q - total','Q - phot','Q - spec','Stars & Outliers','SUM (less totals)'],name='Property')
    col_names = cluster_names_par
    type_par0 = Column([np.sum([phot_only_par,both_par]),np.sum(SF_type_par),np.sum(SF_type_par[1]),np.sum(SF_type_par[0]),np.sum(Q_type_par),np.sum(Q_type_par[1]),np.sum(Q_type_par[0]),np.sum(stars_outliers_par),np.sum([np.sum(SF_type_par),np.sum(Q_type_par),np.sum(stars_outliers_par)])],name='Total')  # total column
    type_par_stats = Table([type_par_names,type_par0])
    for ii in range(len(spec_only_par)):
        type_par_col = Column([np.sum([phot_only_par[ii],both_par[ii]]),(SF_type_par[1][ii]+SF_type_par[0][ii]),SF_type_par[1][ii],SF_type_par[0][ii],(Q_type_par[1][ii]+Q_type_par[0][ii]),Q_type_par[1][ii],Q_type_par[0][ii],stars_outliers_par[ii],np.sum([SF_type_par[1][ii],SF_type_par[0][ii],Q_type_par[1][ii],Q_type_par[0][ii],stars_outliers_par[ii]])],name=col_names[ii])
        type_par_stats.add_column(type_par_col)  # add columns to table one cluster at a time
    #
    print('\nSummary Table 3 - Catalogue by TYPE (PAR): %s'%type_par_stats)
    #
#
print('\n"master_parallel*.py" Section 3 complete.')
#
#
#
#
#
## SECTION (4) : apply MEMBER FILTER: 'member=1' identifies our field sample for the SMF
#
## compute redshift range of galaxies in cluster sample
lower_bound = (min(z_cluster) - z_cutoff[1]) / (1 + z_cutoff[1])
upper_bound = (max(z_cluster) + z_cutoff[1]) / (1 - z_cutoff[1])
z_field_bounds = [lower_bound, upper_bound]
#
a = np.array([0]*6)    # counting array to track total number of objects in each cluster
count_field_sample = 0
count_field_sample_type = np.array([[0]*6]*2)      # row1=SF;  row2=Q, by cluster
count_not_in_field_sample = np.array([0]*6)
count_not_in_parent_sample = np.array([0]*6)
outliers_par = np.array([0]*6)
#
## isolate all galaxies (SF & Q) in the redshift range 0.3 < z < 0.55, for the field sample of the SMF
for counter in range(len(master_cat_par)):
    for ii in range(len(count_field_sample_type[0])):
        if master_cat_par['cluster'][counter] == (ii+1):
            a[ii]+=1
    if master_cat_par['sub'][counter] == 1 or master_cat_par['sub'][counter] == 2:  # only look at Parent sample (spec+phot & phot_only)
        if master_cat_par['z_peak'][counter] < z_field_bounds[1] and master_cat_par['z_peak'][counter] > z_field_bounds[0]:
            master_cat_par['member'][counter] = 1                                   # member=1: preliminary member of FIELD SAMPLE
            count_field_sample+=1
            for ii in range(len(count_field_sample_type[0])):
                if master_cat_par['cluster'][counter] == (ii+1):                    # by cluster
                    if master_cat_par['type'][counter] == 1:                        # SF
                        count_field_sample_type[0][ii]+=1
                    elif master_cat_par['type'][counter] == 2:                      # Q
                        count_field_sample_type[1][ii]+=1
                    else:
                        outliers_par[ii]+=1
                    #####
            ########
        ####
    #####
        elif master_cat_par['z_peak'][counter] > z_field_bounds[1] or master_cat_par['z_peak'][counter] < z_field_bounds[0]:
            master_cat_par['member'][counter] = 2                                   # member=1: redshift outside of field sample
            for ii in range(len(count_not_in_field_sample)):
                if master_cat_par['cluster'][counter] == (ii+1):                    # by cluster
                    count_not_in_field_sample[ii]+=1
    else:
        for ii in range(len(count_not_in_field_sample)):
            if master_cat_par['cluster'][counter] == (ii+1):                    # by cluster
                count_not_in_parent_sample[ii]+=1
#
#
#
#
## SECTION (4.1): SUMMARY table
##  summarize data MEMBER population as segregated above, and display in a table
#
if summary_flag_4 == 1 or adams_flag == 1:
    ## Summarize initial data stats in table
    mem_par_names = Column(['Total catalogue (PAR)','Field mem','Mem - SF','Mem - Q','Mem - Outliers','NIFS','NIPS','SUM'],name='Property')
    col_names = cluster_names_par
    mem_par0 = Column([len(master_cat_par),count_field_sample,np.sum(count_field_sample_type[0]),np.sum(count_field_sample_type[1]),np.sum(outliers_par),np.sum(count_not_in_field_sample),np.sum(count_not_in_parent_sample),np.sum([np.sum(outliers_par),np.sum(count_field_sample_type),np.sum(count_not_in_field_sample),np.sum(count_not_in_parent_sample)])],name='Total')  # total column
    mem_par_stats = Table([mem_par_names,mem_par0])
    for ii in range(len(spec_only_par)):
        mem_par_col = Column([a[ii],(count_field_sample_type[0][ii]+count_field_sample_type[1][ii]),count_field_sample_type[0][ii],count_field_sample_type[1][ii],outliers_par[ii],count_not_in_field_sample[ii],count_not_in_parent_sample[ii],np.sum([count_field_sample_type[0][ii],count_field_sample_type[1][ii],count_not_in_field_sample[ii],count_not_in_parent_sample[ii],outliers_par[ii]])],name=col_names[ii])
        mem_par_stats.add_column(mem_par_col)  # add columns to table one cluster at a time
    #
    print('\nSummary Table 4 - Catalogue by MEMBER (PAR):\n%s'%mem_par_stats,'\nNOTE: "LBLM": Lost Below Limiting Mass of parallel field.\nNOTE: "NIFS": Not In Field Sample (i.e. z>~0.6 or z<~0.25)')
    #
#
#
print('\n"master_parallel*.py" Section 4 complete.')
#
#
#
#
#
## SECTION (5) : determine LIMITING MASS for each cluster - call "data_mass_completeness*_par.py"
#
## call the "data_mass_completeness*_par.py" program, which is a carbon copy of "data_mass_completeness*.py" in each filter, set to the limiting magnitude of each parallel field
#
print('\nBeginning "data_mass_completeness*_par.py"')
#
if limiting_mass_flag == 1:
    exec(open('data_mass_completeness_F160W_par.py').read())      #opens and executes the script
elif limiting_mass_flag == 2:
    exec(open('data_mass_completeness_F814W_par.py').read())      #opens and executes the script
#
#
#
#
#
## SECTION (6) : cut all galaxies below the limiting mass of the parallel field
#
##
#
count_field_sample_mem = np.array([[0]*6]*2)      # row1=SF;  row2=Q, by cluster
SF_field_par_list = [ [], [], [], [], [], [] ]                   # THESE LISTS WILL STORE THE SMF FIELD SAMPLE MASSES
Q_field_par_list = [ [], [], [], [], [], [] ]
lost_below_limiting_mass_par = np.array([0]*6)
other_type_par = np.array([0]*6)
outside_parent_par = np.array([0]*6)
#
counting_array = np.array([0]*10)
#
## isolate all galaxies (SF & Q) in the redshift range 0.3 < z < 0.55, for the field sample of the SMF
for counter in range(len(master_cat_par)):
    counting_array[0]+=1
    if master_cat_par['sub'][counter] == 1 or master_cat_par['sub'][counter] == 2:
        if master_cat_par['member'][counter] == 1:                 # member=1: preliminary member of FIELD SAMPLE
            counting_array[1]+=1
            for cluster in range(len(count_field_sample_type[0])):
                if master_cat_par['cluster'][counter] == (cluster+1):                    # by cluster
                    counting_array[2]+=1
                    if master_cat_par['lmass'][counter] >= limiting_mass_par[cluster]:
                        if master_cat_par['type'][counter] == 1:
                            counting_array[3]+=1# SF
                            SF_field_par_list[cluster].append(master_cat_par['lmass'][counter])
                            count_field_sample_mem[0][cluster]+=1
                        elif master_cat_par['type'][counter] == 2:                      # Q
                            counting_array[4]+=1
                            Q_field_par_list[cluster].append(master_cat_par['lmass'][counter])
                            count_field_sample_mem[1][cluster]+=1
                        else:
                            other_type_par[cluster]+=1
                            counting_array[5]+=1
                    else:
                        master_cat_par['member'][counter] = 3
                        lost_below_limiting_mass_par[cluster]+=1
                        counting_array[6]+=1
                    #####
            ########
        ####
    #####
        elif master_cat_par['z_peak'][counter] > z_field_bounds[1] or master_cat_par['z_peak'][counter] < z_field_bounds[0]:
            master_cat_par['member'][counter] = 2                                   # member=2: redshift outside of field sample
            counting_array[7]+=1
            for cluster in range(len(count_not_in_field_sample)):
                if master_cat_par['cluster'][counter] == (ii+1):                    # by cluster
                    count_not_in_field_sample[cluster]+=1
                    counting_array[8]+=1
    else:
        counting_array[9]+=1
        for cluster in range(len(count_not_in_field_sample)):
            if master_cat_par['cluster'][counter] == (ii+1):                    # by cluster
                outside_parent_par[cluster]+=1
    #
#

## SECTION (6.1): SUMMARY table
##  summarize data MEMBER population as segregated above, and display in a table
#
if summary_flag_5 == 1 or adams_flag == 1:
    ## Summarize initial data stats in table
    mem_par_names = Column(['Parent sample field mem','Parent - SF','Parent - Q','Mem - SF','Mem - Q','Mem - Other type (?)','LBLM','Outliers','SUM'],name='Property')
    col_names = cluster_names_par
    mem_par0 = Column([count_field_sample,np.sum(count_field_sample_type[0]),np.sum(count_field_sample_type[1]),np.sum([count_field_sample_mem[0]]),np.sum([count_field_sample_mem[1]]),np.sum(other_type_par),np.sum(lost_below_limiting_mass_par),np.sum(outliers_par),np.sum([np.sum(count_field_sample_mem),np.sum(lost_below_limiting_mass_par),np.sum(outliers_par)])],name='Total')  # total column
    mem_par_stats = Table([mem_par_names,mem_par0])
    for ii in range(len(spec_only_par)):
        mem_par_col = Column([np.sum([count_field_sample_type[0][ii],count_field_sample_type[1][ii]]),count_field_sample_type[0][ii],count_field_sample_type[1][ii],count_field_sample_mem[0][ii],count_field_sample_mem[1][ii],other_type_par[ii],lost_below_limiting_mass_par[ii],outliers_par[ii],np.sum([count_field_sample_mem[0][ii],count_field_sample_mem[1][ii],other_type_par[ii],lost_below_limiting_mass_par[ii],outliers_par[ii]])],name=col_names[ii])
        mem_par_stats.add_column(mem_par_col)  # add columns to table one cluster at a time
    #
    print('\ncounting_array = [catalogue, members, members_in_clusters, SF_mem_clu, Q_mem_clu, Other_mem_clu, LBLM, non-members, non-members_clu, NIPS] =\m%s'%counting_array)
    #
    print('\nSummary Table 5 - Catalogue by MEMBER (PAR), post limiting mass calculation:\n%s'%mem_par_stats,'\nNOTE: "LBLM": Lost Below Limiting Mass of parallel field')
    #
#
print('\n"master_parallel*.py" Section 6 complete.')
#
#
#
## TIME_FLAG END
#
if time_flag == 1:
    if project_time_flag == 1:
        pass
    else:
        print('Program "master_parallel*.py" took: %s seconds to run.\n\n' % (time.time() - start_time))
#
#
#
#
#
print('\n\n"master_parallel*.py"  terminated successfully.\n')
#
#
#
###################     PROGRAM END     ###################
