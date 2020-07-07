#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 04:32:11 2020

@author: gsarrouh
"""
#
###This script reads in all data for the Hubble Frontier Fields images and prepares data for plotting and analysis. Key information is summarized in the tables: 
###     **"sub_stats","type_stats","SF_spec_stats","Q_spec_stats","SF_phot_stats","Q_phot_stats"**
#
### Data is organized in a single catalogue ("master_cat"), and objects are identified through applying a series of "FILTERS" to designate key populations using a numerical designation for ease of writing code.
#
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
#
#
#
##  Notes from previous versions w/ only 2 rough data sets:
### v6: calculates del_z, z_clusterphot/spec on a line-by-line basis to ensure 
### zplots file indeed plots all desired targets with spectroscopy. 
### it also fixes the problem of having loop counters incrementred before 
### carrying out whatever operation on the desired row w/in master_cat. this 
### probelm affected spec sample, all z-calculations, and and phot sample as 
### well (pervasive).
#
### v7 also replaces all individual counting variables for keeping track of 
### statistics by cluster with single arrays for more efficient computation.
# 
### v7a investigates the difference to phot sample (both for cluster members 
### AND field) of making the photometric sample cut as simply (z_phot - z_cl), 
### instead of (z_phot - z_cl) / (1 + z_cl); change made directly to 
### 'z_clusterphot for 'sub' = 2 (i.e. phot sample) in Section 2
#
###          ****  THIS IS THE ONLY DIFFERENCE B/W v7 * v7A  ****
###     ****may easily be undone by adding back denominator (1 + z_cl)****
#
##  Notes for new version w/ 6 finalized data sets
#
## v1: have undone experiment w/ z_clusterphot as noted above in v7a
## v2: removes the double-segregation of sub-samples, which lead to errors of 
##     mis-classifying targets w/ nodata/spec_only 
## v3: applies cut to photometry of targets in .cat catalogues using column 
##     "use_phot"; if use_phot ==0, don't use phot estimates, if ==1, use.
#
### Section summary:
#
### PROGRAM START
#
### (1) import data into single table, creating KEY TABLE: "master_cat" 
### (1.1)  add filter ("sieves") columns, apply SUB-TYPE FILTER in each 
###        cluster ["nodata", "phot_only","spec_only","both"], 
### (1.2)  add DIAG_FLAG_1: summarize in table "sub_stats"
### (1.3)  convert flux to mag.,
### (2) calculate various del_z's, 
### (2.1)  identify outliers
### (2.2)  compute & DISPLAY OUTLIER FRACTION, SCATTER (i.e. std dev), and MEAN of |del_z|.
### (3) distinguish b/w SF/Q: apply TYPE FILTER
### (3.1)  add DIAG_FLAG_2: summarize in table "type_stats"
### (4) make membership cuts to spec samples (i.e. apply MEMBER FILTER), apply diagnostic to test different redshift cutoff (OUTPUT FILE: /Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/section_4_false_pos_neg_redshift_cuts_*.txt) 
### (4.1) add DIAG_FLAG_3: summarize in table "SF_spec_stats" & "Q_spec_stats"
### (5) make membership cuts to phot samples
### (5.1) add DIAG_FLAG_4: summarize in table "SF_phot_stats" & "Q_phot_stats"
#
### PROGRAM END
#
#
#
###################     PROGRAM START
#
## TIME_FLAG: START
## superior time_flag which supercedes all others and times the entire program
time_flag = 0     # track & print time to execute current section
#
if time_flag == 1:
    start_time = time.time()
#  
# Read in ALL data from WORKING DIRECTORY: NSERC17/HFFtoAdam/working_data
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import astropy
from astropy.table import Table
from astropy.table import Column
import time
#
##SECTION 1: import all data from HFF team, convert flux to luminosity & gather full 
#
## TIME_FLAG_1 START
time_flag_1 = 1     # track & print time to execute current section
#
if time_flag_1 == 1 and time_flag == 0:
    start_time = time.time()
#  
## import catalogues into single table "master_cat"; separate objects for which there is no redshift data (both photo & spec) as "nodata"
#
##create table from redshift ".zout" file
z_macs0416 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/macs0416clu_v3.5.zout',format='ascii')
z_macs1149 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/macs1149clu_v3.5.zout',format='ascii')
z_macs0717 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/macs0717clu_v3.5.zout',format='ascii')
z_abell370 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/abell370clu_v3.5.zout',format='ascii')
z_abell1063 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/abell1063clu_v3.5.zout',format='ascii')   
z_abell2744 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/abell2744clu_v3.5.zout',format='ascii') 
#   
#create table from redshift ".fout" file
f_macs0416 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/macs0416clu_v3.5.fout',format='ascii')
f_macs1149 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/macs1149clu_v3.5.fout',format='ascii')
f_macs0717 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/macs0717clu_v3.5.fout',format='ascii')
f_abell370 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/abell370clu_v3.5.fout',format='ascii')
f_abell1063 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/abell1063clu_v3.5.fout',format='ascii')
f_abell2744 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/abell2744clu_v3.5.fout',format='ascii')
#
##read in the whole bloody catalogue to sift for stars
cat_macs0416 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/hffds_macs0416clu_v3.5.cat',format='ascii')
cat_macs1149 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/hffds_macs1149clu_v3.5.cat',format='ascii')
cat_macs0717 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/hffds_macs0717clu_v3.5.cat',format='ascii')
cat_abell370 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/hffds_abell370clu_v3.5.cat',format='ascii')
cat_abell1063 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/hffds_abell1063clu_v3.5.cat',format='ascii')
cat_abell2744 = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/hffds_abell2744clu_v3.5.cat',format='ascii')
#
##creat table for colours
#macs0416
F0416_u = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/macs0416clu_v3.5.153.rf',format='ascii')
F0416_v = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/macs0416clu_v3.5.155.rf',format='ascii')
F0416_j = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/macs0416clu_v3.5.161.rf',format='ascii')
#macs1149
F1149_u = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/macs1149clu_v3.5.153.rf',format='ascii')
F1149_v = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/macs1149clu_v3.5.155.rf',format='ascii')
F1149_j = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/macs1149clu_v3.5.161.rf',format='ascii')
#macs0717
F0717_u = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/macs0717clu_v3.5.153.rf',format='ascii')
F0717_v = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/macs0717clu_v3.5.155.rf',format='ascii')   
F0717_j = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/macs0717clu_v3.5.161.rf',format='ascii')
#abell370
F370_u = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/abell370clu_v3.5.153.rf',format='ascii')
F370_v = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/abell370clu_v3.5.155.rf',format='ascii')
F370_j = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/abell370clu_v3.5.161.rf',format='ascii')
#abell1063
F1063_u = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/abell1063clu_v3.5.153.rf',format='ascii')
F1063_v = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/abell1063clu_v3.5.155.rf',format='ascii')
F1063_j = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/abell1063clu_v3.5.161.rf',format='ascii')
#abell2744
F2744_u = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/abell2744clu_v3.5.153.rf',format='ascii')
F2744_v = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/abell2744clu_v3.5.155.rf',format='ascii')
F2744_j = Table.read('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/abell2744clu_v3.5.161.rf',format='ascii')
##aggregate into a single table
macs0416 = Table([z_macs0416['id'],z_macs0416['z_peak'],z_macs0416['z_spec'],F0416_u['L153'],F0416_v['L155'],F0416_j['L161'],F0416_u['DM'],f_macs0416['lmass'],f_macs0416['lsfr'],f_macs0416['lssfr'],cat_macs0416['flux_radius'],cat_macs0416['star_flag'],cat_macs0416['use_KS'],cat_macs0416['use_phot'],cat_macs0416['f_F160W'],cat_macs0416['e_F160W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','BCG','use_phot','f_F160W','e_F160W'))
macs1149 = Table([z_macs1149['id'],z_macs1149['z_peak'],z_macs1149['z_spec'],F1149_u['L153'],F1149_v['L155'],F1149_j['L161'],F1149_u['DM'],f_macs1149['lmass'],f_macs1149['lsfr'],f_macs1149['lssfr'],cat_macs1149['flux_radius'],cat_macs1149['star_flag'],cat_macs1149['use_KS'],cat_macs1149['use_phot'],cat_macs1149['f_F160W'],cat_macs1149['e_F160W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','BCG','use_phot','f_F160W','e_F160W'))
macs0717 = Table([z_macs0717['id'],z_macs0717['z_peak'],z_macs0717['z_spec'],F0717_u['L153'],F0717_v['L155'],F0717_j['L161'],F0717_u['DM'],f_macs0717['lmass'],f_macs0717['lsfr'],f_macs0717['lssfr'],cat_macs0717['flux_radius'],cat_macs0717['star_flag'],cat_macs0717['use_KS'],cat_macs0717['use_phot'],cat_macs0717['f_F160W'],cat_macs0717['e_F160W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','BCG','use_phot','f_F160W','e_F160W'))
abell370 = Table([z_abell370['id'],z_abell370['z_peak'],z_abell370['z_spec'],F370_u['L153'],F370_v['L155'],F370_j['L161'],F370_u['DM'],f_abell370['lmass'],f_abell370['lsfr'],f_abell370['lssfr'],cat_abell370['flux_radius'],cat_abell370['star_flag'],cat_abell370['use_KS'],cat_abell370['use_phot'],cat_abell370['f_F160W'],cat_abell370['e_F160W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','BCG','use_phot','f_F160W','e_F160W'))
abell1063 = Table([z_abell1063['id'],z_abell1063['z_peak'],z_abell1063['z_spec'],F1063_u['L153'],F1063_v['L155'],F1063_j['L161'],F1063_u['DM'],f_abell1063['lmass'],f_abell1063['lsfr'],f_abell1063['lssfr'],cat_abell1063['flux_radius'],cat_abell1063['star_flag'],cat_abell1063['use_KS'],cat_abell1063['use_phot'],cat_abell1063['f_F160W'],cat_abell1063['e_F160W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','BCG','use_phot','f_F160W','e_F160W'))
abell2744 = Table([z_abell2744['id'],z_abell2744['z_peak'],z_abell2744['z_spec'],F2744_u['L153'],F2744_v['L155'],F2744_j['L161'],F2744_u['DM'],f_abell2744['lmass'],f_abell2744['lsfr'],f_abell2744['lssfr'],cat_abell2744['flux_radius'],cat_abell2744['star_flag'],cat_abell2744['use_KS'],cat_abell2744['use_phot'],cat_abell2744['f_F160W'],cat_abell2744['e_F160W']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','BCG','use_phot','f_F160W','e_F160W'))
#
## create columns and append to master_cat to identify sub-sample, type, 
## and fraction for catalogue: 
##   type: 0 = stars,  1 = SF,  2 = Q
##   member:  0 = not in spec sub-sample,  1 = secure in cluster,  2 = secure in field
##            3 = false positive,  4 = false negative
## Note: "member" only applies to the spec sub-sample
D1 = Column([1]*len(macs0416),name='cluster')
D2 = Column([2]*len(macs1149),name='cluster')        #cluster designation columns
D3 = Column([3]*len(macs0717),name='cluster')
D4 = Column([4]*len(abell370),name='cluster')
D5 = Column([5]*len(abell1063),name='cluster')
D6 = Column([6]*len(abell2744),name='cluster')
#
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
error = 0
both = np.array([0]*6)
no_data = np.array([0]*6)
stars_sub = np.array([0]*6)
#
for counter in range(len(master_cat)):
    if master_cat['z_spec'][counter] > 0 and master_cat['z_peak'][counter] > 0:  #entries w/ both spectroscopy and photometry  (SPECTROSCOPIC sub-sample)
        if master_cat['use_phot'][counter] ==0:     # use_phot = 0 means bad photometry
            for jj in range(len(spec_only)):
                if master_cat['cluster'][counter] == (jj+1):      # identify # of spec only objects by cluster
                    if master_cat['star_flag'][counter] == 1:
                        master_cat['sub'][counter] = 4          # APPLY FILTER: sub=4 for STARS
                        stars_sub[jj]+=1
                    else:
                        spec_only[(jj)]+=1
            master_cat['sub'][counter] = 3          # APPLY FILTER: sub=3 for objects w/ SPEC ONLY
        elif master_cat['use_phot'][counter] ==1:
            for jj in range(len(both)):
                if master_cat['cluster'][counter] == (jj+1):      # identify # of (spec & phot) objects by cluster
                    if master_cat['star_flag'][counter] == 1:
                        master_cat['sub'][counter] = 4          # APPLY FILTER: sub=4 for STARS
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
                    if master_cat['star_flag'][counter] == 1:
                        master_cat['sub'][counter] = 4          # APPLY FILTER: sub=4 for STARS
                        stars_sub[jj]+=1
                    else:
                        no_data[(jj)]+=1
                        master_cat['sub'][counter] = 0          # APPLY FILTER: sub=0 for objects w/ NEITHER SPEC NOR PHOT
        elif master_cat['use_phot'][counter] ==1:
            for jj in range(len(both)):
                if master_cat['cluster'][counter] == (jj+1):      # identify # of phot only objects by cluster
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
#
#
## SECTION (1.2): SUMMARY table
### MAY NEED TO EDIT ### diag_flag_1
##  summarize data population as segregated above, and display in a table
diag_flag_1 = 1             # 0=off (don't display diagnostic); 1=on (display diagnostic table)
#
if diag_flag_1 == 1:
    ## Summarize initial data stats in table
    sub_names = Column(['total','spec & phot','only phot','spec only','no data','stars'],name='Property')
    col_names = ['macs0416','macs1149','macs0717','abell370','abell1063','abell2744']
    sub0 = Column([np.sum([np.sum(both),np.sum(phot_only),np.sum(spec_only),np.sum(no_data),np.sum(stars_sub)]),np.sum(both),np.sum(phot_only),np.sum(spec_only),np.sum(no_data),np.sum(stars_sub)],name='Total')  # total column
    sub_stats = Table([data_names,data0])
    for ii in range(len(spec_only)):
        sub_col = Column([np.sum([both[ii],phot_only[ii],spec_only[ii],no_data[ii],stars_sub[ii]]),both[ii],phot_only[ii],spec_only[ii],no_data[ii],stars_sub[ii]],name=col_names[ii])  # add columns to table one cluster at a time
        sub_stats.add_column(sub_col)
    #
    print('Catalogue by SUB-type:')
    print(data_stats)
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
##convert flux to magnitude (erroneously labelled as luminosity, e.g. L_u for magnitude in UV), get color indices U-V, V-J, add to table
#
for counter in range(len(master_cat)):
    if master_cat[counter]['sub'] ==0 or master_cat[counter]['sub'] ==3 or master_cat[counter]['sub'] ==4:      # skip objects w/ "no data" (sub=0); "spec only" (sub=3); "stars" (sub=4)
        pass
    else:
        master_cat['L_u'][counter] = -2.5*np.log10(master_cat['u'][counter]) + 25
        master_cat['L_v'][counter] = -2.5*np.log10(master_cat['v'][counter]) + 25
        master_cat['L_j'][counter] = -2.5*np.log10(master_cat['j'][counter]) + 25
        master_cat['uv'][counter] = master_cat['L_u'][counter] - master_cat['L_v'][counter]
        master_cat['vj'][counter] = master_cat['L_v'][counter] - master_cat['L_j'][counter]
#
## TIME_FLAG_1 END
#
if time_flag_1 == 1 and time_flag == 0:
    print("Section 1 took: %s seconds.\n\n" % (time.time() - start_time))
#
#
#
#
#
## SECTION (2): calulate DEL_Z's & separate OUTLIERS
#
## TIME_FLAG_2 START
time_flag_2 = 1     # track & print time to execute current section
#
if time_flag_2 == 1 and time_flag == 0:
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
#
#calucalte del_z, z_clusterspec, z_clusterphot for outlier cut (defined above); these will be used to make cuts (member, field, false pos/neg) to spec sample, from which we will use relative fractions by mass to correct the photometric sample for completeness.
for counter in range(len(master_cat)):
    if master_cat['sub'][counter] == 1:   # sub=1 identifies spec&phot subsample
        master_cat['del_z'][counter] = ((master_cat['z_peak'][counter] - master_cat['z_spec'][counter]) / (1 + master_cat['z_spec'][counter]))
        for ii in range(len(z_cluster)):
            if master_cat['cluster'][counter] == (ii+1):
                master_cat['z_clusterspec'][counter] = ((master_cat['z_spec'][counter] - z_cluster[ii]) / (1 + master_cat['z_spec'][counter]))
                master_cat['z_clusterphot'][counter] = ((master_cat['z_peak'][counter] - z_cluster[ii]) / (1 + master_cat['z_peak'][counter]))
    elif master_cat['sub'][counter] == 2:   # sub=2 identifies phot-only subsample
        for ii in range(len(z_cluster)):
            if master_cat['cluster'][counter] == (ii+1):
                master_cat['z_clusterphot'][counter] = ((master_cat['z_peak'][counter] - z_cluster[ii]) / (1 + master_cat['z_peak'][counter]))
#
#
## SECTION(2.1): separate OUTLIERS from both phot & spec sub-sample, defined as |del_z| < 0.15. apply FILTER TYPE = 3 for outliers;   
#
outliers = np.array([0]*6)      # initialize array to track outliers by cluster, for computing outlier fraction later
sum_delz = []                   # for computing mean |del_z|

#
for counter in range(len(master_cat)):
    if master_cat['sub'][counter] == 1:                        # sub=1 for objects with both spec & phot; total # of such objects tracked by cluster above in the array called "both"
        if np.abs(master_cat['del_z'][counter]) > 0.15:        # |del_z| > 0.15 for outliers, threshold chosen for historical reasons to facilitate comparison with other studies
            master_cat['type'][counter] = 3                  # type = 3 identifies outliers
            for ii in range(len(outliers)):
                if master_cat['cluster'][counter] == (ii+1):   # keep track of outliers by cluster
                    outliers[ii]+=1
        sum_delz.append(np.abs(master_cat['del_z'][counter]))
#
#
## SECTION(2.2): compute & DISPLAY OUTLIER FRACTION, SCATTER (i.e. std dev), and MEAN of |del_z|.
#
delz_mean = np.sum(sum_delz)/len(sum_delz)
delz_scatter = np.std(sum_delz)
print('OUTLIERS: %s' % np.sum(outliers))
print('Outlier fraction: %s' % (np.sum(outliers)/np.sum(both)))
print('|del_z| mean: %s'%delz_mean)
print('|del_z| scatter: %s\n'%delz_scatter)
#
## TIME_FLAG_2 END
#
if time_flag_2 == 1 and time_flag == 0:
    print("Section 2 took: %s seconds.\n\n" % (time.time() - start_time))
#  
#
#
#
#
## SECTION (3): add FILTER TYPE: separate SF/Q for both subsamples based on van der Burg (2013) 
## colour criteria;    filter name: 'type'  IDs below
##   0 = stars;  1 = SF (star-forming);    2 = Q (quiscient);    3 = outliers
#
#
## TIME_FLAG_3 START
time_flag_3 = 1     # track & print time to execute current section
#
if time_flag_3 == 1 and time_flag == 0:
    start_time = time.time()
#  
SF_type = np.array([0]*6)                             # initialize arrays
Q_type = np.array([0]*6)
stars_type = np.array([0]*6)                          # count # of stars by cluster
lost_type = np.array([0]*6)                           # objects lost due no data (sub=0), spec only (sub=3), 
for counter in range(len(master_cat)):
    if master_cat['star_flag'][counter] == 1:          # identify STARS, filter type=0
        master_cat['type'][counter] = 0              
        for ii in range(len(stars_type)):
            if master_cat['cluster'][counter] == (ii+1):   # keep stars of outliers by cluster
                stars_type[ii]+=1
    elif master_cat['sub'][counter]==1 or master_cat['sub'][counter]==2:    # sub=1 for (spec & phot) subsample, sub=2 for phot subsample; i.e. look at all objects with photometry
        if master_cat['uv'][counter] > 1.3 and master_cat['vj'][counter] < 1.6 and master_cat['uv'][counter] > ((0.88*master_cat[counter]['vj']) + 0.6): 
            master_cat['type'][counter] = 2             # identify passive (QUIESCENT) galaxies, type=2
            for ii in range(len(Q_type)):
                if master_cat['cluster'][counter] == (ii+1):   # keep stars of outliers by cluster
                    Q_type[ii]+=1
        else:
            master_cat['type'][counter] = 1             # identify STAR-FORMING galaxies, type=1
            for ii in range(len(SF_type)):
                if master_cat['cluster'][counter] == (ii+1):   # keep stars of outliers by cluster
                    SF_type[ii]+=1
    else:
        for ii in range(len(lost_type)):
            if master_cat['cluster'][counter] == (ii+1):   # keep stars of outliers by cluster
                lost_type[ii]+=1      # objects lost due to spec only (sub=3),no data (sub=0), and possibly outliers (type=3) 
#        
#
## SECTION (3.1): SUMMARY table
### MAY NEED TO EDIT ### diag_flag_2
##  summarize data population as segregated above, and display in a table
diag_flag_2 = 1
#
if diag_flag_2 == 1:
    ## Summarize initial data stats in table
    type_names = Column(['Total','SF','Q','Other'],name='Property')
    col_names = ['macs0416','macs1149','macs0717','abell370','abell1063','abell2744']
    type0 = Column([np.sum([np.sum(SF_type),np.sum(Q_type),np.sum(lost_type),np.sum(stars_type)]),np.sum(SF_type),np.sum(Q_type),np.sum(lost_type+stars_type)],name='Total')  # total column
    type_stats = Table([type_names,type0])
    for ii in range(len(spec_only)):
        type_col = Column([np.sum([SF_type[ii],Q_type[ii],lost_type[ii],stars_type[ii]]),SF_type[ii],Q_type[ii],(lost_type[ii]+stars_type[ii])],name=col_names[ii])
        type_stats.add_column(type_col)  # add columns to table one cluster at a time
    #
    print('Catalogue by TYPE:')
    print(type_stats)
    print('NOTE: "Other" is comprised of objects without photometry (i.e. stars, and objects with either bad photometry or both bad photometry and bad spectroscopy).\n')
#
## TIME_FLAG_3 END
#
if time_flag_3 == 1 and time_flag == 0:
    print("Section 3 took: %s seconds.\n\n" % (time.time() - start_time))
#
#
#
#
#
## SECTION 4: apply MEMBER FILTER based on cuts discussed with AM in 2017 (see comments below), based on cuts made in VDB2013. MEMBER isolates: 0=cluster member (secure); 1=field (secure); 2=false pos; 3=false neg;
#
## CRITERIA
#
## SF: cluster = 0: abs(z_clusterspec) < 0.01 & abs(z_cluster phot) < 0.03; 
##      field = 1: abs(z_clusterspec) > 0.02 & abs(z_cluster phot) > 0.1; 
## false pos = 2: abs(z_clusterspec) > 0.01 & abs(z_cluster phot) < 0.03;
## false neg = 3: abs(z_clusterspec) < 0.01 & abs(z_cluster phot) > 0.03;
## Q: same cutoff for z_clusterspec, cutoff for z_clusterphot > 0.07
#
## Note: this is only done for spec sample; results used to correct photo sample for completeness at end of analysis in file "master_smf_8.py"
#
### SF sub-sample: SF_spec into secure member, field, false pos/neg
#
#
## TIME_FLAG_4 START
time_flag_4 = 1     # track & print time to execute current section
#
if time_flag_4 == 1 and time_flag == 0:
    start_time = time.time()
#
## the following code is part of a diagnostic to test the # of false pos/neg produced by altering the photometric/spectroscopic redshift cutoff. upon completion it will be commented out permenantly.
#
# define cut-offs for SF & Q
z_cutoff = np.linspace(0.01,0.09,9)       # create array from [0.01,0.09] in steps of 0.01; replace in loop below with SF_cutoff & Q_cutoff once cutoffs are determined 
SF_cutoff = [0.01,0.03]      # [spec,phot]
Q_cutoff = [0.01,0.08]
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
##
## open a file to print to
#f = open('/Users/gsarrouh/Documents/Programs/Python/nserc17/HFF_ToAdamFinal/working_data/section_4_false_pos_neg_redshift_cuts_001.txt','w+')
#
#for cutoff in range(len(z_cutoff)):
##
## INDENT HERE for diagnostic 
mem = np.array([[0]*6]*2)       # initialize arrays to track cluster members, field, false pos/neg by cluster
field = np.array([[0]*6]*2)     # row_1=SF; row_2=Q; for all arrays
pos = np.array([[0]*6]*2)
neg = np.array([[0]*6]*2)
other_member = 0                 # track objects outside of (phot + spec) subsample
lost_due_to_buffer = np.array([[0]*6]*2)
#
## The following loop isolates the (spec + phot) sample, i.e. 'sub'=1, and makes the cuts defined above, assigning different classifications to the MEMBER FILTER 
for counter in range(len(master_cat)):
    if master_cat['sub'][counter] == 1:                   # sub=1 identifies subsample with both spec & phot
        if master_cat['type'][counter]==1:                # type=1 identifies SF sample
            if abs(master_cat['z_clusterspec'][counter]) > 0.02 and abs(master_cat['z_clusterphot'][counter]) > 0.1: 
                master_cat['member'][counter] = 1         # member=1 for FIELD
                for ii in range(len(field[0])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                        field[0][ii]+=1
            elif abs(master_cat['z_clusterspec'][counter]) > SF_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) < SF_cutoff[1]: #z_cutoff[cutoff]: #
                master_cat['member'][counter] = 2         # member=2 for FALSE POSITIVE
                for ii in range(len(pos[0])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of false pos by cluster
                        pos[0][ii]+=1
            elif abs(master_cat['z_clusterspec'][counter]) < SF_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) > SF_cutoff[1]: #z_cutoff[cutoff]: #
                master_cat['member'][counter] = 3         # member=3 for FALSE NEGATIVE
                for ii in range(len(neg[0])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of false neg by cluster
                        neg[0][ii]+=1
            elif abs(master_cat['z_clusterspec'][counter]) < SF_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) < SF_cutoff[1]: #z_cutoff[cutoff]: #
                master_cat['member'][counter] = 0         # member=0 for cluster MEMBERS
                for ii in range(len(mem[0])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of cluster members by cluster
                        mem[0][ii]+=1
            else: 
                for ii in range(len(lost_due_to_buffer[0])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of cluster members by cluster
                        lost_due_to_buffer[0][ii]+=1
        elif master_cat['type'][counter]==2:                # type=2 identifies Q sample
            if abs(master_cat['z_clusterspec'][counter]) > 0.02 and abs(master_cat['z_clusterphot'][counter]) > 0.1: 
                master_cat['member'][counter] = 1         # member=1 for FIELD
                for ii in range(len(field[1])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                        field[1][ii]+=1
            elif abs(master_cat['z_clusterspec'][counter]) > Q_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) < Q_cutoff[1]: #z_cutoff[cutoff]: #
                master_cat['member'][counter] = 2         # member=2 for FALSE POSITIVE
                for ii in range(len(pos[1])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of false pos by cluster
                        pos[1][ii]+=1
            elif abs(master_cat['z_clusterspec'][counter]) < Q_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) > Q_cutoff[1]: #z_cutoff[cutoff]: #
                master_cat['member'][counter] = 3         # member=3 for FALSE NEGATIVE
                for ii in range(len(neg[1])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of false neg by cluster
                        neg[1][ii]+=1
            elif abs(master_cat['z_clusterspec'][counter]) < Q_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) < Q_cutoff[1]: #z_cutoff[cutoff]: #
                master_cat['member'][counter] = 0         # member=0 for cluster MEMBERS
                for ii in range(len(mem[1])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of cluster members by cluster
                        mem[1][ii]+=1
            else: 
                for ii in range(len(lost_due_to_buffer[1])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of cluster members by cluster
                        lost_due_to_buffer[1][ii]+=1
    else: other_member+=1
## INDENT HERE
#    #
#    # compute membership acceptance fraction
#    SF_specphot = np.sum([np.sum(mem[0]),np.sum(field[0]),np.sum(pos[0]),np.sum(neg[0])]) # total # of galaxies with both spec & phot - SF
#    Q_specphot = np.sum([np.sum(mem[1]),np.sum(field[1]),np.sum(pos[1]),np.sum(neg[1])])  # total # of galaxies with both spec & phot - Q
#    mem_fraction = np.array([(np.sum(mem[0])/SF_specphot),(np.sum(mem[1])/Q_specphot)]) # [SF acceptance fraction, Q acceptance fraction]
#    print('\nOverall membership fraction: \nSF: %s'%mem_fraction[0],' & Q: %s'%mem_fraction[1],'   for cutoff:  ',str(z_cutoff[cutoff]))
#    # print # of memebrs, false pos/neg, & acceptance fraction for each of SF & Q to a file
#    # FORMAT: row_1 = SF;   row_2 = Q;
#    #   col_1 = z_cut;  col_2 = TYPE;  col_3=members; col_4 = false pos; col_5 = false neg; col_6 = acceptance fraction; col_7 = False pos/neg ratio
#    space = ' '
#    zcut = str(z_cutoff[cutoff])
#    SF_1 = str(np.sum(mem[0]))
#    SF_2 = str(np.sum(pos[0]))
#    SF_3 = str(np.sum(neg[0]))
#    SF_4 = str(mem_fraction[0])
#    SF_5 = str(max((np.sum(pos[0])/np.sum(neg[0])),(np.sum(neg[0])/np.sum(pos[0]))))
#    a = zcut+space+'SF'+space+SF_1+space+SF_2+space+SF_3+space+SF_4+space+SF_5
#    Q_1 = str(np.sum(mem[1]))
#    Q_2 = str(np.sum(pos[1]))
#    Q_3 = str(np.sum(neg[1]))
#    Q_4 = str(mem_fraction[1])
#    Q_5 = str(max((np.sum(pos[1])/np.sum(neg[1])),(np.sum(neg[1])/np.sum(pos[1]))))
#    b = zcut+space+'Q'+space+Q_1+space+Q_2+space+Q_3+space+Q_4+space+Q_5
#    if cutoff == 0:
#        header1 = '### The del_z spec cutoff has been set to > %s for all runs in this document.'%SF_cutoff[0]
#        header2 = '### Columns:  del_z cutoff  Type   Member   False pos.   False neg.   % acceptance   #   False pos/neg ratio\n#\n#\n'
#        writer = '%s\n'%header1+'%s'%header2+'%s\n'%a+'%s\n'%b
#        f.write(writer)
#    else: 
#        writer = '%s\n'%a+'%s\n'%b
#        f.write(writer)
##        
#f.close()
#
#
# END of DIAGNOSTIC loop
#
#    
## SECTION (4.1): SUMMARY table
### MAY NEED TO EDIT ### diag_flag_3
##  summarize data population as segregated above, and display in a table
diag_flag_3 = 1
#
if diag_flag_3 == 1:
    ## Summarize initial data stats in table
    member_names = Column(['Total','Secure member','Secure field','Fasle pos','False neg'],name='Property')
    col_names = ['macs0416','macs1149','macs0717','abell370','abell1063','abell2744']
    # SF table
    member0 = Column([np.sum(both),np.sum(mem[0]),np.sum(field[0]),np.sum(pos[0]),np.sum(neg[0])],name='Total')  # total column
    SF_spec_stats = Table([member_names,member0])
    for ii in range(len(mem[0])):
        col = Column([both[ii],mem[0][ii],field[0][ii],pos[0][ii],neg[0][ii]],name=col_names[ii])
        SF_spec_stats.add_column(col)  # add columns to table one cluster at a time
    #
    # Q table
    member0 = Column([np.sum(both),np.sum(mem[1]),np.sum(field[1]),np.sum(pos[1]),np.sum(neg[1])],name='Total')  # total column
    Q_spec_stats = Table([member_names,member0])
    for ii in range(len(mem[1])):
        col = Column([both[ii],mem[1][ii],field[1][ii],pos[1][ii],neg[1][ii]],name=col_names[ii])
        Q_spec_stats.add_column(col)  # add columns to table one cluster at a time
    #
    print('(SPEC+PHOT) Subsample\nCatalogue by MEMBER - Star-forming:')
    print(SF_spec_stats)
    print('NOTE: Total reported under each cluster is sum of SF+Q.\n')
    print('Catalogue by MEMBER - Quiescent:')
    print(Q_spec_stats)
    print('NOTE: Total reported under each cluster is sum of SF+Q.')
    print('NOTE: Differences b/w Total row and sum of other rows might arise due to the "buffer" zone built in between classifying objects as secure member vs field.\n')
#
SF_specphot = np.sum([np.sum(mem[0]),np.sum(field[0]),np.sum(pos[0]),np.sum(neg[0])])
Q_specphot = np.sum([np.sum(mem[1]),np.sum(field[1]),np.sum(pos[1]),np.sum(neg[1])])
print('\nOverall membership fraction: \nSF: %s'%(np.sum(mem[0])/SF_specphot),' & Q: %s'%(np.sum(mem[1])/Q_specphot),'   for cutoff:\nSF: ',str(SF_cutoff),'    Q: ',str(Q_cutoff))
#
print('\nTotal catalogue length: %s'%len(master_cat))
print('SPEC+PHOT sub-sample: %s' % np.sum(both))
print('SF: %s' % np.sum([mem[0],field[0],pos[0],neg[0]]))
print('Q: %s' % np.sum([mem[1],field[1],pos[1],neg[1]]))
print('Lost due to buffer b/w member & field\nSF: %s'%np.sum(lost_due_to_buffer[0]),';    Q: %s'%np.sum(lost_due_to_buffer[1]))
print('Other (not in (spec + phot) subsample): %s'%other_member)
print('NOTE: Differences b/w Total row and sum of other rows might arise due to the "buffer" zone built in between classifying objects as secure member vs field.\n')                        
#
## TIME_FLAG_4 END
#
if time_flag_4 == 1 and time_flag == 0:
    print("Section 4 took: %s seconds.\n\n" % (time.time() - start_time))
#
#
#
#
#
## SETION (5): make cuts to photometric sample based on defintion for del_z:  del_z = (z_phot - z_cl) < 0.05 (to match cut made above in SF_cutoff[1]/Q_cutoff[1]; apply MEMBER FILTER = 3 for PHOTOMETRIC SUBSAMPLE MEMBERS); same photometric cut made to specroscopic sub-sample. this is a preliminary measure for determining the photometric sample, final corrections will be made by mass bins to match false pos/neg fractions in spec. sample per van der Burg (2013)
#
## apply cut at |z_clusterphot| < SF_cutoff[1]/Q_cutoff[1] to separate cluster members from field for targets with photometry only. store in MEMBER FILTER: 
## 0 = cluster member; 1 = field; 2 = false pos; 3 = false neg; 4 = field outlier
#
## recall: z_clusterphot defined as (z_peak - z_cluster / 1 + z_peak)
#
#
## TIME_FLAG_5 START
time_flag_5 = 1     # track & print time to execute current section
#
if time_flag_5 == 1 and time_flag == 0:
    start_time = time.time()
#  
### Quiescent sub-sample: segregate into secure member, field
# initialize arrays, format is row_1=SF, row_2=Q
mem_phot = np.array([[0]*6]*2)
field_phot = np.array([[0]*6]*2)
other_phot = 0      # objects not in sub=2 phot only subsample
n_phot_only = 0     # number of objects in sub=2 subsample
n_SF = 0
n_Q = 0
lost_due_to_buffer_phot = np.array([[0]*6]*2)    # objects lost due to buffer b/w definition of cluster and field
#
for counter in range(len(master_cat)):
    if master_cat[counter]['sub'] ==2:      # sub=2 identifies phot only subsample
        n_phot_only+=1
        if master_cat[counter]['type'] == 0 or master_cat[counter]['type'] ==3: #skip stars and outliers
            pass       
        elif master_cat[counter]['type'] ==1:       # type=1 identifies SF
            n_SF+=1
            if abs(master_cat[counter]['z_clusterphot']) > SF_cutoff[1]:     # identify field galaxies
                if master_cat[counter]['z_peak'] >0.55 or master_cat[counter]['z_peak'] <0.3:
                    master_cat[counter]['member'] = 4       # memfield outlier
                else:
                    master_cat[counter]['member'] = 1               #phot SF field sample
                    for ii in range(len(field_phot[0])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field galaxies by cluster
                            field_phot[0][ii]+=1
            elif abs(master_cat[counter]['z_clusterphot']) < SF_cutoff[1]:
                master_cat[counter]['member'] = 0           # member=0 is secure cluster member
                for ii in range(len(mem_phot[0])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of field galaxies by cluster
                        mem_phot[0][ii]+=1
            else:
                for ii in range(len(lost_due_to_buffer_phot[0])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of field galaxies by cluster
                        lost_due_to_buffer_phot[0][ii]+=1
        elif master_cat[counter]['type'] ==2:       #Q
            n_Q+=1
            if abs(master_cat[counter]['z_clusterphot']) > Q_cutoff[1]:     # identify field galaxies
                if master_cat[counter]['z_peak'] >0.55 or master_cat[counter]['z_peak'] <0.3:
                    master_cat[counter]['member'] = 4       # memfield outlier
                else:
                    master_cat[counter]['member'] = 1               #phot SF field sample
                    for ii in range(len(field_phot[1])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field galaxies by cluster
                            field_phot[1][ii]+=1
            elif abs(master_cat[counter]['z_clusterphot']) < Q_cutoff[1]:
                master_cat[counter]['member'] = 0           # member=0 is secure cluster member
                for ii in range(len(mem_phot[1])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of field galaxies by cluster
                        mem_phot[1][ii]+=1
            else:
                for ii in range(len(lost_due_to_buffer_phot[1])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of field galaxies by cluster
                        lost_due_to_buffer_phot[1][ii]+=1
    else: 
        other_phot+=1
#                        
## SECTION (5.1): SUMMARY table
### MAY NEED TO EDIT ### diag_flag_4
##  summarize data population as segregated above, and display in a table
diag_flag_4 = 1
#
if diag_flag_4 == 1:
    ## Summarize initial data stats in table
    member_names = Column(['Total','Secure member','Secure field'],name='Property')
    col_names = ['macs0416','macs1149','macs0717','abell370','abell1063','abell2744']
    # SF table
    member0 = Column([np.sum([mem_phot[0],field_phot[0]]),np.sum(mem_phot[0]),np.sum(field_phot[0])],name='Total')  # total column
    SF_phot_stats = Table([member_names,member0])
    for ii in range(len(mem_phot[0])):
        col = Column([np.sum([mem_phot[0][ii],field_phot[0][ii]]),mem_phot[0][ii],field_phot[0][ii]],name=col_names[ii])
        SF_phot_stats.add_column(col)  # add columns to table one cluster at a time
    #
    # Q table
    member0 = Column([np.sum([mem_phot[1],field_phot[1]]),np.sum(mem_phot[1]),np.sum(field_phot[1])],name='Total')  # total column
    Q_phot_stats = Table([member_names,member0])
    for ii in range(len(mem_phot[1])):
        col = Column([np.sum([mem_phot[1][ii],field_phot[1][ii]]),mem_phot[1][ii],field_phot[1][ii]],name=col_names[ii])
        Q_phot_stats.add_column(col)  # add columns to table one cluster at a time
    #
    print('PHOT-ONLY Subsample\nCatalogue by MEMBER - Star-forming:')
    print(SF_phot_stats)
    print('NOTE: Total reported under each cluster is sum of SF+Q.\n')
    print('Catalogue by MEMBER - Quiescent:')
    print(Q_phot_stats)
    print('Lost due to buffer b/w member & field\nSF: %s'%np.sum(lost_due_to_buffer_phot[0]),';    Q: %s'%np.sum(lost_due_to_buffer_phot[1]))
    print('NOTE: Total reported under each cluster is sum of SF+Q.')
    print('NOTE: Differences b/w Total row and sum of other rows might arise due to the "buffer" zone built in between classifying objects as secure member vs field.\n')
#
mem_phot_fraction = np.array([0]*2,dtype='float32')     # to keep track of membership acceptance fraction
mem_phot_fraction[0] = (np.sum(mem_phot[0]) / n_SF)
mem_phot_fraction[1] = (np.sum(mem_phot[1]) / n_Q)
#
print('\nOverall membership fraction: \nSF: %s'%mem_phot_fraction[0],' & Q: %s'%mem_phot_fraction[1],'   for cutoff:\nSF: ',str(SF_cutoff),'    Q: ',str(Q_cutoff))
#
print('\nTotal catalogue length: %s'%len(master_cat))
print('PHOT ONLY sub-sample: %s' %n_phot_only)
print('SF: %s' % np.sum([mem_phot[0],field_phot[0]]))
print('Q: %s' % np.sum([mem_phot[1],field_phot[1]]))
print('Other (not in phot only subsample): %s'%other_phot)
print('NOTE: Differences b/w Total row and sum of other rows might arise due to the "buffer" zone built in between classifying objects as secure member vs field.\n')                        
#
## TIME_FLAG_5 END
#
if time_flag_5 == 1 and time_flag == 0:
    print("Section 5 took: %s seconds.\n\n" % (time.time() - start_time))
#
#
#
#
#                        
## SECTION (6): BCGs. Brightest Cluster Galaxies (BCGs) need to be taken out of the cluster sample as they are unique to overly dense environments and so lack a counterpart in the field against which to make a fair comparison. as such, we remove them from our sample before making the SMF
#
#
## TIME_FLAG_6 START
time_flag_6 = 1     # track & print time to execute current section
#
if time_flag_6 == 1 and time_flag == 0:
    start_time = time.time()
#
#
#####       UNFINISHED        #####
##### INSERT code to model out BCGs
#
#
## TIME_FLAG_6 END
time_flag_6 = 1     # track & print time to execute current section
#
if time_flag_6 == 1 and time_flag == 0:
    print("Section 6 took: %s seconds.\n\n" % (time.time() - start_time))
#
#
#
#
#
## TIME_FLAG END
#
if time_flag == 1:
    print('Program "master_data_6_final.py" took: %s seconds to run.\n\n' % (time.time() - start_time))
#
#
#
#
#                        
###### PROGRAM END ######