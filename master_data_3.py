#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 15:19:59 2017

@author: gsarrouh
"""
##      master_data_3
#
###This script reads in all data for the Hubble Frontier Fields images and 
###prepares data for plotting and analysis. Key information is summarized 
###in the tables: 
###     **"data_stats","macros_type","macros","spec_stats","phot_stats",****
###     **"sample_stats","BCG_stats"
#
###Data is organized in a single catalogue ("master_cat"), and objects are 
###identified through applying a series of filters or "sieves" to designate key populations 
###using a numerical designation for ease of writing code.
#
### NOTE: tagets w/o spec or photometry have been removed prior to creating "master_cat" table
#
## Sieve 1 - 'cluster':  cluster catalogues are designated as follows: 
## 1: macs0416
## 2: macs1149
## 3: macs0717
## 4: abell370
## 5: abell1063
## 6: abell2744
#
## Sieve 2 - 'sub': :   identifies type of data each object has
## 0: no data (no photometry or spectroscopy)
## 1: spectroscopy & photometry
## 2: photometry only 
## 3: spectroscopy only  (there's just 3 in macs0717)   
#
## Sieve 3 - 'type': :   identifies type of data each object has
## 0: star
## 1: star-forming (SF)
## 2: quiescent (Q)  
## 3: outliers (defined as |del_z\ > 0.15)
#
## Sieve 4 - 'member': :   identifies type of data each object has
## 0: secure cluster member
## 1: secure field    <-- this comprises the sample of field galaxies to be 
##                        compared with cluster, at similar redshifts
## 2: false positive
## 3: false positive
## 4: field outliers  <-- these are galaxies well outside the redshift range 
##                        of the 6 clusters
#
## NOTE: targets w/ only photometry are split into 'cluster members' and 'field',
##       no z_spec for determining false pos/neg
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
### (1) import data, convert flux to mag., create master_cat, add filter 
###     ("sieves") columns, identify sub-samples in each cluster ["nodata",
###     "phot_only","spec_only","both"], summarize in "data_stats" 
### (2) identify outliers, calculate various del_z's, summarize 
###     in "macros"
### (3) distinguish b/w SF/Q
### (4) make membership cuts to spec samples, summarize in "spec_stats"
### (5) make membership cuts to phot samples, summarize in "phot_stats"
#
#
###################     Start
# Read in ALL data from WORKING DIRECTORY: NSERC17/HFFtoAdam/working_data
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import astropy
from astropy.table import Table
from astropy.table import Column
#
##Section 1: import all data from HFF team, convert flux to luminosity & gather full 
##catalogue into single table "master_cat"; separate objects for which there is no 
##redshift data (both photo & spec) as "nodata"
#
##create table from redshift ".zout" file
z_macs0416 = Table.read('macs0416clu_v3.5.zout',format='ascii')
z_macs1149 = Table.read('macs1149clu_v3.5.zout',format='ascii')
z_macs0717 = Table.read('macs0717clu_v3.5.zout',format='ascii')
z_abell370 = Table.read('abell370clu_v3.5.zout',format='ascii')
z_abell1063 = Table.read('abell1063clu_v3.5.zout',format='ascii')   
z_abell2744 = Table.read('abell2744clu_v3.5.zout',format='ascii')   
#insert other 4 catalogues
#create table from redshift ".fout" file
f_macs0416 = Table.read('macs0416clu_v3.5.fout',format='ascii')
f_macs1149 = Table.read('macs1149clu_v3.5.fout',format='ascii')
f_macs0717 = Table.read('macs0717clu_v3.5.fout',format='ascii')
f_abell370 = Table.read('abell370clu_v3.5.fout',format='ascii')
f_abell1063 = Table.read('abell1063clu_v3.5.fout',format='ascii')
f_abell2744 = Table.read('abell2744clu_v3.5.fout',format='ascii')
#
##read in the whole bloody catalogue to sift for stars
cat_macs0416 = Table.read('hffds_macs0416clu_v3.5.cat',format='ascii')
cat_macs1149 = Table.read('hffds_macs1149clu_v3.5.cat',format='ascii')
cat_macs0717 = Table.read('hffds_macs0717clu_v3.5.cat',format='ascii')
cat_abell370 = Table.read('hffds_abell370clu_v3.5.cat',format='ascii')
cat_abell1063 = Table.read('hffds_abell1063clu_v3.5.cat',format='ascii')
cat_abell2744 = Table.read('hffds_abell2744clu_v3.5.cat',format='ascii')

##creat table for colours
#macs0416
F0416_u = Table.read('macs0416clu_v3.5.153.rf',format='ascii')
F0416_v = Table.read('macs0416clu_v3.5.155.rf',format='ascii')
F0416_j = Table.read('macs0416clu_v3.5.161.rf',format='ascii')
#macs1149
F1149_u = Table.read('macs1149clu_v3.5.153.rf',format='ascii')
F1149_v = Table.read('macs1149clu_v3.5.155.rf',format='ascii')
F1149_j = Table.read('macs1149clu_v3.5.161.rf',format='ascii')
#macs0717
F0717_u = Table.read('macs0717clu_v3.5.153.rf',format='ascii')
F0717_v = Table.read('macs0717clu_v3.5.155.rf',format='ascii')   
F0717_j = Table.read('macs0717clu_v3.5.161.rf',format='ascii')
#abell370
F370_u = Table.read('abell370clu_v3.5.153.rf',format='ascii')
F370_v = Table.read('abell370clu_v3.5.155.rf',format='ascii')
F370_j = Table.read('abell370clu_v3.5.161.rf',format='ascii')
#abell1063
F1063_u = Table.read('abell1063clu_v3.5.153.rf',format='ascii')
F1063_v = Table.read('abell1063clu_v3.5.155.rf',format='ascii')
F1063_j = Table.read('abell1063clu_v3.5.161.rf',format='ascii')
#abell2744
F2744_u = Table.read('abell2744clu_v3.5.153.rf',format='ascii')
F2744_v = Table.read('abell2744clu_v3.5.155.rf',format='ascii')
F2744_j = Table.read('abell2744clu_v3.5.161.rf',format='ascii')
#
##aggregate into a single table
#REMOVE Z_PHOT ADJUSTMENT FOR FINAL DATA
macs0416 = Table([z_macs0416['id'],z_macs0416['z_peak'],z_macs0416['z_spec'],F0416_u['L153'],F0416_v['L155'],F0416_j['L161'],F0416_u['DM'],f_macs0416['lmass'],f_macs0416['lsfr'],f_macs0416['lssfr'],cat_macs0416['flux_radius'],cat_macs0416['star_flag'],cat_macs0416['use_KS'],cat_macs0416['use_phot']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','BCG','use_phot'))
macs1149 = Table([z_macs1149['id'],z_macs1149['z_peak'],z_macs1149['z_spec'],F1149_u['L153'],F1149_v['L155'],F1149_j['L161'],F1149_u['DM'],f_macs1149['lmass'],f_macs1149['lsfr'],f_macs1149['lssfr'],cat_macs1149['flux_radius'],cat_macs1149['star_flag'],cat_macs1149['use_KS'],cat_macs1149['use_phot']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','BCG','use_phot'))
macs0717 = Table([z_macs0717['id'],z_macs0717['z_peak'],z_macs0717['z_spec'],F0717_u['L153'],F0717_v['L155'],F0717_j['L161'],F0717_u['DM'],f_macs0717['lmass'],f_macs0717['lsfr'],f_macs0717['lssfr'],cat_macs0717['flux_radius'],cat_macs0717['star_flag'],cat_macs0717['use_KS'],cat_macs0717['use_phot']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','BCG','use_phot'))
abell370 = Table([z_abell370['id'],z_abell370['z_peak'],z_abell370['z_spec'],F370_u['L153'],F370_v['L155'],F370_j['L161'],F370_u['DM'],f_abell370['lmass'],f_abell370['lsfr'],f_abell370['lssfr'],cat_abell370['flux_radius'],cat_abell370['star_flag'],cat_abell370['use_KS'],cat_abell370['use_phot']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','BCG','use_phot'))
abell1063 = Table([z_abell1063['id'],z_abell1063['z_peak'],z_abell1063['z_spec'],F1063_u['L153'],F1063_v['L155'],F1063_j['L161'],F1063_u['DM'],f_abell1063['lmass'],f_abell1063['lsfr'],f_abell1063['lssfr'],cat_abell1063['flux_radius'],cat_abell1063['star_flag'],cat_abell1063['use_KS'],cat_abell1063['use_phot']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','BCG','use_phot'))
abell2744 = Table([z_abell2744['id'],z_abell2744['z_peak'],z_abell2744['z_spec'],F2744_u['L153'],F2744_v['L155'],F2744_j['L161'],F2744_u['DM'],f_abell2744['lmass'],f_abell2744['lsfr'],f_abell2744['lssfr'],cat_abell2744['flux_radius'],cat_abell2744['star_flag'],cat_abell2744['use_KS'],cat_abell2744['use_phot']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','BCG','use_phot'))
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
empty_1a = Column([99]*len(macs0416), name='sub', dtype=np.int8)    #identifier coumns ("sieves"). distinct sets for each cluster necessary as each individual catalogue has different number of targets in the field
empty_1b = Column([99]*len(macs0416), name='type', dtype=np.int8)
empty_1c = Column([99]*len(macs0416), name='member', dtype=np.int8)
empty_2a = Column([99]*len(macs1149), name='sub', dtype=np.int8)
empty_2b = Column([99]*len(macs1149), name='type', dtype=np.int8)
empty_2c = Column([99]*len(macs1149), name='member', dtype=np.int8)
empty_3a = Column([99]*len(macs0717), name='sub', dtype=np.int8)    #identifier coumns ("sieves"). distinct sets for each cluster necessary as each individual catalogue has different number of targets in the field
empty_3b = Column([99]*len(macs0717), name='type', dtype=np.int8)
empty_3c = Column([99]*len(macs0717), name='member', dtype=np.int8)
empty_4a = Column([99]*len(abell370), name='sub', dtype=np.int8)
empty_4b = Column([99]*len(abell370), name='type', dtype=np.int8)
empty_4c = Column([99]*len(abell370), name='member', dtype=np.int8)
empty_5a = Column([99]*len(abell1063), name='sub', dtype=np.int8)    #identifier coumns ("sieves"). distinct sets for each cluster necessary as each individual catalogue has different number of targets in the field
empty_5b = Column([99]*len(abell1063), name='type', dtype=np.int8)
empty_5c = Column([99]*len(abell1063), name='member', dtype=np.int8)
empty_6a = Column([99]*len(abell2744), name='sub', dtype=np.int8)
empty_6b = Column([99]*len(abell2744), name='type', dtype=np.int8)
empty_6c = Column([99]*len(abell2744), name='member', dtype=np.int8)

macs0416.add_columns([D1,empty_1a,empty_1b,empty_1c],[0,0,0,0])  
macs1149.add_columns([D2,empty_2a,empty_2b,empty_2c],[0,0,0,0])
macs0717.add_columns([D3,empty_3a,empty_3b,empty_3c],[0,0,0,0])
abell370.add_columns([D4,empty_4a,empty_4b,empty_4c],[0,0,0,0])
abell1063.add_columns([D5,empty_5a,empty_5b,empty_5c],[0,0,0,0])
abell2744.add_columns([D6,empty_6a,empty_6b,empty_6c],[0,0,0,0])
#
# sift for targets without no data, spec only, phot only, and both
spec_only = [[0]*6]
phot_only = [[0]*6]
error = [[0]*6] 
both = [[0]*6]
nodata = [[0]*6]
counter = 0
size = len(macs0416)
## Cluster 1: macs0416:
while counter < size:
    if macs0416[counter]['z_spec'] > 0 and macs0416[counter]['z_peak'] > 0:  #entries w/ both spectroscopy and photometry  (SPECTROSCOPIC sub-sample)
        if macs0416[counter]['use_phot'] ==0:
            spec_only[0][0]+=1
            macs0416[counter]['sub'] = 3
        elif macs0416[counter]['use_phot'] ==1:
            both[0][0] +=1
        else: error[0][0]+=1
        macs0416[counter]['sub'] = 1
    elif macs0416[counter]['z_spec'] > 0 and macs0416[counter]['z_peak'] < 0:  #entries w/ spectroscopy alone
        spec_only[0][0] +=1
        macs0416[counter]['sub'] = 3
    elif macs0416[counter]['z_spec'] < 0 and macs0416[counter]['z_peak'] > 0:  #entries w/ photometry alone (PHOTOMETRIC sub-sample)
        if macs0416[counter]['use_phot'] ==0:
            nodata[0][0]+=1
            macs0416[counter]['sub'] = 0
        elif macs0416[counter]['use_phot'] ==1:
            phot_only[0][0] +=1
            macs0416[counter]['sub'] = 2
    elif macs0416[counter]['z_spec'] < 0 and macs0416[counter]['z_peak'] < 0:  #entries w/ no z estimates at all
        nodata[0][0] +=1
        macs0416[counter]['sub'] = 0
    counter +=1
#
## Cluster 2: macs1149:
counter = 0
size = len(macs1149)
while counter < size:
    if macs1149[counter]['z_spec'] > 0 and macs1149[counter]['z_peak'] > 0:  #entries w/ both spectroscopy and photometry  (SPECTROSCOPIC sub-sample)
        if macs1149[counter]['use_phot'] ==0:
            spec_only[0][1]+=1
            macs1149[counter]['sub'] = 3
        elif macs1149[counter]['use_phot'] ==1:
            both[0][1] +=1
        else: error[0][1]+=1
        macs1149[counter]['sub'] = 1
    elif macs1149[counter]['z_spec'] > 0 and macs1149[counter]['z_peak'] < 0:  #entries w/ spectroscopy alone
        spec_only[0][1] +=1
        macs1149[counter]['sub'] = 3
    elif macs1149[counter]['z_spec'] < 0 and macs1149[counter]['z_peak'] > 0:  #entries w/ photometry alone (PHOTOMETRIC sub-sample)
        if macs1149[counter]['use_phot'] ==0:
            nodata[0][1]+=1
            macs1149[counter]['sub'] = 0
        elif macs1149[counter]['use_phot'] ==1:
            phot_only[0][1] +=1
            macs1149[counter]['sub'] = 2
    elif macs1149[counter]['z_spec'] < 0 and macs1149[counter]['z_peak'] < 0:  #entries w/ no z estimates at all
        nodata[0][1] +=1
        macs1149[counter]['sub'] = 0
    counter +=1
## Cluster 3: macs0717:
counter = 0
size = len(macs0717)
while counter < size:
    if macs0717[counter]['z_spec'] > 0 and macs0717[counter]['z_peak'] > 0:  #entries w/ both spectroscopy and photometry  (SPECTROSCOPIC sub-sample)
        if macs0717[counter]['use_phot'] ==0:
            spec_only[0][2]+=1
            macs0717[counter]['sub'] = 3
        elif macs0717[counter]['use_phot'] ==1:
            both[0][2] +=1
        else: error[0][2]+=1
        macs0717[counter]['sub'] = 1
    elif macs0717[counter]['z_spec'] > 0 and macs0717[counter]['z_peak'] < 0:  #entries w/ spectroscopy alone
        spec_only[0][2] +=1
        macs0717[counter]['sub'] = 3
    elif macs0717[counter]['z_spec'] < 0 and macs0717[counter]['z_peak'] > 0:  #entries w/ photometry alone (PHOTOMETRIC sub-sample)
        if macs0717[counter]['use_phot'] ==0:
            nodata[0][2]+=1
            macs0717[counter]['sub'] = 0
        elif macs0717[counter]['use_phot'] ==1:
            phot_only[0][2] +=1
            macs0717[counter]['sub'] = 2
    elif macs0717[counter]['z_spec'] < 0 and macs0717[counter]['z_peak'] < 0:  #entries w/ no z estimates at all
        nodata[0][2] +=1
        macs0717[counter]['sub'] = 0
    counter +=1
## Cluster 4: abell370:
counter = 0
size = len(abell370)
while counter < size:
    if abell370[counter]['z_spec'] > 0 and abell370[counter]['z_peak'] > 0:  #entries w/ both spectroscopy and photometry  (SPECTROSCOPIC sub-sample)
        if abell370[counter]['use_phot'] ==0:
            spec_only[0][3]+=1
            abell370[counter]['sub'] = 3
        elif abell370[counter]['use_phot'] ==1:
            both[0][3] +=1
        else: error[0][3]+=1
        abell370[counter]['sub'] = 1
    elif abell370[counter]['z_spec'] > 0 and abell370[counter]['z_peak'] < 0:  #entries w/ spectroscopy alone
        spec_only[0][3] +=1
        abell370[counter]['sub'] = 3
    elif abell370[counter]['z_spec'] < 0 and abell370[counter]['z_peak'] > 0:  #entries w/ photometry alone (PHOTOMETRIC sub-sample)
        if abell370[counter]['use_phot'] ==0:
            nodata[0][3]+=1
            abell370[counter]['sub'] = 0
        elif abell370[counter]['use_phot'] ==1:
            phot_only[0][3] +=1
            abell370[counter]['sub'] = 2
    elif abell370[counter]['z_spec'] < 0 and abell370[counter]['z_peak'] < 0:  #entries w/ no z estimates at all
        nodata[0][3] +=1
        abell370[counter]['sub'] = 0
    counter +=1
## Cluster 5: abell1063:
counter = 0
size = len(abell1063)
while counter < size:
    if abell1063[counter]['z_spec'] > 0 and abell1063[counter]['z_peak'] > 0:  #entries w/ both spectroscopy and photometry  (SPECTROSCOPIC sub-sample)
        if abell1063[counter]['use_phot'] ==0:
            spec_only[0][4]+=1
            abell1063[counter]['sub'] = 3
        elif abell1063[counter]['use_phot'] ==1:
            both[0][4] +=1
        else: error[0][4]+=1
        abell1063[counter]['sub'] = 1
    elif abell1063[counter]['z_spec'] > 0 and abell1063[counter]['z_peak'] < 0:  #entries w/ spectroscopy alone
        spec_only[0][4] +=1
        abell1063[counter]['sub'] = 3
    elif abell1063[counter]['z_spec'] < 0 and abell1063[counter]['z_peak'] > 0:  #entries w/ photometry alone (PHOTOMETRIC sub-sample)
        if abell1063[counter]['use_phot'] ==0:
            nodata[0][4]+=1
            abell1063[counter]['sub'] = 0
        elif abell1063[counter]['use_phot'] ==1:
            phot_only[0][4] +=1
            abell1063[counter]['sub'] = 2
    elif abell1063[counter]['z_spec'] < 0 and abell1063[counter]['z_peak'] < 0:  #entries w/ no z estimates at all
        nodata[0][4] +=1
        abell1063[counter]['sub'] = 0
    counter +=1
## Cluster 6: abell2744:
counter = 0
size = len(abell2744)
while counter < size:
    if abell2744[counter]['z_spec'] > 0 and abell2744[counter]['z_peak'] > 0:  #entries w/ both spectroscopy and photometry  (SPECTROSCOPIC sub-sample)
        if abell2744[counter]['use_phot'] ==0:
            spec_only[0][5]+=1
            abell2744[counter]['sub'] = 3
        elif abell2744[counter]['use_phot'] ==1:
            both[0][5] +=1
        else: error[0][5]+=1
        abell2744[counter]['sub'] = 1
    elif abell2744[counter]['z_spec'] > 0 and abell2744[counter]['z_peak'] < 0:  #entries w/ spectroscopy alone
        spec_only[0][5] +=1
        abell2744[counter]['sub'] = 3
    elif abell2744[counter]['z_spec'] < 0 and abell2744[counter]['z_peak'] > 0:  #entries w/ photometry alone (PHOTOMETRIC sub-sample)
        if abell2744[counter]['use_phot'] ==0:
            nodata[0][5]+=1
            abell2744[counter]['sub'] = 0
        elif abell2744[counter]['use_phot'] ==1:
            phot_only[0][5] +=1
            abell2744[counter]['sub'] = 2
    elif abell2744[counter]['z_spec'] < 0 and abell2744[counter]['z_peak'] < 0:  #entries w/ no z estimates at all
        nodata[0][5] +=1
        abell2744[counter]['sub'] = 0
    counter +=1
## Summarize initial data stats in table
data_names = Column(['total','spec & phot','only phot','spec only','no data'],name='Property')
data0 = Column([np.sum([len(cat_macs0416),len(cat_macs1149),len(cat_macs0717),len(cat_abell370),len(cat_abell1063),len(cat_abell2744)]),np.sum([both]),np.sum([phot_only]),np.sum([spec_only]),np.sum([nodata])],name='Total')
data1 = Column([len(cat_macs0416),both[0][0],phot_only[0][0],spec_only[0][0],nodata[0][0]],name='macs0416')
data2 = Column([len(cat_macs1149),both[0][1],phot_only[0][1],spec_only[0][1],nodata[0][1]],name='macs1149')
data3 = Column([len(cat_macs0717),both[0][2],phot_only[0][2],spec_only[0][2],nodata[0][2]],name='macs0717')
data4 = Column([len(cat_abell370),both[0][3],phot_only[0][3],spec_only[0][3],nodata[0][3]],name='abell370')
data5 = Column([len(cat_abell1063),both[0][4],phot_only[0][4],spec_only[0][4],nodata[0][4]],name='abell1063')
data6 = Column([len(cat_abell2744),both[0][5],phot_only[0][5],spec_only[0][5],nodata[0][5]],name='abell2744')
global data
data_stats = Table([data_names,data0,data1,data2,data3,data4,data5,data6])  
#
##add columns for luminosity calculations to each cluster catalogue individually
##for some reason there are problems with simply doing this once to the final master_cat(alogue)
# 1.macs0416
empty_1u = Column([99]*len(macs0416), name='L_u', dtype=np.float64)
empty_1v = Column([99]*len(macs0416), name='L_v', dtype=np.float64)
empty_1j = Column([99]*len(macs0416), name='L_j', dtype=np.float64)
empty_1uv = Column([99]*len(macs0416), name='uv', dtype=np.float64)
empty_1vj = Column([99]*len(macs0416), name='vj', dtype=np.float64)
macs0416.add_columns([empty_1u,empty_1v,empty_1j,empty_1uv,empty_1vj])
# 2.macs1149  
empty_2u = Column([99]*len(macs1149), name='L_u', dtype=np.float64)
empty_2v = Column([99]*len(macs1149), name='L_v', dtype=np.float64)
empty_2j = Column([99]*len(macs1149), name='L_j', dtype=np.float64)
empty_2uv = Column([99]*len(macs1149), name='uv', dtype=np.float64)
empty_2vj = Column([99]*len(macs1149), name='vj', dtype=np.float64)
macs1149.add_columns([empty_2u,empty_2v,empty_2j,empty_2uv,empty_2vj])
# 3.macs0717  
empty_3u = Column([99]*len(macs0717), name='L_u', dtype=np.float64)
empty_3v = Column([99]*len(macs0717), name='L_v', dtype=np.float64)
empty_3j = Column([99]*len(macs0717), name='L_j', dtype=np.float64)
empty_3uv = Column([99]*len(macs0717), name='uv', dtype=np.float64)
empty_3vj = Column([99]*len(macs0717), name='vj', dtype=np.float64)
macs0717.add_columns([empty_3u,empty_3v,empty_3j,empty_3uv,empty_3vj])
# 4.abell370  
empty_4u = Column([99]*len(abell370), name='L_u', dtype=np.float64)
empty_4v = Column([99]*len(abell370), name='L_v', dtype=np.float64)
empty_4j = Column([99]*len(abell370), name='L_j', dtype=np.float64)
empty_4uv = Column([99]*len(abell370), name='uv', dtype=np.float64)
empty_4vj = Column([99]*len(abell370), name='vj', dtype=np.float64)
abell370.add_columns([empty_4u,empty_4v,empty_4j,empty_4uv,empty_4vj])
# 5.abell1063  
empty_5u = Column([99]*len(abell1063), name='L_u', dtype=np.float64)
empty_5v = Column([99]*len(abell1063), name='L_v', dtype=np.float64)
empty_5j = Column([99]*len(abell1063), name='L_j', dtype=np.float64)
empty_5uv = Column([99]*len(abell1063), name='uv', dtype=np.float64)
empty_5vj = Column([99]*len(abell1063), name='vj', dtype=np.float64)
abell1063.add_columns([empty_5u,empty_5v,empty_5j,empty_5uv,empty_5vj])
# 6.abell2744
empty_6u = Column([99]*len(abell2744), name='L_u', dtype=np.float64)
empty_6v = Column([99]*len(abell2744), name='L_v', dtype=np.float64)
empty_6j = Column([99]*len(abell2744), name='L_j', dtype=np.float64)
empty_6uv = Column([99]*len(abell2744), name='uv', dtype=np.float64)
empty_6vj = Column([99]*len(abell2744), name='vj', dtype=np.float64)
abell2744.add_columns([empty_6u,empty_6v,empty_6j,empty_6uv,empty_6vj])
#
#
##convert flux to magnitude (erroneously labelled as luminosity, e.g. L_u for magnitude in UV), get color indices U-V, V-J, add to table
#
####1. macs0416
counter = 0
size = len(macs0416)
while counter < size:
    if macs0416[counter]['sub'] ==0:
        counter +=1 
    else:
        macs0416[counter]['L_u'] = -2.5*np.log10(macs0416[counter]['u']) + 25
        macs0416[counter]['L_v'] = -2.5*np.log10(macs0416[counter]['v']) + 25
        macs0416[counter]['L_j'] = -2.5*np.log10(macs0416[counter]['j']) + 25
        macs0416[counter]['uv'] = macs0416[counter]['L_u'] - macs0416[counter]['L_v']
        macs0416[counter]['vj'] = macs0416[counter]['L_v'] - macs0416[counter]['L_j']
        counter +=1
#
#####2. macs1149
      #cluster designation column
counter = 0
size = len(macs1149)
while counter < size:
    if macs1149[counter]['sub'] ==0:
        counter +=1 
    else:
        macs1149[counter]['L_u'] = -2.5*np.log10(macs1149[counter]['u']) + 25
        macs1149[counter]['L_v'] = -2.5*np.log10(macs1149[counter]['v']) + 25
        macs1149[counter]['L_j'] = -2.5*np.log10(macs1149[counter]['j']) + 25
        macs1149[counter]['uv'] = macs1149[counter]['L_u'] - macs1149[counter]['L_v']
        macs1149[counter]['vj'] = macs1149[counter]['L_v'] - macs1149[counter]['L_j']
        counter +=1
#####3.macs0717
aa = 0
counter = 0
size = len(macs0717)
while counter < size:
    if macs0717[counter]['sub'] ==0 or macs0717[counter]['sub'] ==3:
        aa+=1
        counter +=1 
    else:
        macs0717[counter]['L_u'] = -2.5*np.log10(macs0717[counter]['u']) + 25
        macs0717[counter]['L_v'] = -2.5*np.log10(macs0717[counter]['v']) + 25
        macs0717[counter]['L_j'] = -2.5*np.log10(macs0717[counter]['j']) + 25
        macs0717[counter]['uv'] = macs0717[counter]['L_u'] - macs0717[counter]['L_v']
        macs0717[counter]['vj'] = macs0717[counter]['L_v'] - macs0717[counter]['L_j']
        counter +=1
#####4.abell370
counter = 0
size = len(abell370)
while counter < size:
    if abell370[counter]['sub'] ==0:
        counter +=1 
    else:
        abell370[counter]['L_u'] = -2.5*np.log10(abell370[counter]['u']) + 25
        abell370[counter]['L_v'] = -2.5*np.log10(abell370[counter]['v']) + 25
        abell370[counter]['L_j'] = -2.5*np.log10(abell370[counter]['j']) + 25
        abell370[counter]['uv'] = abell370[counter]['L_u'] - abell370[counter]['L_v']
        abell370[counter]['vj'] = abell370[counter]['L_v'] - abell370[counter]['L_j']
        counter +=1
#####5.abell1063
counter = 0
size = len(abell1063)
while counter < size:
    if abell1063[counter]['sub'] ==0:
        counter +=1 
    else:
        abell1063[counter]['L_u'] = -2.5*np.log10(abell1063[counter]['u']) + 25
        abell1063[counter]['L_v'] = -2.5*np.log10(abell1063[counter]['v']) + 25
        abell1063[counter]['L_j'] = -2.5*np.log10(abell1063[counter]['j']) + 25
        abell1063[counter]['uv'] = abell1063[counter]['L_u'] - abell1063[counter]['L_v']
        abell1063[counter]['vj'] = abell1063[counter]['L_v'] - abell1063[counter]['L_j']
        counter +=1
#####6.abell2744
counter = 0
size = len(abell2744)
while counter < size:
    if abell2744[counter]['sub'] ==0:
        counter +=1 
    else:
        abell2744[counter]['L_u'] = -2.5*np.log10(abell2744[counter]['u']) + 25
        abell2744[counter]['L_v'] = -2.5*np.log10(abell2744[counter]['v']) + 25
        abell2744[counter]['L_j'] = -2.5*np.log10(abell2744[counter]['j']) + 25
        abell2744[counter]['uv'] = abell2744[counter]['L_u'] - abell2744[counter]['L_v']
        abell2744[counter]['vj'] = abell2744[counter]['L_v'] - abell2744[counter]['L_j']
        counter +=1
#
global master_cat
master_cat = Table(np.concatenate((macs0416,macs1149,macs0717,abell370,abell1063,abell2744), axis=0))  #create a master catalogue of all clusters
#
#
### **Section 2: calulate z's & separate outliers
## (i) calculate delta_z for spec subsample; (ii) separate outliers from 
## both phot & spec sub-sample;  
#
#####Note: the photometric redshift used is 'z_peak' from data
#
#(i) calculate delta_z for targets w/ both z_spec & z_phot
master_cat.add_column(Column([-99]*len(master_cat),name='del_z', dtype=np.float64))                    #<-- z_phot - z_spec / 1 + z_spec
master_cat.add_column(Column([-99]*len(master_cat),name='z_clusterspec', dtype=np.float64))   #<-- z_spec - z_cl / 1 + z_spec
master_cat.add_column(Column([-99]*len(master_cat),name='z_clusterphot', dtype=np.float64))   #<-- z_phot - z_cl / 1 + z_phot
#
counter = 0
size = len(master_cat)
while counter < size:                       #calucalte del_z for outlier cut
    if master_cat[counter]['sub'] == 1:
        master_cat[counter]['del_z'] = ((master_cat[counter]['z_peak'] - master_cat[counter]['z_spec']) / (1 + master_cat[counter]['z_spec']))
        counter +=1
    else:
        counter +=1
#
# (ii) separate outliers for |del_z| > 0.15
## "member" == 4 is outliers;
outliers = 0
out = [[0]*6]
spec_subsample = 0
phot_subsample = 0
size = len(master_cat)
counter = 0
while counter < size: 
    if master_cat[counter]['sub'] == 2 and master_cat[counter]['sub'] != 0:
        phot_subsample +=1
    elif master_cat[counter]['sub'] == 1 and master_cat[counter]['sub'] != 0:
        if abs(master_cat[counter]['del_z']) > 0.15: #identify objects w/ del_z > 0.15
            master_cat[counter]['type'] = 3
            outliers += 1
            if master_cat[counter]['cluster'] == 1:
                    out[0][0] +=1
                    both[0][0] -=1
            elif master_cat[counter]['cluster'] == 2:
                    out[0][1] +=1
                    both[0][1] -=1
            elif master_cat[counter]['cluster'] == 3:
                    out[0][2] +=1
                    both[0][2] -=1
            elif master_cat[counter]['cluster'] == 4:
                    out[0][3] +=1
                    both[0][3] -=1
            elif master_cat[counter]['cluster'] == 5:
                    out[0][4] +=1
                    both[0][4] -=1
            elif master_cat[counter]['cluster'] == 6:
                    out[0][5] +=1
                    both[0][5] -=1
        else: 
            spec_subsample +=1          ##spec_subsample includes sub==1 & sub==3
    counter +=1
#
## store cluster redshifts; obtained from https://archive.stsci.edu/prepds/frontier/
z_cluster = [0.396,0.543,0.545,0.375,0.348,0.308]      
#  [macs0416,macs1149,macs0717,abell370,abell1063,abell2744]   <- order of elements in cluster array, i.e. z_macs0717 = z_cluster[2] = 0.545
#
## calculate z_clusters for del_zphot vs del_zspec plot, defined below. these 
## will be used to make cuts to spec sample, from which we will use relative 
## fractions by mass & radial distance to correct the photometric sample for 
## completeness.
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['sub'] == 1:         #spec sample del_z plot for cuts: denominator = (1 + z_spec)
        if master_cat[counter]['cluster'] == 1:
            master_cat[counter]['z_clusterspec'] = ((master_cat[counter]['z_spec'] - z_cluster[0]) / (1 + master_cat[counter]['z_spec']))
            master_cat[counter]['z_clusterphot'] = ((master_cat[counter]['z_peak'] - z_cluster[0]) / (1 + master_cat[counter]['z_peak']))
        elif master_cat[counter]['cluster'] == 2:
            master_cat[counter]['z_clusterspec'] = ((master_cat[counter]['z_spec'] - z_cluster[1]) / (1 + master_cat[counter]['z_spec']))
            master_cat[counter]['z_clusterphot'] = ((master_cat[counter]['z_peak'] - z_cluster[1]) / (1 + master_cat[counter]['z_peak']))
        elif master_cat[counter]['cluster'] == 3:
            master_cat[counter]['z_clusterspec'] = (master_cat[counter]['z_spec'] - z_cluster[2]) / (1 + master_cat[counter]['z_spec'])
            master_cat[counter]['z_clusterphot'] = (master_cat[counter]['z_peak'] - z_cluster[2]) / (1 + master_cat[counter]['z_peak'])
        elif master_cat[counter]['cluster'] == 4:
            master_cat[counter]['z_clusterspec'] = (master_cat[counter]['z_spec'] - z_cluster[3]) / (1 + master_cat[counter]['z_spec'])
            master_cat[counter]['z_clusterphot'] = (master_cat[counter]['z_peak'] - z_cluster[3]) / (1 + master_cat[counter]['z_peak'])
        elif master_cat[counter]['cluster'] == 5:
            master_cat[counter]['z_clusterspec'] = (master_cat[counter]['z_spec'] - z_cluster[4]) / (1 + master_cat[counter]['z_spec'])
            master_cat[counter]['z_clusterphot'] = (master_cat[counter]['z_peak'] - z_cluster[4]) / (1 + master_cat[counter]['z_peak'])
        elif master_cat[counter]['cluster'] == 6:
            master_cat[counter]['z_clusterspec'] = (master_cat[counter]['z_spec'] - z_cluster[5]) / (1 + master_cat[counter]['z_spec'])
            master_cat[counter]['z_clusterphot'] = (master_cat[counter]['z_peak'] - z_cluster[5]) / (1 + master_cat[counter]['z_peak'])
        counter +=1
    elif master_cat[counter]['sub'] == 2:       #phot sample del_z plot for cuts defined in v7a to be (z_phot - z_cl)
        if master_cat[counter]['cluster'] == 1:
            master_cat[counter]['z_clusterphot'] = (master_cat[counter]['z_peak'] - z_cluster[0]) / (1 + z_cluster[0])
        elif master_cat[counter]['cluster'] == 2:
            master_cat[counter]['z_clusterphot'] = (master_cat[counter]['z_peak'] - z_cluster[1]) / (1 + z_cluster[1])
        elif master_cat[counter]['cluster'] == 3:
            master_cat[counter]['z_clusterphot'] = (master_cat[counter]['z_peak'] - z_cluster[2]) / (1 + z_cluster[2])
        elif master_cat[counter]['cluster'] == 4:
            master_cat[counter]['z_clusterphot'] = (master_cat[counter]['z_peak'] - z_cluster[3]) / (1 + z_cluster[3])
        elif master_cat[counter]['cluster'] == 5:
            master_cat[counter]['z_clusterphot'] = (master_cat[counter]['z_peak'] - z_cluster[4]) / (1 + z_cluster[4])
        elif master_cat[counter]['cluster'] == 6:
            master_cat[counter]['z_clusterphot'] = (master_cat[counter]['z_peak'] - z_cluster[5]) / (1 + z_cluster[5]) 
        counter +=1
    else: 
        counter +=1
#
## Section 3: separate SF/Q for both subsamples based on van der Burg (2013) 
## colour criteria;    filter name: 'type'  IDs below
##   0 = stars;  1 = SF (star-forming);    2 = Q (quiscient);    3 = outliers
#
##(i) sift for stars in “star_flag” from catalogue; 0: galaxies, 1: stars, 2: fainter than 25 magnitude
star = [[0]*6]
star_type = [[0]*3]
size = len(master_cat)
counter = 0
while counter < size: 
    if master_cat[counter]['star_flag'] == 1: #identify objects with a "star flag" = 1 for stars; gal=0, mag>25 = 2
        master_cat[counter]['type'] = 0
        if master_cat[counter]['sub'] == 0:     #keep track of which subsample stars came from
            star_type[0][0] +=1
        elif master_cat[counter]['sub'] == 1:   #spec stars
            spec_subsample -=1
            star_type[0][1] +=1
        elif master_cat[counter]['sub'] ==2:    #phot stars
            phot_subsample -=1
            star_type[0][2] +=1
        if master_cat[counter]['cluster'] == 1:
            star[0][0] +=1
        elif master_cat[counter]['cluster'] == 2:
            star[0][1] +=1
        elif master_cat[counter]['cluster'] == 3:
            star[0][2] +=1    
        elif master_cat[counter]['cluster'] == 4:
            star[0][3] +=1    
        elif master_cat[counter]['cluster'] == 5:
            star[0][4] +=1    
        elif master_cat[counter]['cluster'] == 6:
            star[0][5] +=1
        counter +=1
    else:
        counter +=1
#
## SF/Q segragation
spec_SF = 0
spec_Q = 0
phot_SF = 0
phot_Q = 0
#stat counting variables
SF_spec = [[0]*6]
SF_phot = [[0]*6]
Q_spec = [[0]*6]
Q_phot = [[0]*6]
check1 = 0  #check for stars
check2 = 0  #check that rows NOT counted should equal sum(nodata,out,spec_only)
size = len(master_cat)
counter = 0
while counter < size: 
    if master_cat[counter]['sub'] == 1 and master_cat[counter]['type'] != 3 and master_cat[counter]['type'] != 0:        #spec subsample, skipping stars & outliers
        if master_cat[counter]['uv'] > 1.3 and master_cat[counter]['vj'] < 1.6 and master_cat[counter]['uv'] > ((0.88*master_cat[counter]['vj']) + 0.6): #identify passive (quiescent)
            master_cat[counter]['type'] = 2
#            spec_Q +=1
            if master_cat[counter]['cluster'] == 1:
                Q_spec[0][0] +=1
            elif master_cat[counter]['cluster'] == 2:
                Q_spec[0][1] +=1
            elif master_cat[counter]['cluster'] == 3:
                Q_spec[0][2] +=1    
            elif master_cat[counter]['cluster'] == 4:
                Q_spec[0][3] +=1    
            elif master_cat[counter]['cluster'] == 5:
                Q_spec[0][4] +=1    
            elif master_cat[counter]['cluster'] == 6:
                Q_spec[0][5] +=1   
        else:
            master_cat[counter]['type'] = 1         #if not Q, then must be SF
#            spec_SF +=1
            if master_cat[counter]['cluster'] == 1:
                SF_spec[0][0] +=1
            elif master_cat[counter]['cluster'] == 2:
                SF_spec[0][1] +=1
            elif master_cat[counter]['cluster'] == 3:
                SF_spec[0][2] +=1    
            elif master_cat[counter]['cluster'] == 4:
                SF_spec[0][3] +=1    
            elif master_cat[counter]['cluster'] == 5:
                SF_spec[0][4] +=1    
            elif master_cat[counter]['cluster'] == 6:
                SF_spec[0][5] +=1  
    elif master_cat[counter]['sub'] == 2 and master_cat[counter]['type'] != 3 and master_cat[counter]['type'] != 0:               #phot subsample
        if master_cat[counter]['uv'] > 1.3 and master_cat[counter]['vj'] < 1.6 and master_cat[counter]['uv'] > ((0.88*master_cat[counter]['vj'] + 0.6)): #identify passive (quiescent)
            master_cat[counter]['type'] = 2
#            phot_Q +=1
            if master_cat[counter]['cluster'] == 1:
                Q_phot[0][0] +=1
            elif master_cat[counter]['cluster'] == 2:
                Q_phot[0][1] +=1
            elif master_cat[counter]['cluster'] == 3:
                Q_phot[0][2] +=1    
            elif master_cat[counter]['cluster'] == 4:
                Q_phot[0][3] +=1    
            elif master_cat[counter]['cluster'] == 5:
                Q_phot[0][4] +=1    
            elif master_cat[counter]['cluster'] == 6:
                Q_phot[0][5] +=1     
        else:
            master_cat[counter]['type'] = 1
#            phot_SF +=1
            if master_cat[counter]['cluster'] == 1:
                SF_phot[0][0] +=1
            elif master_cat[counter]['cluster'] == 2:
                SF_phot[0][1] +=1
            elif master_cat[counter]['cluster'] == 3:
                SF_phot[0][2] +=1    
            elif master_cat[counter]['cluster'] == 4:
                SF_phot[0][3] +=1    
            elif master_cat[counter]['cluster'] == 5:
                SF_phot[0][4] +=1    
            elif master_cat[counter]['cluster'] == 6:
                SF_phot[0][5] +=1   
    elif master_cat[counter]['type'] == 0:               #stars
        check1 +=1
    else: 
        check2 +=1      #check that entries not listed above are outliers & nodata ONLY i.e. check should equal np.sum(nodata)+np.sum(out)
    counter +=1
#
#
macro_type = Column(['total','SF','Q','stars','outliers'])
tot = Column([np.sum(both)-np.sum(spec_only)+outliers+star_type[0][1]+np.sum(phot_only),np.sum(SF_spec)+np.sum(SF_phot),np.sum(Q_spec)+np.sum(Q_phot),np.sum(star_type),outliers])
sp = Column([np.sum(both)-np.sum(spec_only)+outliers+star_type[0][1],np.sum(SF_spec),np.sum(Q_spec),star_type[0][0],outliers])
ph = Column([np.sum(phot_only),np.sum(SF_phot),np.sum(Q_phot),star_type[0][1],0])
macros_type = Table([macro_type,tot,sp,ph],names=('Property','Total','spec','phot'))
#
#
##   **Section 4: apply delta_z cuts to determine membership of z_phot sample based on cut
## criteria discussed w/ AM: |del_z phot| < 0.05 & |del_z spec| < 0.01, in the process isolate 
## secure field, false pos/neg items from. ID column to identify populations is
## called "member"         
## separate  secure field/member & false pos/neg using criteria: 
## cluster = 0: abs(z_clusterspec) < 0.01 & abs(z_cluster phot) < 0.05; 
## field = 1: abs(z_clusterspec) > 0.01 & abs(z_cluster phot) > 0.05; 
## false pos = 2: abs(z_clusterspec) > 0.01 & abs(z_cluster phot) < 0.05;
## false neg = 3: abs(z_clusterspec) < 0.01 & abs(z_cluster phot) > 0.05; 
## Note: this is only done for spec sample; results used to correct photo sample for completeness at end of analysis
#
### SF sub-sample: SF_spec into secure member, field, false pos/neg
###
#
smem = [[0]*6]
qmem = [[0]*6]
sfield = [[0]*6]
qfield = [[0]*6]
spos = [[0]*6]
qpos = [[0]*6]
sneg = [[0]*6]
qneg = [[0]*6]
other = 0
size = len(master_cat)
counter = 0
while counter < size:
    if master_cat[counter]['sub'] == 1:       # ['sub'] == 1: spec
        if master_cat[counter]['type'] == 0 or master_cat[counter]['type'] ==3:    #['type'] == 0: stars,    ['type'] == 3: spec_only
            other +=1
        elif master_cat[counter]['type'] == 1:      # '[type'] == 1:   look only at SF spec population
            if abs(master_cat[counter]['z_clusterspec']) > 0.01 and abs(master_cat[counter]['z_clusterphot']) > 0.05:  # secure field
                master_cat[counter]['member'] = 1   #field
                if master_cat[counter]['cluster'] == 1:
                    sfield[0][0] +=1
                elif master_cat[counter]['cluster'] == 2:
                    sfield[0][1] +=1
                elif master_cat[counter]['cluster'] == 3:
                    sfield[0][2] +=1    
                elif master_cat[counter]['cluster'] == 4:
                    sfield[0][3] +=1    
                elif master_cat[counter]['cluster'] == 5:
                    sfield[0][4] +=1    
                elif master_cat[counter]['cluster'] == 6:
                    sfield[0][5] +=1  
            elif abs(master_cat[counter]['z_clusterspec']) > 0.01 and abs(master_cat[counter]['z_clusterphot']) < 0.05:  #  false positive
                master_cat[counter]['member'] = 2   #false positive
                if master_cat[counter]['cluster'] == 1:
                    spos[0][0] +=1
                elif master_cat[counter]['cluster'] == 2:
                    spos[0][1] +=1
                elif master_cat[counter]['cluster'] == 3:
                    spos[0][2] +=1    
                elif master_cat[counter]['cluster'] == 4:
                    spos[0][3] +=1    
                elif master_cat[counter]['cluster'] == 5:
                    spos[0][4] +=1    
                elif master_cat[counter]['cluster'] == 6:
                    spos[0][5] +=1
            elif abs(master_cat[counter]['z_clusterspec']) < 0.01 and abs(master_cat[counter]['z_clusterphot']) > 0.05:  #  false negative
                master_cat[counter]['member'] = 3   #false negative
                if master_cat[counter]['cluster'] == 1:
                    sneg[0][0] +=1
                elif master_cat[counter]['cluster'] == 2:
                    sneg[0][1] +=1
                elif master_cat[counter]['cluster'] == 3:
                    sneg[0][2] +=1    
                elif master_cat[counter]['cluster'] == 4:
                    sneg[0][3] +=1    
                elif master_cat[counter]['cluster'] == 5:
                    sneg[0][4] +=1    
                elif master_cat[counter]['cluster'] == 6:
                    sneg[0][5] +=1
            else:
                master_cat[counter]['member'] = 0              #secure cluster member
                if master_cat[counter]['cluster'] == 1:
                    smem[0][0] +=1
                elif master_cat[counter]['cluster'] == 2:
                    smem[0][1] +=1
                elif master_cat[counter]['cluster'] == 3:
                    smem[0][2] +=1    
                elif master_cat[counter]['cluster'] == 4:
                    smem[0][3] +=1    
                elif master_cat[counter]['cluster'] == 5:
                    smem[0][4] +=1    
                elif master_cat[counter]['cluster'] == 6:
                    smem[0][5] +=1 
        elif master_cat[counter]['type'] == 2:          
            if abs(master_cat[counter]['z_clusterspec']) > 0.01 and abs(master_cat[counter]['z_clusterphot']) > 0.05:  # secure field
                master_cat[counter]['member'] = 1       #field
                if master_cat[counter]['cluster'] == 1:
                    qfield[0][0] +=1
                elif master_cat[counter]['cluster'] == 2:
                    qfield[0][1] +=1
                elif master_cat[counter]['cluster'] == 3:
                    qfield[0][2] +=1    
                elif master_cat[counter]['cluster'] == 4:
                    qfield[0][3] +=1    
                elif master_cat[counter]['cluster'] == 5:
                    qfield[0][4] +=1    
                elif master_cat[counter]['cluster'] == 6:
                    qfield[0][5] +=1  
            elif abs(master_cat[counter]['z_clusterspec']) > 0.01 and abs(master_cat[counter]['z_clusterphot']) < 0.05:  #  false positive
                master_cat[counter]['member'] = 2       #false positive
                if master_cat[counter]['cluster'] == 1:
                    qpos[0][0] +=1
                elif master_cat[counter]['cluster'] == 2:
                    qpos[0][1] +=1
                elif master_cat[counter]['cluster'] == 3:
                    qpos[0][2] +=1    
                elif master_cat[counter]['cluster'] == 4:
                    qpos[0][3] +=1    
                elif master_cat[counter]['cluster'] == 5:
                    qpos[0][4] +=1    
                elif master_cat[counter]['cluster'] == 6:
                    qpos[0][5] +=1
            elif abs(master_cat[counter]['z_clusterspec']) < 0.01 and abs(master_cat[counter]['z_clusterphot']) > 0.05:  #  false negative
                master_cat[counter]['member'] = 3       #false negative
                if master_cat[counter]['cluster'] == 1:
                    qneg[0][0] +=1
                elif master_cat[counter]['cluster'] == 2:
                    qneg[0][1] +=1
                elif master_cat[counter]['cluster'] == 3:
                    qneg[0][2] +=1    
                elif master_cat[counter]['cluster'] == 4:
                    qneg[0][3] +=1    
                elif master_cat[counter]['cluster'] == 5:
                    qneg[0][4] +=1    
                elif master_cat[counter]['cluster'] == 6:
                    qneg[0][5] +=1
            else:
                master_cat[counter]['member'] = 0              #secure cluster member
                if master_cat[counter]['cluster'] == 1:
                    qmem[0][0] +=1
                elif master_cat[counter]['cluster'] == 2:
                    qmem[0][1] +=1
                elif master_cat[counter]['cluster'] == 3:
                    qmem[0][2] +=1    
                elif master_cat[counter]['cluster'] == 4:
                    qmem[0][3] +=1    
                elif master_cat[counter]['cluster'] == 5:
                    qmem[0][4] +=1    
                elif master_cat[counter]['cluster'] == 6:
                    qmem[0][5] +=1
        counter +=1
    else:     
        counter +=1
#
## summarize key macro information by cluster in a table
#
macro_names = Column(['total','no data','spec total(no out)','phot total','stars','SF total','SF spec','SF phot','Q total','Q spec','Q phot','outliers (spec)'])
total = Column([len(master_cat),np.sum(nodata),spec_subsample,phot_subsample,np.sum(star),np.sum(SF_spec)+np.sum(SF_phot),np.sum(SF_spec),np.sum(SF_phot),np.sum(Q_spec)+np.sum(Q_phot),np.sum(Q_spec),np.sum(Q_phot),outliers])
cl1 = Column([len(cat_macs0416),nodata[0][0],both[0][0]+spec_only[0][0],phot_only[0][0],star[0][0],SF_spec[0][0]+SF_phot[0][0],SF_spec[0][0],SF_phot[0][0],Q_spec[0][0]+Q_phot[0][0],Q_spec[0][0],Q_phot[0][0],out[0][0]])
cl2 = Column([len(cat_macs1149),nodata[0][1],both[0][1]+spec_only[0][1],phot_only[0][1],star[0][1],SF_spec[0][1]+SF_phot[0][1],SF_spec[0][1],SF_phot[0][1],Q_spec[0][1]+Q_phot[0][1],Q_spec[0][1],Q_phot[0][1],out[0][1]])
cl3 = Column([len(cat_macs0717),nodata[0][2],both[0][2]+spec_only[0][2],phot_only[0][2],star[0][2],SF_spec[0][2]+SF_phot[0][2],SF_spec[0][2],SF_phot[0][2],Q_spec[0][2]+Q_phot[0][2],Q_spec[0][2],Q_phot[0][2],out[0][2]])
cl4 = Column([len(cat_abell370),nodata[0][3],both[0][3]+spec_only[0][3],phot_only[0][3],star[0][3],SF_spec[0][3]+SF_phot[0][3],SF_spec[0][3],SF_phot[0][3],Q_spec[0][3]+Q_phot[0][3],Q_spec[0][3],Q_phot[0][3],out[0][3]])
cl5 = Column([len(cat_abell1063),nodata[0][4],both[0][4]+spec_only[0][4],phot_only[0][4],star[0][4],SF_spec[0][4]+SF_phot[0][4],SF_spec[0][4],SF_phot[0][4],Q_spec[0][4]+Q_phot[0][4],Q_spec[0][4],Q_phot[0][4],out[0][4]])
cl6 = Column([len(cat_abell2744),nodata[0][5],both[0][5]+spec_only[0][5],phot_only[0][5],star[0][5],SF_spec[0][5]+SF_phot[0][5],SF_spec[0][5],SF_phot[0][5],Q_spec[0][5]+Q_phot[0][5],Q_spec[0][5],Q_phot[0][5],out[0][5]])
macros1 = Table([macro_names,total,cl1,cl2],names=('Property','Full catalogue','1.macs0416','2.macs1149'))
macros2 = Table([macro_names,cl3,cl4,cl5,cl6],names=('Property','3.macs0717','4.abell370','5.abell1063','6.abel2744'))
#
##summarize spec sample stats in table
#
spec_names = Column(['secure cluster - SF','secure field - SF','false positive - SF','false negative - SF','total SF','secure cluster - Q','secure field - Q','false positive - Q','false negative - Q','total Q','Total'])
total_mems = Column([np.sum(smem),np.sum(sfield),np.sum(spos),np.sum(sneg),np.sum([smem,sfield,spos,sneg]),np.sum(qmem),np.sum(qfield),np.sum(qpos),np.sum(qneg),np.sum([qmem,qfield,qpos,qneg]),np.sum([smem,sfield,spos,sneg,qmem,qfield,qpos,qneg])])
cl1_mems = Column([smem[0][0],sfield[0][0],spos[0][0],sneg[0][0],np.sum([smem[0][0],sfield[0][0],spos[0][0],sneg[0][0]]),qmem[0][0],qfield[0][0],qpos[0][0],qneg[0][0],np.sum([qmem[0][0],qfield[0][0],qpos[0][0],qneg[0][0]]),np.sum([smem[0][0],sfield[0][0],spos[0][0],sneg[0][0],qmem[0][0],qfield[0][0],qpos[0][0],qneg[0][0]])])
cl2_mems = Column([smem[0][1],sfield[0][1],spos[0][1],sneg[0][1],np.sum([smem[0][1],sfield[0][1],spos[0][1],sneg[0][1]]),qmem[0][1],qfield[0][1],qpos[0][1],qneg[0][1],np.sum([qmem[0][1],qfield[0][1],qpos[0][1],qneg[0][1]]),np.sum([smem[0][1],sfield[0][1],spos[0][1],sneg[0][1],qmem[0][1],qfield[0][1],qpos[0][1],qneg[0][1]])])
cl3_mems = Column([smem[0][2],sfield[0][2],spos[0][2],sneg[0][2],np.sum([smem[0][2],sfield[0][2],spos[0][2],sneg[0][2]]),qmem[0][2],qfield[0][2],qpos[0][2],qneg[0][2],np.sum([qmem[0][2],qfield[0][2],qpos[0][2],qneg[0][2]]),np.sum([smem[0][2],sfield[0][2],spos[0][2],sneg[0][2],qmem[0][2],qfield[0][2],qpos[0][2],qneg[0][2]])])
cl4_mems = Column([smem[0][3],sfield[0][3],spos[0][3],sneg[0][3],np.sum([smem[0][3],sfield[0][3],spos[0][3],sneg[0][3]]),qmem[0][3],qfield[0][3],qpos[0][3],qneg[0][3],np.sum([qmem[0][3],qfield[0][3],qpos[0][3],qneg[0][3]]),np.sum([smem[0][3],sfield[0][3],spos[0][3],sneg[0][3],qmem[0][3],qfield[0][3],qpos[0][3],qneg[0][3]])])
cl5_mems = Column([smem[0][4],sfield[0][4],spos[0][4],sneg[0][4],np.sum([smem[0][4],sfield[0][4],spos[0][4],sneg[0][4]]),qmem[0][4],qfield[0][4],qpos[0][4],qneg[0][4],np.sum([qmem[0][4],qfield[0][4],qpos[0][4],qneg[0][4]]),np.sum([smem[0][4],sfield[0][4],spos[0][4],sneg[0][4],qmem[0][4],qfield[0][4],qpos[0][4],qneg[0][4]])])
cl6_mems = Column([smem[0][5],sfield[0][5],spos[0][5],sneg[0][5],np.sum([smem[0][5],sfield[0][5],spos[0][5],sneg[0][5]]),qmem[0][5],qfield[0][5],qpos[0][5],qneg[0][5],np.sum([qmem[0][5],qfield[0][5],qpos[0][5],qneg[0][5]]),np.sum([smem[0][5],sfield[0][5],spos[0][5],sneg[0][5],qmem[0][5],qfield[0][5],qpos[0][5],qneg[0][5]])])

spec_stats1 = Table([spec_names,total_mems,cl1_mems,cl2_mems],names=('Property','Full catalogue','1.macs0416','2.macs1149'))
spec_stats2 = Table([spec_names,cl3_mems,cl4_mems,cl5_mems,cl6_mems],names=('Property','3.macs0717','4.abell370','5.abell1063','6.abel2744'))
#
#
#
#
##Section 5: make cuts to photometric sample based on defintion for del_z:
##del_z = (z_phot - z_cl) < 0.05; same photometric cut made to 
##specroscopic sub-sample. this is a preliminary measure for determining the 
##photometric sample, final corrections will be made by radial/mass bins to 
##match outlier fractions in spec. sample per van der Burg (2013)
#
## (i) apply cut at |z_clusterphot| < 0.05 to separate cluster from field for 
## targets with photometry only. store in column 'members': 0 = member, 1 = field
## recall: z_clusterphot defined as (z_peak - z_cluster / 1 + z_cluster)
#
### Quiescent sub-sample: segregate into secure member, field
phot_mems = [[0]*6]
phot_memq = [[0]*6]
phot_fields = [[0]*6]
phot_fieldq = [[0]*6]
other = 0       #track phot stars & outliers; should be 89  (spec has 47, total of 136 per 'macros table')
size = len(master_cat)
counter = 0
while counter < size:
    if master_cat[counter]['sub'] ==2:
        if master_cat[counter]['type'] == 0 or master_cat[counter]['type'] ==3: #skip stars and outliers
            other +=1       
        elif master_cat[counter]['type'] ==1:       # SF
            if abs(master_cat[counter]['z_clusterphot']) > 0.05:     
                if master_cat[counter]['z_peak'] >0.5 or master_cat[counter]['z_peak'] <0.3:
                    master_cat[counter]['member'] = 4
                else:
                    master_cat[counter]['member'] = 1               #phot SF field sample
                    if master_cat[counter]['cluster'] == 1:
                        phot_fields[0][0] +=1
                    elif master_cat[counter]['cluster'] == 2:
                        phot_fields[0][1] +=1
                    elif master_cat[counter]['cluster'] == 3:
                        phot_fields[0][2] +=1    
                    elif master_cat[counter]['cluster'] == 4:
                        phot_fields[0][3] +=1    
                    elif master_cat[counter]['cluster'] == 5:
                        phot_fields[0][4] +=1    
                    elif master_cat[counter]['cluster'] == 6:
                        phot_fields[0][5] +=1
            else:
                master_cat[counter]['member'] = 0
                if master_cat[counter]['cluster'] == 1:
                    phot_mems[0][0] +=1
                elif master_cat[counter]['cluster'] == 2:
                    phot_mems[0][1] +=1
                elif master_cat[counter]['cluster'] == 3:
                    phot_mems[0][2] +=1    
                elif master_cat[counter]['cluster'] == 4:
                    phot_mems[0][3] +=1    
                elif master_cat[counter]['cluster'] == 5:
                    phot_mems[0][4] +=1    
                elif master_cat[counter]['cluster'] == 6:
                    phot_mems[0][5] +=1
        elif master_cat[counter]['type'] ==2:       #Q
            if abs(master_cat[counter]['z_clusterphot']) > 0.05:     #phot field
                if master_cat[counter]['z_peak'] >0.5 or master_cat[counter]['z_peak'] <0.3:
                    master_cat[counter]['member'] = 4
                else:
                    master_cat[counter]['member'] = 1       #phot Q field sample
                    if master_cat[counter]['cluster'] == 1:
                        phot_fieldq[0][0] +=1
                    elif master_cat[counter]['cluster'] == 2:
                        phot_fieldq[0][1] +=1
                    elif master_cat[counter]['cluster'] == 3:
                        phot_fieldq[0][2] +=1    
                    elif master_cat[counter]['cluster'] == 4:
                        phot_fieldq[0][3] +=1    
                    elif master_cat[counter]['cluster'] == 5:
                        phot_fieldq[0][4] +=1    
                    elif master_cat[counter]['cluster'] == 6:
                        phot_fieldq[0][5] +=1
            else:
                master_cat[counter]['member'] = 0       #cluster member, phot
                if master_cat[counter]['cluster'] == 1:
                    phot_memq[0][0] +=1
                elif master_cat[counter]['cluster'] == 2:
                    phot_memq[0][1] +=1
                elif master_cat[counter]['cluster'] == 3:
                    phot_memq[0][2] +=1    
                elif master_cat[counter]['cluster'] == 4:
                    phot_memq[0][3] +=1    
                elif master_cat[counter]['cluster'] == 5:
                    phot_memq[0][4] +=1    
                elif master_cat[counter]['cluster'] == 6:
                    phot_memq[0][5] +=1
    counter +=1
#
## make cut to photometric field population to isolate galaxies of similar 
## redshift range as cluster sample. clusters range from 0.308 to 0.545; select 
## field sample in range from 0.3 to 0.55
#aa=0
#bb=0
#counter = 0
#size = len(master_cat)
#while counter < size:
##    if master_cat[counter]['member'] == 1 and master_cat[counter]['sub'] ==2:
 ##       bb+=1
  #      if master_cat[counter]['z_peak'] >0.5 or master_cat[counter]['z_peak'] <0.3:
#            master_cat[counter]['member'] = 4
#            aa+=1
#    counter +=1
## summarize photometric sample info into table
#
phot_names = Column(['secure cluster - SF','secure field - SF','total SF','secure cluster - Q','secure field - Q','total Q','Total'])
total_memp = Column([np.sum(phot_mems),np.sum(phot_fields),np.sum([phot_mems,phot_fields]),np.sum(phot_memq),np.sum(phot_fieldq),np.sum([phot_memq,phot_fieldq]),np.sum([phot_mems,phot_fields,phot_memq,phot_fieldq])])
cl1_memp = Column([phot_mems[0][0],phot_fields[0][0],np.sum([phot_mems[0][0],phot_fields[0][0]]),phot_memq[0][0],phot_fieldq[0][0],np.sum([phot_memq[0][0],phot_fieldq[0][0]]),np.sum([phot_mems[0][0],phot_fields[0][0],phot_memq[0][0],phot_fieldq[0][0]])])
cl2_memp = Column([phot_mems[0][1],phot_fields[0][1],np.sum([phot_mems[0][1],phot_fields[0][1]]),phot_memq[0][1],phot_fieldq[0][1],np.sum([phot_memq[0][1],phot_fieldq[0][1]]),np.sum([phot_mems[0][1],phot_fields[0][1],phot_memq[0][1],phot_fieldq[0][1]])])
cl3_memp = Column([phot_mems[0][2],phot_fields[0][2],np.sum([phot_mems[0][2],phot_fields[0][2]]),phot_memq[0][2],phot_fieldq[0][2],np.sum([phot_memq[0][2],phot_fieldq[0][2]]),np.sum([phot_mems[0][2],phot_fields[0][2],phot_memq[0][2],phot_fieldq[0][2]])])
cl4_memp = Column([phot_mems[0][3],phot_fields[0][3],np.sum([phot_mems[0][3],phot_fields[0][3]]),phot_memq[0][3],phot_fieldq[0][3],np.sum([phot_memq[0][3],phot_fieldq[0][3]]),np.sum([phot_mems[0][3],phot_fields[0][3],phot_memq[0][3],phot_fieldq[0][3]])])
cl5_memp = Column([phot_mems[0][4],phot_fields[0][4],np.sum([phot_mems[0][4],phot_fields[0][4]]),phot_memq[0][4],phot_fieldq[0][4],np.sum([phot_memq[0][4],phot_fieldq[0][4]]),np.sum([phot_mems[0][4],phot_fields[0][4],phot_memq[0][4],phot_fieldq[0][4]])])
cl6_memp = Column([phot_mems[0][5],phot_fields[0][5],np.sum([phot_mems[0][5],phot_fields[0][5]]),phot_memq[0][5],phot_fieldq[0][5],np.sum([phot_memq[0][5],phot_fieldq[0][5]]),np.sum([phot_mems[0][5],phot_fields[0][5],phot_memq[0][5],phot_fieldq[0][5]])])
phot_stats1 = Table([phot_names,total_memp,cl1_memp,cl2_memp],names=('Property','Full catalogue','1.macs0416','2.macs1149'))
phot_stats2 = Table([phot_names,cl3_memp,cl4_memp,cl5_memp,cl6_memp],names=('Property','3.macs0717','4.abell370','5.abell1063','6.abel2744'))
#
## identify BCGs in master_cat
BCG = Table(names=('cluster','sub','type','member','id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','BCG','use_phot','L_u','L_v','L_j','uv','vj','del_z','z_clusterspec','z_clusterphot'))
#
BCG_cl = [[0]*6]
BCG_sub= [[0]*4]
BCG_mem=[[0]*4]
size = len(master_cat)
counter = 0
BCG_counter = 0
while counter < size: 
    if master_cat[counter]['BCG'] == 2: #per .cat documentation, "use_{KS/CH1/CH2}" ==2 identifies BCGs. have confirmed each tracks the same # of BCGs per cluster
        BCG.add_row(master_cat[counter])
        if master_cat[counter]['cluster'] == 1:
            BCG_cl[0][0] +=1
        elif master_cat[counter]['cluster'] == 2:
            BCG_cl[0][1] +=1
        elif master_cat[counter]['cluster'] == 3:
            BCG_cl[0][2] +=1
        elif master_cat[counter]['cluster'] == 4:
            BCG_cl[0][3] +=1
        elif master_cat[counter]['cluster'] == 5:
            BCG_cl[0][4] +=1
        elif master_cat[counter]['cluster'] == 6:
            BCG_cl[0][5] +=1
        if master_cat[counter]['sub'] ==0:
            BCG_sub[0][0]+=1
        elif master_cat[counter]['sub'] ==1:
            BCG_sub[0][1]+=1
        elif master_cat[counter]['sub'] ==2:
            BCG_sub[0][2]+=1
        elif master_cat[counter]['sub'] ==3:
            BCG_sub[0][3]+=1
        if master_cat[counter]['member'] == 0:
            BCG_mem[0][0]+=1
        elif master_cat[counter]['member'] == 1:
            BCG_mem[0][1]+=1
        elif master_cat[counter]['member'] == 2:
            BCG_mem[0][2]+=1
        elif master_cat[counter]['member'] == 3:
            BCG_mem[0][3]+=1
    counter +=1
#
## 
## summarize all member/field data in signle table - "sample_stats"
na = 'n/a'
sample_names = Column(['SF total','SF phot','SF spec','false pos','false neg','fraction false pos/neg','Q total','Q phot','Q spec','false pos','false neg','fraction false pos/neg', 'Total'])
members = Column([np.sum(phot_mems)+np.sum(smem),np.sum(phot_mems),np.sum(smem),np.sum(spos),np.sum(sneg),np.sum(spos)/np.sum(sneg),np.sum(phot_memq)+np.sum(qmem),np.sum(phot_memq),np.sum(qmem),np.sum(qpos),np.sum(qneg),np.sum(qpos)/np.sum(qneg),np.sum(phot_mems)+np.sum(smem)+np.sum(phot_memq)+np.sum(qmem)])
field = Column([np.sum(phot_fields)+np.sum(sfield),np.sum(phot_fields),np.sum(sfield),na,na,na,np.sum(phot_fieldq)+np.sum(qfield),np.sum(phot_fieldq),np.sum(qfield),na,na,na,np.sum(phot_fields)+np.sum(sfield)+np.sum(phot_fieldq)+np.sum(qfield)])
sample_stats = Table([sample_names,members,field],names=('Property','Selected members','"Secure" field'))
## Display all stat summary tables:
#       data_stats
#       macros_type
#       macros1/2
#       spec_stats1/2
#       phot_stats
#       sample_stats
#
###### END
