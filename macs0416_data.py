#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 15:19:59 2017

@author: gsarrouh
"""
#
###This script reads in data related to MACS0416 and prepares data for plotting
###and analysis. Key information is summarized in the tables: "macros",
#
###v2 has updated section descriptions to improve readability; and includes cuts 
###made to del_z for spec sample
#
###v4 removes the magnitude freq distribution which was written for comparison 
###w/ Jauzac et al 2014; and selects stars using the "star_flag" data; removes 
###class_star columns to save on space
#
# Read in macs0416 data from WORKING DIRECTORY: NSERC17/HFFtoAdam/working_data
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import astropy
from astropy.table import Table
from astropy.table import Column
#
##Section 1: import all data from HFF team, convert flux to luminosity & gather full 
##catalogue into single table "macs0416_full"; separate objects for which there is no 
##redshift data (both photo & spec) from main working catalogue "macs0416"
#
#create table from redshift ".zout" file
global z_macs0416
z_macs0416 = Table.read('abell370clu_v3.5.zout',format='ascii')

#create table from redshift ".fout" file
global f_macs0416
f_macs0416 = Table.read('abell370clu_v3.5.fout',format='ascii')

#read in the whole bloody catalogue to sift for stars
global cat_macs0416
cat_macs0416 = Table.read('hffds_abell370clu_v3.5.cat',format='ascii')

#creat table for colours
dat_u = Table.read('abell370clu_v3.5.153.rf',format='ascii')
dat_v = Table.read('abell370clu_v3.5.155.rf',format='ascii')
dat_j = Table.read('abell370clu_v3.5.161.rf',format='ascii')

#aggregate into a single table
global macs0416
macs0416 = Table([z_macs0416['id'],z_macs0416['z_peak'],z_macs0416['z_spec'],dat_u['L153'],dat_v['L155'],dat_j['L161'],dat_u['DM'],f_macs0416['lmass'],f_macs0416['lsfr'],f_macs0416['lssfr'],cat_macs0416['flux_radius'],cat_macs0416['star_flag']], names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag'))
#macs0416['z_peak'] = macs0416['z_peak'] + 0.05  #adjust for systemic bias in z_peak

# sift for targets without z_phot (there are no targets w/ z_spec w/o z_phot, so this gets all targets w/ no relaiable data)
no_data = Table(names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag'))
nodata = 0
size = len(macs0416)-1
counter = 0
while counter < size: 
    if macs0416[counter]['z_peak'] < 0: #identify entries with no z_spec estimate
        no_data.add_row(macs0416[counter])
        macs0416.remove_row(counter) 
        nodata +=1
        size -=1
    else:
        counter +=1
#
#convert flux to luminosity, get color indices U-V, V-J, add to table
L_u = -2.5*np.log10(macs0416['u']) + 25
L_v = -2.5*np.log10(macs0416['v']) + 25
L_j = -2.5*np.log10(macs0416['j']) + 25

L_u = Column(L_u, name='L_u')
L_v = Column(L_v, name='L_v')
L_j = Column(L_j, name='L_j')
#
uv = Column(L_u - L_v, name='uv')
vj = Column(L_v - L_j, name='vj')
#
macs0416.add_columns([uv,vj,L_u,L_v,L_j]) #append U-V, V-J to table
global macs0416_full
macs0416_full = Table(macs0416)  #create a copy table with all data in catalogue, prior to segregation
#
##Section 2: begin filtering catalogue for sample membership:
##(i) stars; (ii) z_phot only in "onlypeak"; (iii) both z_spec & z_phot (i.e. z_peak) in "macs0416"
##
# sift for stars in "class_star_F814W" from catalogue; cutoff for star is >0.9
global stars
stars = Table(names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','uv','vj','L_u','L_v','L_j'))
size = len(macs0416)
counter = 0
while counter < size: 
    if macs0416[counter]['star_flag'] == 1: #identify objects with a "star flag" = 1 for stars; gal=0, mag>25 = 2
        stars.add_row(macs0416[counter])
        macs0416.remove_row(counter)
        size -=1
    else:
        counter +=1
#
#sift for tagets w/o z_spec, store targets w/ only z_peak in separate table
global onlypeak
onlypeak = Table(names=('id','z_peak','z_spec','u','v','j','DM','lmass','lsfr','lssfr','flux_radius','star_flag','uv','vj','L_u','L_v','L_j'))
size = len(macs0416)
counter = 0
while counter < size: 
    if macs0416[counter]['z_spec'] <0: #identify entries with no z_spec estimate
        onlypeak.add_row(macs0416[counter])
        macs0416.remove_row(counter) 
        size -=1
    else:
        counter +=1
#
#NOTE:   "onlypeak" now contains the sub-sample of targets which have photometry ONLY
#NOTE:   "macs0416" now contains the SPEC sub-sample of targets with both spectroscopy & photometry
#
#
##   **Section 3: separate SF/Q galaxies for both phot & spec sub-samples:
##(i) det delta_z & del_z(cluster) for both phot & spec; (ii) separate SF/Q 
##based on van der Burg (2013) colour criteria; (iii) separate outliers from 
##both phot & spec sub-sample; 
#
#(i) calculate delta_z for targets w/ both z_spec & z_phot
del_z = abs(macs0416['z_peak'] - macs0416['z_spec']) / (1 + macs0416['z_spec'])
#
z_cluster = 0.375   #https://archive.stsci.edu/prepds/frontier/
#
z_clusterphot = (macs0416['z_peak'] - z_cluster) / (1 + macs0416['z_peak'])
#z_clusterphot = z_clusterphot + 0.05 #correct for systemic error
z_clusterspec = (macs0416['z_spec'] - z_cluster) / (1 + macs0416['z_spec'])
#
#make tables for SF/passive plot data; for both spec & phot sub-samples
global SF_spec
global Q_spec
global SF_phot
global Q_phot
global SF_spec_out
global Q_spec_out
#
SF_spec = Table([macs0416['id'],macs0416['z_spec'],macs0416['z_peak'],del_z, z_clusterphot, z_clusterspec, macs0416['uv'], macs0416['vj'], macs0416['lmass'],macs0416['lsfr'],macs0416['lssfr'],macs0416['L_u'],macs0416['L_v'],macs0416['L_j'],macs0416['star_flag']],names=('id','z_spec','z_peak','del_z','z_clusterphot','z_clusterspec','uv','vj','lmass','lsfr','lssfr','L_u','L_v','L_j','star_flag'))
Q_spec = Table(names=('id','z_spec','z_peak','del_z','z_clusterphot','z_clusterspec','uv','vj','lmass','lsfr','lssfr','L_u','L_v','L_j','star_flag'))
SF_phot = Table([onlypeak['id'],onlypeak['z_peak'], onlypeak['uv'], onlypeak['vj'], onlypeak['lmass'],onlypeak['lsfr'],onlypeak['lssfr'],onlypeak['L_u'],onlypeak['L_v'],onlypeak['L_j'],onlypeak['star_flag']],names=('id','z_peak','uv','vj','lmass','lsfr','lssfr','L_u','L_v','L_j','star_flag'))
Q_phot = Table(names=('id','z_peak','uv','vj','lmass','lsfr','lssfr','L_u','L_v','L_j','star_flag'))
SF_spec_out = Table(names=('id','z_spec','z_peak','del_z','z_clusterphot','z_clusterspec','uv','vj','lmass','lsfr','lssfr','L_u','L_v','L_j','star_flag'))
Q_spec_out = Table(names=('id','z_spec','z_peak','del_z','z_clusterphot','z_clusterspec','uv','vj','lmass','lsfr','lssfr','L_u','L_v','L_j','star_flag'))
#
#(ii) separate SF/passive based on criteria from vanderBurg 2013
#spec sub-sample
size = len(SF_spec)
counter = 0
while counter < size: 
    if SF_spec[counter]['uv'] > 1.3 and SF_spec[counter]['vj'] < 1.6 and SF_spec[counter]['uv'] > ((0.88*SF_spec[counter]['vj']) + 0.6): #identify objects 
        Q_spec.add_row(SF_spec[counter])
        SF_spec.remove_row(counter)
        size -=1
    else:
        counter +=1
#
#phot sub-sample
size = len(SF_phot)-1
counter = 0
while counter < size: 
    if SF_phot[counter]['uv'] > 1.3 and SF_phot[counter]['vj'] < 1.6 and SF_phot[counter]['uv'] > ((0.88*SF_phot[counter]['vj']) + 0.6): #identify objects 
        Q_phot.add_row(SF_phot[counter])
        SF_phot.remove_row(counter)
        size -=1
    else:
        counter +=1
#
#(iii) separate SF outliers for |del_z| > 0.15 (spec sub-sample)
size = len(SF_spec)
counter = 0
while counter < size: 
    if SF_spec[counter]['del_z'] > 0.15: #identify objects w/ del_z > 0.15
        SF_spec_out.add_row(SF_spec[counter])
        SF_spec.remove_row(counter)
        size -=1
    else:
        counter +=1
#
#separate Passive outliers for |del_z| > 0.15 (spec sub-sample)
size = len(Q_spec)
counter = 0
while counter < size: 
    if Q_spec[counter]['del_z'] > 0.15: #identify objects w/ del_z > 0.15
        Q_spec_out.add_row(Q_spec[counter])
        Q_spec.remove_row(counter)
        size -=1
    else:
        counter +=1
# 
####create summary table of macro stats
out_frac_SF = len(SF_spec_out)/(len(SF_spec) + len(SF_spec_out))
out_frac_Q = len(Q_spec_out)/(len(Q_spec) + len(Q_spec_out))
out_frac_total = (len(SF_spec_out)+len(Q_spec_out))/(len(SF_spec) + len(SF_spec_out) + len(Q_spec) + len(Q_spec_out))
macro_names = Column(['full catalogue','no data','stars','total spec','total photo','spec members','spec SF','spec Q','spec outliers','total outlier fraction','spec SF outliers','SF outlier fraction','spec Q outliers','Q outliers fraction','total_spec','phot SF','phot Q','total_phot'])
macro = Column([len(cat_macs0416),nodata,len(stars),len(macs0416),len(onlypeak),len(SF_spec)+len(Q_spec),len(SF_spec),len(Q_spec),len(SF_spec_out)+len(Q_spec_out),out_frac_total,len(SF_spec_out),out_frac_SF,len(Q_spec_out),out_frac_Q,len(macs0416),len(SF_phot), len(Q_phot), len(onlypeak)])
global macros0416
macros0416 = Table([macro_names,macro],names=('Property','Value'))
#
##   **Section 4: apply delta_z cuts to determine membership of z_phot sample based on cut
##criteria discussed w/ AM: |del_z phot| < 0.05 & |del_z phot| < 0.01, in the process isolate 
##secure field, false pos/neg items from         
#separate  secure field/member & false pos/neg using criteria: 
#mem: abs(z_clusterspec) < 0.01 & abs(z_cluster phot) < 0.05; keep in SF_spec_ table
#field: abs(z_clusterspec) > 0.01 & abs(z_cluster phot) > 0.05; remove to SF_field
#false pos: abs(z_clusterspec) > 0.01 & abs(z_cluster phot) < 0.05; remove to SF_fpos
#false neg: abs(z_clusterspec) < 0.01 & abs(z_cluster phot) > 0.05; remove to SF_fneg
#Note: this is only done for spec sample; results used to correct photo sample for completeness at end of analysis
#Note: SF_spec will be listing of spectroscopically confirmed members (SF), Q_spec is spec. confirmed quiescient galaxies

global SF_field     ##create new tables to sort spec sample
global SF_fpos
global SF_fneg
global Q_field 
global Q_fpos
global Q_fneg
SF_field = Table(names=('id','z_spec','z_peak','del_z','z_clusterphot','z_clusterspec','uv','vj','lmass','lsfr','lssfr','L_u','L_v','L_j','star_flag'))
SF_fpos = Table(names=('id','z_spec','z_peak','del_z','z_clusterphot','z_clusterspec','uv','vj','lmass','lsfr','lssfr','L_u','L_v','L_j','star_flag'))
SF_fneg = Table(names=('id','z_spec','z_peak','del_z','z_clusterphot','z_clusterspec','uv','vj','lmass','lsfr','lssfr','L_u','L_v','L_j','star_flag'))
Q_field = Table(names=('id','z_spec','z_peak','del_z','z_clusterphot','z_clusterspec','uv','vj','lmass','lsfr','lssfr','L_u','L_v','L_j','star_flag'))
Q_fpos = Table(names=('id','z_spec','z_peak','del_z','z_clusterphot','z_clusterspec','uv','vj','lmass','lsfr','lssfr','L_u','L_v','L_j','star_flag'))
Q_fneg = Table(names=('id','z_spec','z_peak','del_z','z_clusterphot','z_clusterspec','uv','vj','lmass','lsfr','lssfr','L_u','L_v','L_j','star_flag'))
  
###SF sub-sample: SF_spec into secure member, field, false pos/neg
size = len(SF_spec)
counter = 0
while counter < size:
    if abs(SF_spec[counter]['z_clusterspec']) > 0.01 and abs(SF_spec[counter]['z_clusterphot']) > 0.05:  # secure field
        SF_field.add_row(SF_spec[counter])
        SF_spec.remove_row(counter)
        size -=1
    elif abs(SF_spec[counter]['z_clusterspec']) > 0.01 and abs(SF_spec[counter]['z_clusterphot']) < 0.05:  #  false positive
        SF_fpos.add_row(SF_spec[counter])
        SF_spec.remove_row(counter)
        size -=1
    elif abs(SF_spec[counter]['z_clusterspec']) < 0.01 and abs(SF_spec[counter]['z_clusterphot']) > 0.05:  #  false negative
        SF_fneg.add_row(SF_spec[counter])
        SF_spec.remove_row(counter)
        size -=1
    else:
        counter +=1
#
###Passive sub-sample: Q_spec into secure member, field, false pos/neg
size = len(Q_spec)
counter = 0
while counter < size:
    if abs(Q_spec[counter]['z_clusterspec']) > 0.01 and abs(Q_spec[counter]['z_clusterphot']) > 0.05:  # secure field
        Q_field.add_row(Q_spec[counter])
        Q_spec.remove_row(counter)
        size -=1
    elif abs(Q_spec[counter]['z_clusterspec']) > 0.01 and abs(Q_spec[counter]['z_clusterphot']) < 0.05:  #  false positive
        Q_fpos.add_row(Q_spec[counter])
        Q_spec.remove_row(counter)
        size -=1
    elif abs(Q_spec[counter]['z_clusterspec']) < 0.01 and abs(Q_spec[counter]['z_clusterphot']) > 0.05:  #  false negative
        Q_fneg.add_row(Q_spec[counter])
        Q_spec.remove_row(counter)
        size -=1
    else:
        counter +=1
#
##summarize sample stats in table
#SF_spec_mean_delz = np.mean(SF_spec['del_z'])
#Q_spec_mean_delz = np.mean(Q_spec['del_z'])
#SF_spec_mean_z = np.mean(SF_spec['z_spec'])
#Q_spec_mean_z = np.mean(Q_spec['z_spec'])
#spec_mean_z = np.mean(np.concatenate((SF_spec['z_spec'],Q_spec['z_spec'])))
#spec_scatter_z = np.std(np.concatenate((SF_spec['z_spec'],Q_spec['z_spec'])))
#spec_mean_delz = np.mean(np.concatenate((SF_spec['del_z'],Q_spec['del_z'])))
#spec_scatter_delz = np.std(np.concatenate((SF_spec['del_z'],Q_spec['del_z'])))
#spec_names = Column(['Total spec','total SF','total Q','Outliers','SF members','SF field','SF false pos','SF false neg','Q mem','Q field','Q false pos','Q false neg','cluster z','SF members mean z_spec','total member mean z_spec','total scatter','delta z mean - all members','SF mean delta z','delta z scatter - all members']) #'Q members mean z_spec' directly after 'SF mean delta z'; 'Q mean delta z' after 'SF mean delta z'
#spec_stat = Column([len(macs0416),len(SF_spec)+len(SF_field)+len(SF_fpos)+len(SF_fneg),len(Q_spec)+len(Q_field)+len(Q_fpos)+len(Q_fneg),len(SF_spec_out)+len(Q_spec_out),len(SF_spec),len(SF_field),len(SF_fpos),len(SF_fneg),len(Q_spec),len(Q_field),len(Q_fpos),len(Q_fneg),z_cluster,SF_spec_mean_z,spec_mean_z,spec_scatter_z,spec_mean_delz,SF_spec_mean_delz,spec_scatter_delz]) #Q_spec_mean_z,  Q_spec_mean_delz,
#global spec_stats0416
#spec_stats0416 = Table([spec_names,spec_stat],names=('Property','Value'))
#   
#
#
##Section 5: make cuts to photometric sample based on defintion for del_z:
##del_z = (z_phot - z_cl) / (1 + z_cl) < 0.05; same photometric cut made to 
##specroscopic sub-sample. this is a preliminary measure for determining the 
##photometric sample, final corrections will be made by radial/mass bins to 
##match outlier fractions in spec. sample per van der Burg (2013)
#
#(i) calculate del_z for photo sample wrt to cluster
del_z_SF_phot = (SF_phot['z_peak'] - z_cluster) / (1 + z_cluster)
del_z_SF_phot = Column(del_z_SF_phot, name='del_z_phot')
del_z_Q_phot = (Q_phot['z_peak'] - z_cluster) / (1 + z_cluster)
del_z_Q_phot = Column(del_z_Q_phot, name='del_z_phot')
SF_phot.add_column(del_z_SF_phot)
Q_phot.add_column(del_z_Q_phot)
#(ii) apply cut at |del_z| < 0.05 to separate cluster from field
global SF_phot_out
global Q_phot_out
SF_phot_out = Table(names=('id','z_peak','uv','vj','lmass','lsfr','lssfr','L_u','L_v','L_j','star_flag','del_z_phot'))
Q_phot_out = Table(names=('id','z_peak','uv','vj','lmass','lsfr','lssfr','L_u','L_v','L_j','star_flag','del_z_phot'))
##############   good up til here   ###############
#SF sub-sample
size = len(SF_phot)
counter = 0
while counter < size:
    if abs(SF_phot[counter]['del_z_phot']) > 0.05:
        SF_phot_out.add_row(SF_phot[counter])
        SF_phot.remove_row(counter)
        size -=1
    else:
        counter +=1
#Passive sub-sample
size = len(Q_phot)
counter = 0
while counter < size:
    if abs(Q_phot[counter]['del_z_phot']) > 0.05:
        Q_phot_out.add_row(Q_phot[counter])
        Q_phot.remove_row(counter)
        size -=1
    else:
        counter +=1
#
#summarize macro info into table
#SF_phot_mean_delz = np.mean(SF_phot['del_z_phot'])
#Q_phot_mean_delz = np.mean(Q_phot['del_z_phot'])
#SF_phot_mean_z = np.mean(SF_phot['z_peak'])
#Q_phot_mean_z = np.mean(Q_phot['z_peak'])
#phot_mean_z = np.mean(np.concatenate((SF_phot['z_peak'],Q_phot['z_peak'])))
#phot_scatter_z = np.std(np.concatenate((SF_phot['z_peak'],Q_phot['z_peak'])))
#phot_mean_delz = np.mean(np.concatenate((SF_phot['del_z_phot'],Q_phot['del_z_phot'])))
#phot_scatter_delz = np.std(np.concatenate((SF_phot['del_z_phot'],Q_phot['del_z_phot'])))
#phot_mean_out = np.mean(np.concatenate((SF_phot_out['z_peak'],Q_phot_out['z_peak'])))
#phot_names = Column(['Total phot','SF members','Q mem','Outliers','cluster z','SF mean z_phot','Q mean z_phot','total member mean z_phot','total scatter','outliers mean z_phot','delta z mean - all members','SF mean delta z','Q mean delta z','delta z scatter - all members'])
#phot_stat = Column([len(SF_phot)+len(SF_phot_out)+len(Q_phot)+len(Q_phot_out),len(SF_phot),len(Q_phot),len(SF_phot_out)+len(Q_phot_out),z_cluster,SF_phot_mean_z,Q_phot_mean_z,phot_mean_z,phot_scatter_z,phot_mean_out,phot_mean_delz,SF_phot_mean_delz,Q_phot_mean_delz,phot_scatter_delz])
#global phot_stats0416
#phot_stats0416 = Table([phot_names,phot_stat],names=('Property','Value'))
###

#calculate mean, scatter of sample
#mem_mean_spec_delz = 
#mem_std_spec_delz = np.std(np.concatenate((SF_spec['del_z'],Q_spec['del_z'])))
#tot_mean_spec = np.mean(macs0416['z_spec'])
#tot_std_spec = np.std(macs0416['z_spec'])
#mem_mean_spec = np.mean(np.concatenate((SF_spec['z_spec'],Q_spec['z_spec'])))
#mem_std_spec = np.std(np.concatenate((SF_spec['z_spec'],Q_spec['z_spec'])))
#mean_phot = np.mean(onlypeak['z_peak'])
#med_phot = np.mean(onlypeak['z_peak']) 
#std_phot = np.std(onlypeak['z_peak'])
#stat_names = Column(['total catalogue','stars','total spec sub-sample (mem + out)','total photo sub-sample','no data','spec members','spec SF','spec Q','spec outliers','spec SF outliers','spec Q outliers','spec total','total z_spec mean','total z_spec std dev','mem del_z mean ','mem del_z std dev', 'mem z_spec mean','mem z_spec std dev','phot SF', 'phot Q', 'phot total', 'z_phot mean', 'z_phot median', 'z_phot std dev'])
#stat = Column([len(cat_macs0416),len(stars),len(macs0416),len(onlypeak),nodata,len(SF_spec)+len(Q_spec),len(SF_spec),len(Q_spec),len(SF_spec_out)+len(Q_spec_out),len(SF_spec_out),len(Q_spec_out),len(macs0416), tot_mean_spec, tot_std_spec, mem_mean_spec_delz, mem_std_spec_delz, mem_mean_spec, mem_std_spec, len(SF_phot), len(Q_phot), len(onlypeak), mean_phot, med_phot, std_phot])
#global stats
#stats = Table([stat_names,stat],names=('Property','Value'))

#stats.more()
###### END
