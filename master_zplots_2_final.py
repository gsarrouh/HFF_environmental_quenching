#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#Created on Fri Jun 26 05:20:20 2020
#
#@author: gsarrouh
#
### WHAT THIS PROGRAM DOES:
### This script reads in all data for the Hubble Frontier Fields images and prepares data for plotting and analysis. Key information is summarized in the tables: 
#
#  
### This program creates plots for the master data catalogue, containing:
### Fig. 1: photometric v specroscopic redshift;   Fig. 2: delta-z (member/field/false pos/false neg segregation);   Fig. 3: UVJ diagram( "colour-colour" plot)
#
#
#
### Section summary:
#
### PROGRAM START
#
### (0) FLAG: preliminary section to set which plots get created
### (1) add PLOT_FLAG_1: create plot comparing photometric v specroscopic redshift (Fig. 1)
### (2) add PLOT_FLAG_2: create "delta-z" plot, scatterplot of member/field/false pos/false neg by type (SF/Q - Fig. 2)
### (3) add PLOT_FLAG_3: create UVJ diagram per VDB et al. 2013 to segregate between type (SF/Q - i.e. a "colour-colour" plot)
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
# this next line is specific to jupyter notebook, and allows for Figure editting in a GUI, instead of inline
###%matplotlib qt
#  
# import modules
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import astropy
from astropy.table import Table
from astropy.table import Column
import time
#
#import master_data_7_final.py
#master_data_7_final.init()
#
#
#
## SECTION (0): set PLOT_FLAGS    # 0=off (i.e. don't make plot), 1=on (i.e. make plot)
### MAY NEED TO EDIT ### diag_flag_1
#
plot_flag_1 = 1           # Fig.1 - z_spec v z_phot plot
plot_flag_2 = 0           # Fig.2 - del_z plot
plot_flag_3 = 0           # Fig.3 - UVJ diagram
#
#
#
## SECTION (1): create plot comparing photometric ('z_peak') v specroscopic redshift ('z_spec') (Fig. 1)
#
if plot_flag_1 == 1:
    plt.close()
    plt.clf()
    #
    ## PLOT 1: make z_spec vs z_phot plot
    #
    zees = plt.figure(num=1)
    zvz = zees.add_subplot(111)
    ## plot spec subsample using [ 'sub'==1 (phot+spec subsample) & 'sub'!=4 (stars) ] as identifier for loop
    mem_scatter = []               # list del_z values for caluclating stats later
    outlier_scatter = []
    SF_scatter = []
    Q_scatter = []
    stars = 0             # count total # of outliers
    count_mem = 0               # count the rest
    for counter in range(len(master_cat)):
        if master_cat['sub'][counter] == 1:    # sub=1 for (spec+phot) subsample
            if master_cat['type'][counter] == 3:    # type=3 for outliers
                outz = zvz.scatter(master_cat['z_spec'][counter],master_cat['z_peak'][counter],c='r', marker='v', linewidths = 0)
                outlier_scatter.append(np.abs(master_cat['del_z'][counter]))
                #print('Outlier: %s'%master_cat['type'][counter])
            elif master_cat['type'][counter] !=0 :
                memz = zvz.scatter(master_cat['z_spec'][counter],master_cat['z_peak'][counter],c='b', marker='^', linewidths = 0)
                count_mem+=1
                if master_cat['type'][counter] == 1:         # type=1 for SF
                    SF_scatter.append(np.abs(master_cat['del_z'][counter]))
                elif master_cat['type'][counter] == 2:       # type=1 for Q
                    Q_scatter.append(np.abs(master_cat['del_z'][counter]))
                mem_scatter.append(np.abs(master_cat['del_z'][counter]))
                #print('Non-outlier: %s'%master_cat['type'][counter])
            else:
                if master_cat['type'][counter] == 0:      #type=0 for stars
                    stars+=1
    #
    ## compute some stats
    SF_scatter = np.array(SF_scatter)
    Q_scatter = np.array(Q_scatter)
    mem_scatter = np.array(mem_scatter)
    outlier_scatter = np.array(outlier_scatter)             # convert to array for math operations
    abs_scatter = np.concatenate([mem_scatter,outlier_scatter])    # aggregate into a single array
    mean_abs_delz = np.array([np.mean(mem_scatter),np.mean(outlier_scatter),np.mean(abs_scatter)])   # compute mean of |del_z|; [SF, Q, total]
    std_abs_delz = np.array([np.std(mem_scatter),np.std(outlier_scatter),np.std(abs_scatter)])   # compute std dev of |del_z|; [SF, Q, total]
    outlier_fraction = len(outlier_scatter) / (len(outlier_scatter)+len(mem_scatter))
    #
    #
    ### MAY NEED TO EDIT ### diag_flag_1
    ##  compare outliers & |del_z| stats identified in master_data_7_final.py with those identified above
    diag_flag_1 = 1             # 0=off (don't display diagnostic); 1=on (display diagnostic table)
    #
    if diag_flag_1 == 1:
        print('Data preparation file finds the following\nOUTLIERS total: %s' % np.sum(outliers),'\nOutlier fraction: %s' % (np.sum(outliers)/np.sum(both)),'\n|del_z| mean: %s'%delz_mean,'\n|del_z| scatter: %s\n'%delz_scatter)
        print('This plotting program found the following\nOUTLIERS total: %s' %len(outlier_scatter),'\nOutlier fraction: %s' %outlier_fraction,'\n|del_z| mean: %s'%mean_abs_delz[2],'\n|del_z| scatter: %s\n'%std_abs_delz[2])
        #print('FULL MEAN: %s'%np.mean(full_scatter),'\nFULL SCATTER: %s'%np.std(full_scatter))
        print('DIFFERENCES\nOutlier total: %s'%(np.sum(outliers)-count_outlier),'\nOutlier fraction: %s'%((np.sum(outliers)/np.sum(both))-outlier_fraction),'\n|del_z| mean: %s'%(delz_mean-mean_abs_delz[2]),'\n|del_z| scatter: %s\n'%(delz_scatter-std_abs_delz[2]))
        print('# of stars: %s'%stars,'\nNon-outlier count: %s'%len(outlier_scatter),'\nNon-outlier count: %s'%count_mem)
        print('SF scatter: %s'%np.std(SF_scatter),'\nQ scatter: %s'%np.std(Q_scatter))

              
              
    ## construct a string to plot the outlier fraction & std dev ('scatter')
    string = 'Mean |$\Delta$z|: %.3f'%mean_abs_delz[2]+'\n$\sigma_{z}$ = %.3f'%std_abs_delz[2]+'\nOutliers: ~%.3f'%outlier_fraction
    ## add text to plot
    plt.text(0.05,0.7,string)
    plt.plot([0,2],[0,2],':k', linewidth=1)
    plt.xlabel("$z_{spec}$")
    plt.xscale('linear')
    plt.xlim(0,1)
    plt.ylabel("$z_{phot}$")
    plt.yscale('linear')
    plt.ylim(0,1)
    plt.title("$z_{spec} vs. z_{phot}$")
    plt.legend((memz,outz),('spec population','outliers'),scatterpoints=1,loc='upper left', frameon=False)
    #plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
    plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on')
    plt.minorticks_on()
    plt.show()
    #
#
#
#
#
#
## SECTION (2): del_z(spec) vs del_z(phot)
#
if plot_flag_2 == 1:
    delzees = plt.figure(num=2)
    delzSF = delzees.add_subplot(121) #SF plot
    #
    SF_membership = ([0]*4) # [(memebrs),(field),(false pos),(false neg)]
    Q_membership = ([0]*4) # [(memebrs),(field),(false pos),(false neg)]
    counter = 0
    size = len(master_cat)
    while counter < size:
        if master_cat[counter]['sub'] == 1 and master_cat[counter]['type']==1:      #only look at SF spec sub-sample, type 1 = SF
            if master_cat[counter]['member'] == 0:        #secure cluster member
                SF_membership[0] +=1
                Smem = delzSF.scatter(master_cat[counter]['z_clusterspec'],master_cat[counter]['z_clusterphot'],c='g', marker='+', linewidths = 0)
            elif master_cat[counter]['member'] == 1:        #secure field
                SF_membership[1] +=1
                Sfield = delzSF.scatter(master_cat[counter]['z_clusterspec'],master_cat[counter]['z_clusterphot'],c='b', marker='+', linewidths = 0)
            elif master_cat[counter]['member'] == 2:        #false pos
                SF_membership[2] +=1
                Spos = delzSF.scatter(master_cat[counter]['z_clusterspec'],master_cat[counter]['z_clusterphot'],c='r', marker='x', linewidths = 0)
            elif master_cat[counter]['member'] == 3:        #false neg
                SF_membership[3] +=1
                Sneg = delzSF.scatter(master_cat[counter]['z_clusterspec'],master_cat[counter]['z_clusterphot'],c='m', marker='x', linewidths = 0)
        counter +=1
    #
    plt.xlabel('$(z_{spec} - z_{cluster})/(1+z_{spec})$')
    plt.xscale('linear')
    plt.xlim(-0.125,0.125)
    plt.tick_params(axis='both', which='both',direction='in',color='k')
    plt.ylabel('$(z_{phot} - z_{cluster})/(1+z_{phot})$')
    delzSF.yaxis.tick_left()
    plt.yscale('linear')
    plt.ylim(-0.15,0.15)
    lin = delzSF.plot([-0.5,1],[0.03,0.03],':k', linewidth=1)  # horizontal cuts
    lin = delzSF.plot([-0.5,1],[-0.03,-0.03],':k', linewidth=1)
    lin = delzSF.plot([-0.01,-0.01],[-0.5,1],':k', linewidth=1)  #vertical cuts
    lin = delzSF.plot([0.01,0.01],[-0.5,1],':k', linewidth=1)  
    plt.legend((Smem,Sfield,Spos,Sneg),('secure member','secure field','false positive','false negative'),scatterpoints=1,loc='upper left', fontsize=8, frameon=False)
    plt.title('Star-Forming')
    plt.minorticks_on()
    plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=False,labelleft=True)
    #
    #
    #
    delzPass = delzees.add_subplot(122) #Passive plot
    counter = 0
    size = len(master_cat)
    while counter < size:
        if master_cat[counter]['sub'] == 1 and master_cat[counter]['type']==2:      #only look at Q spec sub-sample, type 2 = Q
            if master_cat[counter]['member'] == 0:        #secure cluster member
                Q_membership[0] +=1
                Qmem = delzPass.scatter(master_cat[counter]['z_clusterspec'],master_cat[counter]['z_clusterphot'],c='g', marker='+', linewidths = 0)
            elif master_cat[counter]['member'] == 1:        #secure field
                Q_membership[1] +=1
                Qfield = delzPass.scatter(master_cat[counter]['z_clusterspec'],master_cat[counter]['z_clusterphot'],c='b', marker='+', linewidths = 0)
            elif master_cat[counter]['member'] == 2:        #false pos
                Q_membership[2] +=1
                Qpos = delzPass.scatter(master_cat[counter]['z_clusterspec'],master_cat[counter]['z_clusterphot'],c='r', marker='x', linewidths = 0)
            elif master_cat[counter]['member'] == 3:        #false neg
                Q_membership[3] +=1
                Qneg = delzPass.scatter(master_cat[counter]['z_clusterspec'],master_cat[counter]['z_clusterphot'],c='m', marker='x', linewidths = 0)
        counter +=1
    plt.xlabel('$(z_{spec} - z_{cluster})/(1+z_{spec})$')
    plt.xscale('linear')
    plt.xlim(-0.125,0.125)
    plt.tick_params(axis='both', which='both',direction='in',color='k',top=True)
    plt.ylabel('$(z_{phot} - z_{cluster})/(1+z_{phot})$')
    #delzPass.yaxis.tick_right()
    delzPass.yaxis.set_label_position("right")
    plt.yscale('linear')
    plt.ylim(-0.15,0.15)
    lin = delzPass.plot([-0.5,1],[0.07,0.07],':k', linewidth=1)  # horizontal cuts
    lin = delzPass.plot([-0.5,1],[-0.07,-0.07],':k', linewidth=1)
    lin = delzPass.plot([-0.01,-0.01],[-0.5,1],':k', linewidth=1)  #vertical cuts
    lin = delzPass.plot([0.01,0.01],[-0.5,1],':k', linewidth=1)  
    #plt.legend((SFz,Passz),('del_z < 0.1','del_z > 0.1'),scatterpoints=1,loc='lower left', fontsize=10)
    plt.title('Passive')
    plt.minorticks_on()
    plt.tick_params(axis='y', which='both', direction='in',color='k',top=True,right=True,labelright=True,labelleft=False)
    #
    plt.show()
    #
#
#
#
#
# SECTION 3:
#
## colour-colour plots; modify to add in photometric sub-sample
#
if plot_flag_3 ==1:
    check = 0
    aa = 0
    bb = 0          #counting variables to check plotting code
    cc = 0
    colcol = plt.figure(num=3)
    cvc = colcol.add_subplot(111)
    counter = 0
    size = len(master_cat)
    while counter < size:
        if master_cat[counter]['member'] ==0:
            if master_cat[counter]['type'] ==1:
                aa +=1      #count for SF
                SF = cvc.scatter(master_cat[counter]['vj'],master_cat[counter]['uv'], c='b', marker='*', linewidths=0)
            elif master_cat[counter]['type']==2:
                bb +=1    #count for Q
                Q = cvc.scatter(master_cat[counter]['vj'],master_cat[counter]['uv'], c='r', marker='.', linewidths=0)
            else:
                cc +=1
        else:
            check +=1
        counter +=1  
    bounds = cvc.plot([0,0.8],[1.3,1.3],'-k',[0.8,1.6],[1.3,2],'-k',[1.6,1.6],[2,2.5],'-k', linewidth=1) #overlay boundary cutoff for SF/Passive
    plt.xscale('linear')
    plt.xlabel('$(V-J)_{rest}$')
    plt.xlim(0,2)
    plt.yscale('linear')
    plt.ylabel('$(U-V)_{rest}$')
    plt.ylim(0,2.5)
    #plt.title('V-J vs U-V')
    plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on')
    plt.minorticks_on()
    plt.text(0.25,2,'Passive',fontsize=10)
    plt.text(1.5,0.75,'Star-Forming',fontsize=10)
    #plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')
    #
    plt.show()
#
#
#
#
#
## TIME_FLAG END
#
if time_flag == 1:
    print('Program "master_zplots_2_final.py" took: %s seconds to run.\n\n' % (time.time() - start_time))
#
#
#
#
#                        
###### PROGRAM END ######