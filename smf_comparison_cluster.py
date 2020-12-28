# Created on Sun Dec 18 0842:24 2020
# Author: Ghassan T Sarrouh
#
########################### smf_comparison_cluster.py ###########################
#
### WHAT THIS PROGRAM DOES:
### creates figures comparing the SMF result produced in "master_smf_9_final.py"; comparison data values have been manually entered from the reference provided
#
#
#
### PROGRAM START
#
### (0)     import modules, define functions, FLAGS;
### (1)     import DATA POINTS from literature
### (1.1)
### (1.2)   VISUALIZE
### (2)     import best-fitting SCHECHTER PARAMETERS
### (2.1)
### (2.2)
### (2.3)   VISUALIZE
### (3)     paneled plot comparing INDIVIDUAL STUDIES to our results
#
### PROGRAM END
#
#
## SECTION (0): import MODULES, set FLAGS, define FUNCTIONS
#
## Import modules
import pandas as pd
import numpy as np
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
from astropy.table import Table
from astropy.table import Column
#
#
## FLAGS
vdb2018_plot_flag = 1       # 1== R < 2*R_500;  2== R < 0.5*R_500

plot_flag_1 = 0             # unassigned
plot_flag_2 = 1
plot_flag_3 = 0             # block out old code

#
#
## FUNCTIONS
#

#
#
#
#
## SECTION (1):Fig. 1: compare various SMF DATA POINTS with my SCHECHTER FUNCTIONS
#
## SECTION (1.1): data entry (manual)
#
#
#
## from van der Burg 2018:
#
# x vector for plotting
x_vdb2018 = np.arange(9.55,12.01,0.01)
#
if vdb2018_plot_flag == 1:
    # R < 2R_500
    # SF population
    SF_M_star_vdb2018, SF_phi1_vdb2018, SF_alpha1_vdb2018 = 10.50,SFphi_mcmc,-1.02
    SF_M_star_upper_err_vdb2018, SF_phi1_upper_err_vdb2018, SF_alpha1_upper_err_vdb2018 = 0.04,SFphi_mcmc_upper_err,0.06
    SF_M_star_lower_err_vdb2018, SF_phi1_lower_err_vdb2018, SF_alpha1_lower_err_vdb2018 = 0.04,SFphi_mcmc_lower_err,0.06
    # Q populaQion
    Q_M_star_vdb2018, Q_phi1_vdb2018, Q_alpha1_vdb2018 = 10.81,Qphi_mcmc,-0.83
    Q_M_star_upper_err_vdb2018, Q_phi1_upper_err_vdb2018, Q_alpha1_upper_err_vdb2018 = 0.01,Qphi_mcmc_upper_err,0.03
    Q_M_star_lower_err_vdb2018, Q_phi1_lower_err_vdb2018, Q_alpha1_lower_err_vdb2018 = 0.02,Qphi_mcmc_lower_err,0.02
    # T population
    T_M_star_vdb2018, T_phi1_vdb2018, T_alpha1_vdb2018 = 10.81,Tphi_mcmc,-0.91
    T_M_star_upper_err_vdb2018, T_phi1_upper_err_vdb2018, T_alpha1_upper_err_vdb2018 = 0.02,Tphi_mcmc_upper_err,0.02
    T_M_star_lower_err_vdb2018, T_phi1_lower_err_vdb2018, T_alpha1_lower_err_vdb2018 = 0.02,Tphi_mcmc_lower_err,0.02
#
elif vdb2018_plot_flag == 2:
    # R < 0.5R_500
    # SF population
    SF_M_star_vdb2018, SF_phi1_vdb2018, SF_alpha1_vdb2018 = 10.69,SFphi_mcmc,-1.11
    SF_M_star_upper_err_vdb2018, SF_phi1_upper_err_vdb2018, SF_alpha1_upper_err_vdb2018 = 0.13,SFphi_mcmc_upper_err,0.15
    SF_M_star_lower_err_vdb2018, SF_phi1_lower_err_vdb2018, SF_alpha1_lower_err_vdb2018 = 0.11,SFphi_mcmc_lower_err,0.14
    # Q populaQion
    Q_M_star_vdb2018, Q_phi1_vdb2018, Q_alpha1_vdb2018 = 10.82,Qphi_mcmc,-0.78
    Q_M_star_upper_err_vdb2018, Q_phi1_upper_err_vdb2018, Q_alpha1_upper_err_vdb2018 = 0.03,Qphi_mcmc_upper_err,0.04
    Q_M_star_lower_err_vdb2018, Q_phi1_lower_err_vdb2018, Q_alpha1_lower_err_vdb2018 = 0.03,Qphi_mcmc_lower_err,0.03
    # T population
    T_M_star_vdb2018, T_phi1_vdb2018, T_alpha1_vdb2018 = 10.81,Tphi_mcmc,-0.81
    T_M_star_upper_err_vdb2018, T_phi1_upper_err_vdb2018, T_alpha1_upper_err_vdb2018 = 0.02,Tphi_mcmc_upper_err,0.03
    T_M_star_lower_err_vdb2018, T_phi1_lower_err_vdb2018, T_alpha1_lower_err_vdb2018 = 0.02,Tphi_mcmc_lower_err,0.02
    #
#
## SECTION (2.2): build SCHECHTER MODELS
#
## construct upper/lower uncertainty limits for MY FIT (Sarrouh & Muzzin 2021)
# SF
SF_model_upper =  np.log(10) * np.exp(-10**(x_plot-(SFM_star_mcmc+SFM_star_mcmc_upper_err))) * ((SFphi_mcmc+SFphi_mcmc_upper_err) * (10**(x_plot-(SFM_star_mcmc+SFM_star_mcmc_upper_err)))**(1+(SFalpha_mcmc-SFalpha_mcmc_upper_err)))
SF_model_lower =  np.log(10) * np.exp(-10**(x_plot-(SFM_star_mcmc-SFM_star_mcmc_lower_err))) * ((SFphi_mcmc-SFphi_mcmc_lower_err) * (10**(x_plot-(SFM_star_mcmc-SFM_star_mcmc_lower_err)))**(1+(SFalpha_mcmc+SFalpha_mcmc_lower_err)))
# Q
Q_model_upper =  np.log(10) * np.exp(-10**(x_plot-(QM_star_mcmc+QM_star_mcmc_upper_err))) * ((Qphi_mcmc+Qphi_mcmc_upper_err) * (10**(x_plot-(QM_star_mcmc+QM_star_mcmc_upper_err)))**(1+(Qalpha_mcmc-Qalpha_mcmc_upper_err)))
Q_model_lower =  np.log(10) * np.exp(-10**(x_plot-(QM_star_mcmc-QM_star_mcmc_lower_err))) * ((Qphi_mcmc-Qphi_mcmc_lower_err) * (10**(x_plot-(QM_star_mcmc-QM_star_mcmc_lower_err)))**(1+(Qalpha_mcmc+Qalpha_mcmc_lower_err)))
# TOTAL
total_model_upper =  np.log(10) * np.exp(-10**(x_plot-(TM_star_mcmc+TM_star_mcmc_upper_err))) * ((Tphi_mcmc+Tphi_mcmc_upper_err) * (10**(x_plot-(TM_star_mcmc+TM_star_mcmc_upper_err)))**(1+(Talpha_mcmc-Talpha_mcmc_upper_err)))
total_model_lower =  np.log(10) * np.exp(-10**(x_plot-(TM_star_mcmc-TM_star_mcmc_lower_err))) * ((Tphi_mcmc-Tphi_mcmc_lower_err) * (10**(x_plot-(TM_star_mcmc-TM_star_mcmc_lower_err)))**(1+(Talpha_mcmc+Talpha_mcmc_lower_err)))
#
#
# build MODELS FROM LITERATURE for SMF, SMF_upper_error_bound & SMF_lower_error_bound
#
#
# van der burg+2018
#SMFs
# SF_model_vdb2018 = np.log(10) * SF_phi1_vdb2018 * ( 10**((x_vdb2018-SF_M_star_vdb2018) * (1+SF_alpha1_vdb2018)) ) * np.exp(-10**(x_vdb2018-SF_M_star_vdb2018))
SF_model_vdb2018 = np.log(10) * SF_phi1_vdb2018 * ( 10**((x_vdb2018-SF_M_star_vdb2018) * (1+SF_alpha1_vdb2018)) ) * np.exp(-10**(x_vdb2018-SF_M_star_vdb2018))
Q_model_vdb2018 = np.log(10) * Q_phi1_vdb2018 * (10**((x_vdb2018-Q_M_star_vdb2018)*(1+Q_alpha1_vdb2018))) * np.exp(-10**(x_vdb2018-Q_M_star_vdb2018))
T_model_vdb2018 = np.log(10) * T_phi1_vdb2018 * (10**((x_vdb2018-T_M_star_vdb2018)*(1+T_alpha1_vdb2018))) * np.exp(-10**(x_vdb2018-T_M_star_vdb2018))
# upper_error_bound
SF_model_upper_err_vdb2018 = np.log(10) * (SF_phi1_vdb2018+SF_phi1_upper_err_vdb2018) * (10**((x_vdb2018-(SF_M_star_vdb2018+SF_M_star_upper_err_vdb2018))*(1+(SF_alpha1_vdb2018-SF_alpha1_upper_err_vdb2018)))) * np.exp(-10**(x_vdb2018-(SF_M_star_vdb2018+SF_M_star_upper_err_vdb2018)))
Q_model_upper_err_vdb2018 = np.log(10) * (Q_phi1_vdb2018+Q_phi1_upper_err_vdb2018) * (10**((x_vdb2018-(Q_M_star_vdb2018+Q_M_star_upper_err_vdb2018))*(1+(Q_alpha1_vdb2018-Q_alpha1_upper_err_vdb2018)))) * np.exp(-10**(x_vdb2018-(Q_M_star_vdb2018+Q_M_star_upper_err_vdb2018)))
T_model_upper_err_vdb2018 = np.log(10) * (T_phi1_vdb2018+T_phi1_upper_err_vdb2018) * (10**((x_vdb2018-(T_M_star_vdb2018+T_M_star_upper_err_vdb2018))*(1+(T_alpha1_vdb2018-T_alpha1_upper_err_vdb2018)))) * np.exp(-10**(x_vdb2018-(T_M_star_vdb2018+T_M_star_upper_err_vdb2018)))
# lower_error_bound
SF_model_lower_err_vdb2018 = np.log(10) * (SF_phi1_vdb2018-SF_phi1_lower_err_vdb2018) * (10**((x_vdb2018-(SF_M_star_vdb2018-SF_M_star_lower_err_vdb2018))*(1+(SF_alpha1_vdb2018+SF_alpha1_lower_err_vdb2018)))) * np.exp(-10**(x_vdb2018-(SF_M_star_vdb2018-SF_M_star_lower_err_vdb2018)))
Q_model_lower_err_vdb2018 = np.log(10) * (Q_phi1_vdb2018-Q_phi1_lower_err_vdb2018) * (10**((x_vdb2018-(Q_M_star_vdb2018-Q_M_star_lower_err_vdb2018))*(1+(Q_alpha1_vdb2018+Q_alpha1_lower_err_vdb2018)))) * np.exp(-10**(x_vdb2018-(Q_M_star_vdb2018-Q_M_star_lower_err_vdb2018)))
T_model_lower_err_vdb2018 = np.log(10) * (T_phi1_vdb2018-T_phi1_lower_err_vdb2018) * (10**((x_vdb2018-(T_M_star_vdb2018-T_M_star_lower_err_vdb2018))*(1+(T_alpha1_vdb2018+T_alpha1_lower_err_vdb2018)))) * np.exp(-10**(x_vdb2018-(T_M_star_vdb2018-T_M_star_lower_err_vdb2018)))
#


#
#
#
## SCETION (2.3): Visualization
#
if plot_flag_2 == 1:
    # Plot: compare environments by population (i.e. plot Cluster vs Field for Total, SF, & Q)
    # plt.close()
    SMF = plt.figure()
    gs = gridspec.GridSpec(1,3, wspace=0, hspace=0, width_ratios=[1,1,1])   #make a tiled-plot
    # Total population
    ax0 = plt.subplot(gs[0])
    #Plot Schechter fits:
    # my fits
    ax0.plot(x_plot,total_model_mcmc_plot,color='black',linestyle='-',label='this work: 0.25 < z < 0.75', linewidth=1.0)
    ax0.plot(x_plot,total_model_upper, color='grey', linestyle='--',linewidth=1.0,alpha=0.5)
    ax0.plot(x_plot,total_model_lower, color='grey', linestyle='--',linewidth=1.0,alpha=0.5)
    ax0.fill_between(x_plot,total_model_mcmc_plot,total_model_upper,facecolor='grey',interpolate=True,alpha=0.1)
    ax0.fill_between(x_plot,total_model_mcmc_plot,total_model_lower,facecolor='grey',interpolate=True,alpha=0.1)
    # vdb+2018
    ax0.plot(x_vdb2018,T_model_vdb2018, color='magenta', linestyle='-',label='vdB+2018: 0.5 < z < 0.7', linewidth=1.0)
    ax0.plot(x_vdb2018,T_model_upper_err_vdb2018, color='magenta', linestyle='--',linewidth=1.0,alpha=0.3)
    ax0.plot(x_vdb2018,T_model_lower_err_vdb2018, color='magenta', linestyle='--', linewidth=1.0,alpha=0.3)
    ax0.fill_between(x_vdb2018,T_model_vdb2018,T_model_upper_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
    ax0.fill_between(x_vdb2018,T_model_vdb2018,T_model_lower_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
    #
    ax0.set_xscale('linear')
    ax0.set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    ax0.set_xlim(7.8,12.5)
    ax0.set_yscale('log')
    ax0.set_ylim(5e-4,0.8)
    ax0.minorticks_on()
    ax0.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=False,labelleft=True,labelbottom=True,labelsize=18)
    ax0.yaxis.set_label_position("left")
    ax0.set_ylabel('${\phi}$ [$Mpc^{-3}$ $dex^{-1}$]',fontsize=20)
    ax0.set_title('Total')
    # ax0.legend(scatterpoints=1,loc='upper right', frameon=False, fontsize = 'medium')
    # ax0.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')
    #
    # SF population
    ax1 = plt.subplot(gs[1])
    # MY FITS
    ax1.plot(x_plot,SF_model_mcmc_plot,color='b',linestyle='-',linewidth=1.0)
    ax1.plot(x_plot,SF_model_upper,color='b',linestyle='--',linewidth=1.0,alpha=0.5)
    ax1.plot(x_plot,SF_model_lower,color='b',linestyle='--',linewidth=1.0,alpha=0.5)
    ax1.fill_between(x_plot,SF_model_mcmc_plot,SF_model_upper,facecolor='b',interpolate=True,alpha=0.1)
    ax1.fill_between(x_plot,SF_model_mcmc_plot,SF_model_lower,facecolor='b',interpolate=True,alpha=0.1)
    # vdb+2018
    ax1.plot(x_vdb2018,SF_model_vdb2018, color='magenta', linestyle='-',label='vdB+2018: 0.5 < z < 0.7', linewidth=1.0)
    ax1.plot(x_vdb2018,SF_model_upper_err_vdb2018, color='magenta', linestyle='--',linewidth=1.0,alpha=0.3)
    ax1.plot(x_vdb2018,SF_model_lower_err_vdb2018, color='magenta', linestyle='--', linewidth=1.0,alpha=0.3)
    ax1.fill_between(x_vdb2018,SF_model_vdb2018,SF_model_upper_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
    ax1.fill_between(x_vdb2018,SF_model_vdb2018,SF_model_lower_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
    # annunziatella+2014

    #
    ax1.set_xscale('linear')
    ax1.set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    ax1.set_xlim(7.8,12.5)
    ax1.set_yscale('log')
    ax1.set_ylim(5e-4,0.8)
    ax1.minorticks_on()
    ax1.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=False,labelleft=False,labelbottom=True,labelsize=18)
    # ax0.yaxis.set_label_position("left")
    # ax0.ylabel('???')
    ax1.set_title('Star-Forming')
    # ax1.legend(scatterpoints=1,loc='lower left', frameon=False, fontsize = 'x-small')
    # ax1.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')#
    # Q population
    #Plot Schechter fits:
    ax2 = plt.subplot(gs[2])
    ax2.plot(x_plot,Q_model_mcmc_plot,color='r',linestyle='-',label='this work: 0.25 < z < 0.75', linewidth=1.0)
    ax2.plot(x_plot,Q_model_upper,color='r',linestyle='--',linewidth=1.0,alpha=0.5)
    ax2.plot(x_plot,Q_model_lower,color='r',linestyle='--',linewidth=1.0,alpha=0.5)
    ax2.fill_between(x_plot,Q_model_mcmc_plot,Q_model_upper,facecolor='r',interpolate=True,alpha=0.1)
    ax2.fill_between(x_plot,Q_model_mcmc_plot,Q_model_lower,facecolor='r',interpolate=True,alpha=0.1)
    # vdb+2018
    ax2.plot(x_vdb2018,Q_model_vdb2018, color='magenta', linestyle='-',label='vdB+2018: 0.5 < z < 0.7', linewidth=1.0)
    ax2.plot(x_vdb2018,Q_model_upper_err_vdb2018, color='magenta', linestyle='--',linewidth=1.0,alpha=0.3)
    ax2.plot(x_vdb2018,Q_model_lower_err_vdb2018, color='magenta', linestyle='--', linewidth=1.0,alpha=0.3)
    ax2.fill_between(x_vdb2018,Q_model_vdb2018,Q_model_upper_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
    ax2.fill_between(x_vdb2018,Q_model_vdb2018,Q_model_lower_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
    # annunziatella+2014

    #
    ax2.set_xscale('linear')
    ax2.set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    ax2.set_xlim(7.8,12.5)
    ax2.set_yscale('log')
    ax2.set_ylim(5e-4,0.8)
    ax2.minorticks_on()
    ax2.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=True, labelleft=False,labelbottom=True,labelsize=18)
    ax2.yaxis.set_label_position("right")
    ax2.set_ylabel('${\phi}$ [$Mpc^{-3}$ $dex^{-1}$]',fontsize=20)
    ax2.set_title('Quiescent')
    ax2.legend(scatterpoints=1,loc='upper right', frameon=False, fontsize = 'medium')
    # ax2.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')#
    #
    plt.show()
    #
#
#
#
## SECTION (3): paneled plot comparing individual studies from literature to our results
#
if plot_flag_3 == 99:
    #
    ## TOTAL population
    fig, axs = plt.subplots(2,2, sharex=True,sharey=True,gridspec_kw={'hspace': 0,'wspace': 0})
    #
    fig.suptitle('Field: Total',fontsize=30)
    # Panel 1: Baldry+2012
    # my fits
    axs[0,0].plot(x_plot,total_model_field_mcmc_plot_double,color='black',linestyle='-',label='this work: 0.25 < z < 0.75', linewidth=1.0)
    axs[0,0].plot(x_plot,total_model_field_upper, color='black', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].plot(x_plot,total_model_field_lower, color='black', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_upper,facecolor='black',interpolate=True,alpha=0.1)
    axs[0,0].fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_lower,facecolor='black',interpolate=True,alpha=0.1)
    # baldry+2012 data pts + schechter fits
    axs[0,0].errorbar(x_baldry2012_data,total_smf_baldry2012,yerr=total_err_baldry2012,marker='s',color='navy',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[0,0].plot(x_baldry2012,T_model_baldry2012, color='navy', linestyle='-',label='baldry+2012: z < 0.06\n(GAMA)', linewidth=1.0)
    axs[0,0].plot(x_baldry2012,T_model_upper_err_baldry2012, color='navy', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].plot(x_baldry2012,T_model_lower_err_baldry2012, color='navy', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,0].fill_between(x_baldry2012,T_model_baldry2012,T_model_upper_err_baldry2012,facecolor='navy',interpolate=True,alpha=0.1)
    axs[0,0].fill_between(x_baldry2012,T_model_baldry2012,T_model_lower_err_baldry2012,facecolor='navy',interpolate=True,alpha=0.1)
    # fig parameters
    axs[0,0].set_xscale('linear')
    axs[0,0].set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    axs[0,0].set_xlim(7.8,12.5)
    axs[0,0].set_yscale('log')
    axs[0,0].set_ylim(1e-6,0.2)
    axs[0,0].minorticks_on()
    axs[0,0].tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=False, labelleft=True,labelbottom=False,labelsize=18)
    axs[0,0].yaxis.set_label_position("right")
    axs[0,0].set_ylabel('${\phi}$ [$Mpc^{-3}$ $dex^{-1}$]',fontsize=20)
    # axs[0,0].set_title('Quiescent')
    axs[0,0].legend(scatterpoints=1,loc='upper right', frameon=False, fontsize = 'medium')
    #
    # Panel 2: Muzzin+2013
    # my fits
    axs[0,1].plot(x_plot,total_model_field_mcmc_plot_double,color='black',linestyle='-',label='this work: 0.25 < z < 0.75', linewidth=1.0)
    axs[0,1].plot(x_plot,total_model_field_upper, color='black', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_plot,total_model_field_lower, color='black', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_upper,facecolor='black',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_lower,facecolor='black',interpolate=True,alpha=0.1)
    # muzzin+2013 data + fits
    axs[0,1].plot(x_muzzin2013_lo_z,T_model_muzzin_lo_z, color='goldenrod', linestyle='-',label='muzzin+2013: 0.2 < z < 0.5', linewidth=1.0)
    axs[0,1].plot(x_muzzin2013_lo_z,T_model_muzzin_upper_lo_z, color='goldenrod', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_muzzin2013_lo_z,T_model_muzzin_lower_lo_z, color='goldenrod', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_muzzin2013_lo_z,T_model_muzzin_lo_z,T_model_muzzin_upper_lo_z,facecolor='goldenrod',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_muzzin2013_lo_z,T_model_muzzin_lo_z,T_model_muzzin_lower_lo_z,facecolor='goldenrod',interpolate=True,alpha=0.1)
    axs[0,1].plot(x_muzzin2013_hi_z,T_model_muzzin_hi_z, color='darkgoldenrod', linestyle='-',label='muzzin+2013: 0.5 < z < 1.0\n(COSMOS/UltraVISTA)', linewidth=1.0)
    axs[0,1].plot(x_muzzin2013_hi_z,T_model_muzzin_upper_hi_z, color='darkgoldenrod', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_muzzin2013_hi_z,T_model_muzzin_lower_hi_z, color='darkgoldenrod', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_muzzin2013_hi_z,T_model_muzzin_hi_z,T_model_muzzin_upper_hi_z,facecolor='darkgoldenrod',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_muzzin2013_hi_z,T_model_muzzin_hi_z,T_model_muzzin_lower_hi_z,facecolor='darkgoldenrod',interpolate=True,alpha=0.1)
    # fig parameters
    axs[0,1].set_xscale('linear')
    axs[0,1].set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    axs[0,1].set_xlim(7.8,12.5)
    axs[0,1].set_yscale('log')
    axs[0,1].set_ylim(1e-6,0.2)
    axs[0,1].minorticks_on()
    axs[0,1].tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=True, labelleft=False,labelbottom=False,labelsize=18)
    axs[0,1].yaxis.set_label_position("right")
    axs[0,1].set_ylabel('${\phi}$ [$Mpc^{-3}$ $dex^{-1}$]',fontsize=20)
    # axs[0,1].set_title('Quiescent')
    axs[0,1].legend(scatterpoints=1,loc='upper right', frameon=False, fontsize = 'medium')
    #
    # Panel 3: Tomczak+2014
    # my fits
    axs[1,0].plot(x_plot,total_model_field_mcmc_plot_double,color='black',linestyle='-',label='this work: 0.25 < z < 0.75', linewidth=1.0)
    axs[1,0].plot(x_plot,total_model_field_upper, color='black', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_plot,total_model_field_lower, color='black', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_upper,facecolor='black',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_lower,facecolor='black',interpolate=True,alpha=0.1)
    # tomczak+2014 data + fits
    axs[1,0].errorbar(x_tom2014_data,total_smf_tom2014_lo_z,yerr=total_err_tom2014_lo_z,marker='<',color='lime',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[1,0].errorbar(x_tom2014_data,total_smf_tom2014_hi_z,yerr=total_err_tom2014_hi_z,marker='>',color='darkgreen',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[1,0].plot(x_tomczak2014_lo_z,T_model_tomczak_lo_z, color='lime', linestyle='-',label='tomczak+2014: 0.2 < z < 0.5', linewidth=1.0)
    axs[1,0].plot(x_tomczak2014_lo_z,T_model_tomczak_upper_lo_z, color='lime', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_tomczak2014_lo_z,T_model_tomczak_lower_lo_z, color='lime', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_tomczak2014_lo_z,T_model_tomczak_lo_z,T_model_tomczak_upper_lo_z,facecolor='lime',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_tomczak2014_lo_z,T_model_tomczak_lo_z,T_model_tomczak_lower_lo_z,facecolor='lime',interpolate=True,alpha=0.1)
    axs[1,0].plot(x_tomczak2014_hi_z,T_model_tomczak_hi_z, color='darkgreen', linestyle='-',label='tomczak+2014: 0.5 < z < 0.75\n(ZFOURGE/CANDELS)', linewidth=1.0)
    axs[1,0].plot(x_tomczak2014_hi_z,T_model_tomczak_upper_hi_z, color='darkgreen', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_tomczak2014_hi_z,T_model_tomczak_lower_hi_z, color='darkgreen', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_tomczak2014_hi_z,T_model_tomczak_hi_z,T_model_tomczak_upper_hi_z,facecolor='darkgreen',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_tomczak2014_hi_z,T_model_tomczak_hi_z,T_model_tomczak_lower_hi_z,facecolor='darkgreen',interpolate=True,alpha=0.1)
    # fig parameters
    axs[1,0].set_xscale('linear')
    axs[1,0].set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    axs[1,0].set_xlim(7.8,12.5)
    axs[1,0].set_yscale('log')
    axs[1,0].set_ylim(1e-6,0.2)
    axs[1,0].minorticks_on()
    axs[1,0].tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=False, labelleft=True,labelbottom=True,labelsize=18)
    axs[1,0].yaxis.set_label_position("right")
    axs[1,0].set_ylabel('${\phi}$ [$Mpc^{-3}$ $dex^{-1}$]',fontsize=20)
    # axs[1,0].set_title('Quiescent')
    axs[1,0].legend(scatterpoints=1,loc='upper right', frameon=False, fontsize = 'medium')
    #
    # Panel 4: mcleod+2020
    # my fits
    axs[1,1].plot(x_plot,total_model_field_mcmc_plot_double,color='black',linestyle='-',label='this work: 0.25 < z < 0.75', linewidth=1.0)
    axs[1,1].plot(x_plot,total_model_field_upper, color='black', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,1].plot(x_plot,total_model_field_lower, color='black', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,1].fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_upper,facecolor='black',interpolate=True,alpha=0.1)
    axs[1,1].fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_lower,facecolor='black',interpolate=True,alpha=0.1)
    # mcleod+2020 data pts + schechter fits
    axs[1,1].errorbar(x_mcleod2020_data,total_smf_mcleod2020,yerr=total_err_mcleod2020,marker='s',color='magenta',label='mcleod+2020: 0.25 < z < 0.75\n(UKIDSS UDS/COSMOS/CFHTLS-D1)',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    # axs[1,1].plot(x_mcleod2020,T_model_mcleod2020, color='magenta', linestyle='-',label='mcleod+2020: 0.25 < z < 0.75\n(UKIDSS UDS/COSMOS/CFHTLS-D1)', linewidth=1.0)
    # axs[1,1].plot(x_mcleod2020,T_model_upper_err_mcleod2020, color='magenta', linestyle='--',linewidth=1.0,alpha=0.3)
    # axs[1,1].plot(x_mcleod2020,T_model_lower_err_mcleod2020, color='magenta', linestyle='--', linewidth=1.0,alpha=0.3)
    # axs[1,1].fill_between(x_mcleod2020,T_model_mcleod2020,T_model_upper_err_mcleod2020,facecolor='magenta',interpolate=True,alpha=0.1)
    # axs[1,1].fill_between(x_mcleod2020,T_model_mcleod2020,T_model_lower_err_mcleod2020,facecolor='magenta',interpolate=True,alpha=0.1)
    # fig parameters
    axs[1,1].set_xscale('linear')
    axs[1,1].set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    axs[1,1].set_xlim(7.8,12.5)
    axs[1,1].set_yscale('log')
    axs[1,1].set_ylim(1e-6,0.2)
    axs[1,1].minorticks_on()
    axs[1,1].tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=True, labelleft=False,labelbottom=True,labelsize=18)
    axs[1,1].yaxis.set_label_position("right")
    axs[1,1].set_ylabel('${\phi}$ [$Mpc^{-3}$ $dex^{-1}$]',fontsize=20)
    # axs[1,1].set_title('Quiescent')
    axs[1,1].legend(scatterpoints=1,loc='upper right', frameon=False, fontsize = 'medium')
    #
    #
    #
    #
    ## SF population
    #
    fig, axs = plt.subplots(2,2, sharex=True,sharey=True,gridspec_kw={'hspace': 0,'wspace': 0})
    #
    fig.suptitle('Field: Star-forming',fontsize=30)
    # Panel 1: Baldry+2012
    # my fits
    axs[0,0].plot(x_plot,SF_model_field_mcmc_plot,color='b',linestyle='-',label='this work: 0.25 < z < 0.75', linewidth=1.0)
    axs[0,0].plot(x_plot,SF_model_field_upper, color='b', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].plot(x_plot,SF_model_field_lower, color='b', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_upper,facecolor='b',interpolate=True,alpha=0.1)
    axs[0,0].fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_lower,facecolor='b',interpolate=True,alpha=0.1)
    # baldry+2012 data pts + schechter fits
    # axs[0,0].errorbar(x_baldry2012_data,total_smf_baldry2012,yerr=total_err_baldry2012,marker='s',color='navy',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[0,0].errorbar(x_baldry2012_data_SF,SF_smf_baldry2012,yerr=SF_err_baldry2012,marker='s',color='navy',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[0,0].plot(x_baldry2012,SF_model_baldry2012, color='navy', linestyle='-',label='baldry+2012: z < 0.06\n(GAMA)', linewidth=1.0)
    axs[0,0].plot(x_baldry2012,SF_model_upper_err_baldry2012, color='navy', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].plot(x_baldry2012,SF_model_lower_err_baldry2012, color='navy', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,0].fill_between(x_baldry2012,SF_model_baldry2012,SF_model_upper_err_baldry2012,facecolor='navy',interpolate=True,alpha=0.1)
    axs[0,0].fill_between(x_baldry2012,SF_model_baldry2012,SF_model_lower_err_baldry2012,facecolor='navy',interpolate=True,alpha=0.1)
    # fig parameters
    axs[0,0].set_xscale('linear')
    axs[0,0].set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    axs[0,0].set_xlim(7.8,12.5)
    axs[0,0].set_yscale('log')
    axs[0,0].set_ylim(1e-6,0.2)
    axs[0,0].minorticks_on()
    axs[0,0].tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=False, labelleft=True,labelbottom=False,labelsize=18)
    axs[0,0].yaxis.set_label_position("right")
    axs[0,0].set_ylabel('${\phi}$ [$Mpc^{-3}$ $dex^{-1}$]',fontsize=20)
    # axs[0,0].set_title('Quiescent')
    axs[0,0].legend(scatterpoints=1,loc='upper right', frameon=False, fontsize = 'medium')
    #
    # Panel 2: Muzzin+2013
    # my fits
    axs[0,1].plot(x_plot,SF_model_field_mcmc_plot,color='b',linestyle='-',label='this work: 0.25 < z < 0.75', linewidth=1.0)
    axs[0,1].plot(x_plot,SF_model_field_upper, color='b', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_plot,SF_model_field_lower, color='b', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_upper,facecolor='b',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_lower,facecolor='b',interpolate=True,alpha=0.1)
    # muzzin+2013 data + fits
    axs[0,1].plot(x_muzzin2013_lo_z,SF_model_muzzin_lo_z, color='goldenrod', linestyle='-',label='muzzin+2013: 0.2 < z < 0.5', linewidth=1.0)
    axs[0,1].plot(x_muzzin2013_lo_z,SF_model_muzzin_upper_lo_z, color='goldenrod', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_muzzin2013_lo_z,SF_model_muzzin_lower_lo_z, color='goldenrod', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_muzzin2013_lo_z,SF_model_muzzin_lo_z,SF_model_muzzin_upper_lo_z,facecolor='goldenrod',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_muzzin2013_lo_z,SF_model_muzzin_lo_z,SF_model_muzzin_lower_lo_z,facecolor='goldenrod',interpolate=True,alpha=0.1)
    axs[0,1].plot(x_muzzin2013_hi_z,SF_model_muzzin_hi_z, color='darkgoldenrod', linestyle='-',label='muzzin+2013: 0.5 < z < 1.0\n(COSMOS/UltraVISTA)', linewidth=1.0)
    axs[0,1].plot(x_muzzin2013_hi_z,SF_model_muzzin_upper_hi_z, color='darkgoldenrod', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_muzzin2013_hi_z,SF_model_muzzin_lower_hi_z, color='darkgoldenrod', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_muzzin2013_hi_z,SF_model_muzzin_hi_z,SF_model_muzzin_upper_hi_z,facecolor='darkgoldenrod',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_muzzin2013_hi_z,SF_model_muzzin_hi_z,SF_model_muzzin_lower_hi_z,facecolor='darkgoldenrod',interpolate=True,alpha=0.1)
    # fig parameters
    axs[0,1].set_xscale('linear')
    axs[0,1].set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    axs[0,1].set_xlim(7.8,12.5)
    axs[0,1].set_yscale('log')
    axs[0,1].set_ylim(1e-6,0.2)
    axs[0,1].minorticks_on()
    axs[0,1].tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=True, labelleft=False,labelbottom=False,labelsize=18)
    axs[0,1].yaxis.set_label_position("right")
    axs[0,1].set_ylabel('${\phi}$ [$Mpc^{-3}$ $dex^{-1}$]',fontsize=20)
    # axs[0,1].set_title('Quiescent')
    axs[0,1].legend(scatterpoints=1,loc='upper right', frameon=False, fontsize = 'medium')
    #
    # Panel 3: Tomczak+2014
    # my fits
    axs[1,0].plot(x_plot,SF_model_field_mcmc_plot,color='b',linestyle='-',label='this work: 0.25 < z < 0.75', linewidth=1.0)
    axs[1,0].plot(x_plot,SF_model_field_upper, color='b', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_plot,SF_model_field_lower, color='b', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_upper,facecolor='b',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_lower,facecolor='b',interpolate=True,alpha=0.1)
    # tomczak+2014 data + fits
    axs[1,0].errorbar(x_tom2014_data,SF_smf_tom2014_lo_z,yerr=SF_err_tom2014_lo_z,marker='<',color='lime',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[1,0].errorbar(x_tom2014_data,SF_smf_tom2014_hi_z,yerr=SF_err_tom2014_hi_z,marker='>',color='darkgreen',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[1,0].plot(x_tomczak2014_lo_z,SF_model_tomczak_lo_z, color='lime', linestyle='-',label='tomczak+2014: 0.2 < z < 0.5', linewidth=1.0)
    axs[1,0].plot(x_tomczak2014_lo_z,SF_model_tomczak_upper_lo_z, color='lime', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_tomczak2014_lo_z,SF_model_tomczak_lower_lo_z, color='lime', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_tomczak2014_lo_z,SF_model_tomczak_lo_z,SF_model_tomczak_upper_lo_z,facecolor='lime',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_tomczak2014_lo_z,SF_model_tomczak_lo_z,SF_model_tomczak_lower_lo_z,facecolor='lime',interpolate=True,alpha=0.1)
    axs[1,0].plot(x_tomczak2014_hi_z,SF_model_tomczak_hi_z, color='darkgreen', linestyle='-',label='tomczak+2014: 0.5 < z < 0.75\n(ZFOURGE/CANDELS)', linewidth=1.0)
    axs[1,0].plot(x_tomczak2014_hi_z,SF_model_tomczak_upper_hi_z, color='darkgreen', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_tomczak2014_hi_z,SF_model_tomczak_lower_hi_z, color='darkgreen', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_tomczak2014_hi_z,SF_model_tomczak_hi_z,SF_model_tomczak_upper_hi_z,facecolor='darkgreen',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_tomczak2014_hi_z,SF_model_tomczak_hi_z,SF_model_tomczak_lower_hi_z,facecolor='darkgreen',interpolate=True,alpha=0.1)
    # fig parameters
    axs[1,0].set_xscale('linear')
    axs[1,0].set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    axs[1,0].set_xlim(7.8,12.5)
    axs[1,0].set_yscale('log')
    axs[1,0].set_ylim(1e-6,0.2)
    axs[1,0].minorticks_on()
    axs[1,0].tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=False, labelleft=True,labelbottom=True,labelsize=18)
    axs[1,0].yaxis.set_label_position("right")
    axs[1,0].set_ylabel('${\phi}$ [$Mpc^{-3}$ $dex^{-1}$]',fontsize=20)
    # axs[1,0].set_title('Quiescent')
    axs[1,0].legend(scatterpoints=1,loc='upper right', frameon=False, fontsize = 'medium')
    #
    # Panel 4: mcleod+2020
    # my fits
    axs[1,1].plot(x_plot,SF_model_field_mcmc_plot,color='b',linestyle='-',label='this work: 0.25 < z < 0.75', linewidth=1.0)
    axs[1,1].plot(x_plot,SF_model_field_upper, color='b', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,1].plot(x_plot,SF_model_field_lower, color='b', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,1].fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_upper,facecolor='b',interpolate=True,alpha=0.1)
    axs[1,1].fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_lower,facecolor='b',interpolate=True,alpha=0.1)
    # mcleod+2020 data pts + schechter fits
    axs[1,1].errorbar(x_mcleod2020_data,SF_smf_mcleod2020,yerr=total_err_mcleod2020,marker='s',color='magenta',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[1,1].plot(x_mcleod2020,SF_model_mcleod2020, color='magenta', linestyle='-',label='mcleod+2020: 0.25 < z < 0.75\n(UKIDSS UDS/COSMOS/CFHTLS-D1)', linewidth=1.0)
    axs[1,1].plot(x_mcleod2020,SF_model_upper_err_mcleod2020, color='magenta', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,1].plot(x_mcleod2020,SF_model_lower_err_mcleod2020, color='magenta', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,1].fill_between(x_mcleod2020,SF_model_mcleod2020,SF_model_upper_err_mcleod2020,facecolor='magenta',interpolate=True,alpha=0.1)
    axs[1,1].fill_between(x_mcleod2020,SF_model_mcleod2020,SF_model_lower_err_mcleod2020,facecolor='magenta',interpolate=True,alpha=0.1)
    # fig parameters
    axs[1,1].set_xscale('linear')
    axs[1,1].set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    axs[1,1].set_xlim(7.8,12.5)
    axs[1,1].set_yscale('log')
    axs[1,1].set_ylim(1e-6,0.2)
    axs[1,1].minorticks_on()
    axs[1,1].tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=True, labelleft=False,labelbottom=True,labelsize=18)
    axs[1,1].yaxis.set_label_position("right")
    axs[1,1].set_ylabel('${\phi}$ [$Mpc^{-3}$ $dex^{-1}$]',fontsize=20)
    # axs[1,1].set_title('Quiescent')
    axs[1,1].legend(scatterpoints=1,loc='upper right', frameon=False, fontsize = 'medium')
    #
    #
    #
    ## Q population
    #
    fig, axs = plt.subplots(2,2, sharex=True,sharey=True,gridspec_kw={'hspace': 0,'wspace': 0})
    #
    fig.suptitle('Field: Quiescent',fontsize=30)
    # Panel 1: Baldry+2012
    # my fits
    axs[0,0].plot(x_plot,Q_model_field_mcmc_plot_double,color='r',linestyle='-',label='this work: 0.25 < z < 0.75', linewidth=1.0)
    axs[0,0].plot(x_plot,Q_model_field_upper, color='r', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].plot(x_plot,Q_model_field_lower, color='r', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_upper,facecolor='r',interpolate=True,alpha=0.1)
    axs[0,0].fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_lower,facecolor='r',interpolate=True,alpha=0.1)
    # baldry+2012 data pts + schechter fits
    # axs[0,0].errorbar(x_baldry2012_data,total_smf_baldry2012,yerr=total_err_baldry2012,marker='s',color='navy',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[0,0].errorbar(x_baldry2012_data_Q,Q_smf_baldry2012,yerr=Q_err_baldry2012,marker='s',color='navy',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[0,0].plot(x_baldry2012_Q,Q_model_baldry2012, color='navy', linestyle='-',label='baldry+2012: z < 0.06\n(GAMA)', linewidth=1.0)
    axs[0,0].plot(x_baldry2012_Q,Q_model_upper_err_baldry2012, color='navy', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].plot(x_baldry2012_Q,Q_model_lower_err_baldry2012, color='navy', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,0].fill_between(x_baldry2012_Q,Q_model_baldry2012,Q_model_upper_err_baldry2012,facecolor='navy',interpolate=True,alpha=0.1)
    axs[0,0].fill_between(x_baldry2012_Q,Q_model_baldry2012,Q_model_lower_err_baldry2012,facecolor='navy',interpolate=True,alpha=0.1)
    # fig parameters
    axs[0,0].set_xscale('linear')
    axs[0,0].set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    axs[0,0].set_xlim(7.8,12.5)
    axs[0,0].set_yscale('log')
    axs[0,0].set_ylim(1e-6,0.2)
    axs[0,0].minorticks_on()
    axs[0,0].tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=False, labelleft=True,labelbottom=False,labelsize=18)
    axs[0,0].yaxis.set_label_position("right")
    axs[0,0].set_ylabel('${\phi}$ [$Mpc^{-3}$ $dex^{-1}$]',fontsize=20)
    # axs[0,0].set_title('Quiescent')
    axs[0,0].legend(scatterpoints=1,loc='upper right', frameon=False, fontsize = 'medium')
    #
    # Panel 2: Muzzin+2013
    # my fits
    axs[0,1].plot(x_plot,Q_model_field_mcmc_plot_double,color='r',linestyle='-',label='this work: 0.25 < z < 0.75', linewidth=1.0)
    axs[0,1].plot(x_plot,Q_model_field_upper, color='r', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_plot,Q_model_field_lower, color='r', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_upper,facecolor='r',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_lower,facecolor='r',interpolate=True,alpha=0.1)
    # muzzin+2013 data + fits
    axs[0,1].plot(x_muzzin2013_lo_z,Q_model_muzzin_lo_z, color='goldenrod', linestyle='-',label='muzzin+2013: 0.2 < z < 0.5', linewidth=1.0)
    axs[0,1].plot(x_muzzin2013_lo_z,Q_model_muzzin_upper_lo_z, color='goldenrod', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_muzzin2013_lo_z,Q_model_muzzin_lower_lo_z, color='goldenrod', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_muzzin2013_lo_z,Q_model_muzzin_lo_z,Q_model_muzzin_upper_lo_z,facecolor='goldenrod',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_muzzin2013_lo_z,Q_model_muzzin_lo_z,Q_model_muzzin_lower_lo_z,facecolor='goldenrod',interpolate=True,alpha=0.1)
    axs[0,1].plot(x_muzzin2013_hi_z,Q_model_muzzin_hi_z, color='darkgoldenrod', linestyle='-',label='muzzin+2013: 0.5 < z < 1.0\n(COSMOS/UltraVISTA)', linewidth=1.0)
    axs[0,1].plot(x_muzzin2013_hi_z,Q_model_muzzin_upper_hi_z, color='darkgoldenrod', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_muzzin2013_hi_z,Q_model_muzzin_lower_hi_z, color='darkgoldenrod', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_muzzin2013_hi_z,Q_model_muzzin_hi_z,Q_model_muzzin_upper_hi_z,facecolor='darkgoldenrod',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_muzzin2013_hi_z,Q_model_muzzin_hi_z,Q_model_muzzin_lower_hi_z,facecolor='darkgoldenrod',interpolate=True,alpha=0.1)
    # fig parameters
    axs[0,1].set_xscale('linear')
    axs[0,1].set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    axs[0,1].set_xlim(7.8,12.5)
    axs[0,1].set_yscale('log')
    axs[0,1].set_ylim(1e-6,0.2)
    axs[0,1].minorticks_on()
    axs[0,1].tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=True, labelleft=False,labelbottom=False,labelsize=18)
    axs[0,1].yaxis.set_label_position("right")
    axs[0,1].set_ylabel('${\phi}$ [$Mpc^{-3}$ $dex^{-1}$]',fontsize=20)
    # axs[0,1].set_title('Quiescent')
    axs[0,1].legend(scatterpoints=1,loc='upper right', frameon=False, fontsize = 'medium')
    #
    # Panel 3: Tomczak+2014
    # my fits
    axs[1,0].plot(x_plot,Q_model_field_mcmc_plot_double,color='r',linestyle='-',label='this work: 0.25 < z < 0.75', linewidth=1.0)
    axs[1,0].plot(x_plot,Q_model_field_upper, color='r', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_plot,Q_model_field_lower, color='r', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_upper,facecolor='r',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_lower,facecolor='r',interpolate=True,alpha=0.1)
    # tomczak+2014 data + fits
    axs[1,0].errorbar(x_tom2014_data,Q_smf_tom2014_lo_z,yerr=Q_err_tom2014_lo_z,marker='<',color='lime',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[1,0].errorbar(x_tom2014_data,Q_smf_tom2014_hi_z,yerr=Q_err_tom2014_hi_z,marker='>',color='darkgreen',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[1,0].plot(x_tomczak2014_lo_z_Q,Q_model_tomczak_lo_z, color='lime', linestyle='-',label='tomczak+2014: 0.2 < z < 0.5', linewidth=1.0)
    axs[1,0].plot(x_tomczak2014_lo_z_Q,Q_model_tomczak_upper_lo_z, color='lime', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_tomczak2014_lo_z_Q,Q_model_tomczak_lower_lo_z, color='lime', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_tomczak2014_lo_z_Q,Q_model_tomczak_lo_z,Q_model_tomczak_upper_lo_z,facecolor='lime',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_tomczak2014_lo_z_Q,Q_model_tomczak_lo_z,Q_model_tomczak_lower_lo_z,facecolor='lime',interpolate=True,alpha=0.1)
    axs[1,0].plot(x_tomczak2014_hi_z_Q,Q_model_tomczak_hi_z, color='darkgreen', linestyle='-',label='tomczak+2014: 0.5 < z < 0.75\n(ZFOURGE/CANDELS)', linewidth=1.0)
    axs[1,0].plot(x_tomczak2014_hi_z_Q,Q_model_tomczak_upper_hi_z, color='darkgreen', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_tomczak2014_hi_z_Q,Q_model_tomczak_lower_hi_z, color='darkgreen', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_tomczak2014_hi_z_Q,Q_model_tomczak_hi_z,Q_model_tomczak_upper_hi_z,facecolor='darkgreen',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_tomczak2014_hi_z_Q,Q_model_tomczak_hi_z,Q_model_tomczak_lower_hi_z,facecolor='darkgreen',interpolate=True,alpha=0.1)
    # fig parameters
    axs[1,0].set_xscale('linear')
    axs[1,0].set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    axs[1,0].set_xlim(7.8,12.5)
    axs[1,0].set_yscale('log')
    axs[1,0].set_ylim(1e-6,0.2)
    axs[1,0].minorticks_on()
    axs[1,0].tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=False, labelleft=True,labelbottom=True,labelsize=18)
    axs[1,0].yaxis.set_label_position("right")
    axs[1,0].set_ylabel('${\phi}$ [$Mpc^{-3}$ $dex^{-1}$]',fontsize=20)
    # axs[1,0].set_title('Quiescent')
    axs[1,0].legend(scatterpoints=1,loc='upper right', frameon=False, fontsize = 'medium')
    #
    # Panel 4: mcleod+2020
    # my fits
    axs[1,1].plot(x_plot,Q_model_field_mcmc_plot_double,color='r',linestyle='-',label='this work: 0.25 < z < 0.75', linewidth=1.0)
    axs[1,1].plot(x_plot,Q_model_field_upper, color='r', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,1].plot(x_plot,Q_model_field_lower, color='r', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,1].fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_upper,facecolor='r',interpolate=True,alpha=0.1)
    axs[1,1].fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_lower,facecolor='r',interpolate=True,alpha=0.1)
    # mcleod+2020 data pts + schechter fits
    axs[1,1].errorbar(x_mcleod2020_data,Q_smf_mcleod2020,yerr=total_err_mcleod2020,marker='s',color='magenta',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[1,1].plot(x_mcleod2020_Q,Q_model_mcleod2020, color='magenta', linestyle='-',label='mcleod+2020: 0.25 < z < 0.75\n(UKIDSS UDS/COSMOS/CFHTLS-D1)', linewidth=1.0)
    axs[1,1].plot(x_mcleod2020_Q,Q_model_upper_err_mcleod2020, color='magenta', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,1].plot(x_mcleod2020_Q,Q_model_lower_err_mcleod2020, color='magenta', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,1].fill_between(x_mcleod2020_Q,Q_model_mcleod2020,Q_model_upper_err_mcleod2020,facecolor='magenta',interpolate=True,alpha=0.1)
    axs[1,1].fill_between(x_mcleod2020_Q,Q_model_mcleod2020,Q_model_lower_err_mcleod2020,facecolor='magenta',interpolate=True,alpha=0.1)
    # fig parameters
    axs[1,1].set_xscale('linear')
    axs[1,1].set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    axs[1,1].set_xlim(7.8,12.5)
    axs[1,1].set_yscale('log')
    axs[1,1].set_ylim(1e-6,0.2)
    axs[1,1].minorticks_on()
    axs[1,1].tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=True, labelleft=False,labelbottom=True,labelsize=18)
    axs[1,1].yaxis.set_label_position("right")
    axs[1,1].set_ylabel('${\phi}$ [$Mpc^{-3}$ $dex^{-1}$]',fontsize=20)
    # axs[1,1].set_title('Quiescent')
    axs[1,1].legend(scatterpoints=1,loc='upper right', frameon=False, fontsize = 'medium')
    #

    #
    #
    plt.show()
#
#
#
#
print('\n\n"smf_comparison_cluster.py"  terminated successfully.\n')
#
#
#
###################     PROGRAM END     ###################
