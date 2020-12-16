# Created on Sun Dec 13 13:55:31 2020
# Author: Ghassan T Sarrouh
#
########################### smf_comparison.py ###########################
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
plot_flag_1 = 0
plot_flag_2 = 0
plot_flag_3 = 1

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
## from Tomczak et al. 2014: (lo_z = 0.2 < z < 0.5; hi_z = 0.5 < z < 0.75)
#
# x-values (mass bin midpoints)
x_tom2014_data = np.arange(8,11.75,0.25)
#
# SF population
# lo z
SF_smf_tom2014_lo_z = np.array([1.42,1.59,1.76,1.91,2.08,2.06,2.17,2.25,2.36,2.50,2.63,2.91,3.43,4.39,0])
SF_lower_err_tom2014_lo_z = np.array([0.06,0.06,0.07,0.07,0.08,0.07,0.07,0.08,0.08,0.08,0.09,0.10,0.13,0.30,0])
SF_upper_err_tom2014_lo_z = np.array([0.07,0.07,0.08,0.08,0.09,0.08,0.09,0.10,0.10,0.09,0.11,0.12,0.18,0.41,0])
# in unit of log(phi) - i.e. 10^ - (SF_smf_tom2014)
SF_upper_err_tom2014_lo_z = np.abs(10**(-1 * (SF_smf_tom2014_lo_z - SF_upper_err_tom2014_lo_z)) - 10**(-1 * (SF_smf_tom2014_lo_z)))
SF_lower_err_tom2014_lo_z = np.abs(10**(-1 * (SF_smf_tom2014_lo_z + SF_lower_err_tom2014_lo_z)) - 10**(-1 * (SF_smf_tom2014_lo_z)))
SF_smf_tom2014_lo_z = 10**(-1 * SF_smf_tom2014_lo_z)
SF_err_tom2014_lo_z = np.array([SF_lower_err_tom2014_lo_z, SF_upper_err_tom2014_lo_z])
# hiz
SF_smf_tom2014_hi_z = np.array([0,1.6,1.67,1.83,1.92,2.09,2.19,2.28,2.39,2.55,2.76,3.00,3.46,4.3,0])
SF_lower_err_tom2014_hi_z = np.array([0,0.06,0.05,0.06,0.06,0.06,0.07,0.06,0.07,0.07,0.08,0.08,0.10,0.20,0])
SF_upper_err_tom2014_hi_z = np.array([0,0.07,0.06,0.06,0.07,0.07,0.08,0.07,0.08,0.08,0.09,0.10,0.13,0.25,0])
# in unit of log(phi) - i.e. 10^ - (SF_smf_tom2014)
SF_upper_err_tom2014_hi_z = np.abs(10**(-1 * (SF_smf_tom2014_hi_z - SF_upper_err_tom2014_hi_z)) - 10**(-1 * (SF_smf_tom2014_hi_z)))
SF_lower_err_tom2014_hi_z = np.abs(10**(-1 * (SF_smf_tom2014_hi_z + SF_lower_err_tom2014_hi_z)) - 10**(-1 * (SF_smf_tom2014_hi_z)))
SF_smf_tom2014_hi_z = 10**(-1 * SF_smf_tom2014_hi_z)
SF_err_tom2014_hi_z = np.array([SF_lower_err_tom2014_hi_z, SF_upper_err_tom2014_hi_z])
# Q population
# lo z
Q_smf_tom2014_lo_z = np.array([0,2.41,2.62,2.82,2.96,2.96,2.98,2.91,2.86,2.78,2.80,2.76,3.07,3.52,0])
Q_lower_err_tom2014_lo_z = np.array([0,0.08,0.10,0.12,0.14,0.08,0.09,0.09,0.09,0.08,0.09,0.09,0.12,0.14,0])
Q_upper_err_tom2014_lo_z = np.array([0,0.10,0.11,0.14,0.16,0.10,0.10,0.11,0.11,0.10,0.11,0.12,0.16,0.19,0])
# in unit of log(phi) - i.e. 10^ - (SF_smf_tom2014)
Q_upper_err_tom2014_lo_z = np.abs(10**(-1 * (Q_smf_tom2014_lo_z - Q_upper_err_tom2014_lo_z)) - 10**(-1 * (Q_smf_tom2014_lo_z)))
Q_lower_err_tom2014_lo_z = np.abs(10**(-1 * (Q_smf_tom2014_lo_z + Q_lower_err_tom2014_lo_z)) - 10**(-1 * (Q_smf_tom2014_lo_z)))
Q_smf_tom2014_lo_z = 10**(-1 * Q_smf_tom2014_lo_z)
Q_err_tom2014_lo_z = np.array([Q_lower_err_tom2014_lo_z, Q_upper_err_tom2014_lo_z])
# hi z
Q_smf_tom2014_hi_z = np.array([0,0,2.42,2.58,2.77,2.75,2.94,2.99,2.83,2.78,2.75,2.75,2.93,3.37,0])
Q_lower_err_tom2014_hi_z = np.array([0,0,0.07,0.07,0.09,0.09,0.10,0.07,0.07,0.07,0.08,0.08,0.09,0.11,0])
Q_upper_err_tom2014_hi_z = np.array([0,0,0.08,0.08,0.10,0.10,0.11,0.08,0.08,0.09,0.09,0.10,0.11,0.14,0])
# in unit of log(phi) - i.e. 10^ - (SF_smf_tom2014)
Q_upper_err_tom2014_hi_z = np.abs(10**(-1 * (Q_smf_tom2014_hi_z - Q_upper_err_tom2014_hi_z)) - 10**(-1 * (Q_smf_tom2014_hi_z)))
Q_lower_err_tom2014_hi_z = np.abs(10**(-1 * (Q_smf_tom2014_hi_z + Q_lower_err_tom2014_hi_z)) - 10**(-1 * (Q_smf_tom2014_hi_z)))
Q_smf_tom2014_hi_z = 10**(-1 * Q_smf_tom2014_hi_z)
Q_err_tom2014_hi_z = np.array([Q_lower_err_tom2014_hi_z, Q_upper_err_tom2014_hi_z])
# TOTAL population
# lo z
total_smf_tom2014_lo_z = np.array([1.37,1.53,1.71,1.86,2.03,2.01,2.10,2.17,2.24,2.31,2.41,2.53,2.91,3.46,0])
total_lower_err_tom2014_lo_z = np.array([0.06,0.06,0.07,0.07,0.08,0.07,0.07,0.08,0.08,0.08,0.08,0.09,0.11,0.14,0])
total_upper_err_tom2014_lo_z = np.array([0.07,0.07,0.08,0.08,0.09,0.08,0.09,0.10,0.10,0.09,0.10,0.11,0.15,0.18,0])
# in unit of log(phi) - i.e. 10^ - (SF_smf_tom2014)
total_upper_err_tom2014_lo_z = np.abs(10**(-1 * (total_smf_tom2014_lo_z - total_upper_err_tom2014_lo_z)) - 10**(-1 * (total_smf_tom2014_lo_z)))
total_lower_err_tom2014_lo_z = np.abs(10**(-1 * (total_smf_tom2014_lo_z + total_lower_err_tom2014_lo_z)) - 10**(-1 * (total_smf_tom2014_lo_z)))
total_smf_tom2014_lo_z = 10**(-1 * total_smf_tom2014_lo_z)
total_err_tom2014_lo_z = np.array([total_lower_err_tom2014_lo_z, total_upper_err_tom2014_lo_z])
# hi z
total_smf_tom2014_hi_z = np.array([0,1.53,1.60,1.76,1.86,2.00,2.12,2.21,2.25,2.35,2.45,2.55,2.82,3.32,0])
total_lower_err_tom2014_hi_z = np.array([0,0.06,0.05,0.06,0.06,0.06,0.07,0.06,0.06,0.07,0.07,0.08,0.09,0.10,0])
total_upper_err_tom2014_hi_z = np.array([0,0.07,0.06,0.06,0.07,0.07,0.08,0.07,0.08,0.08,0.09,0.09,0.11,0.13,0])
# in unit of log(phi) - i.e. 10^ - (SF_smf_tom2014)
total_upper_err_tom2014_hi_z = np.abs(10**(-1 * (total_smf_tom2014_hi_z - total_upper_err_tom2014_hi_z)) - 10**(-1 * (total_smf_tom2014_hi_z)))
total_lower_err_tom2014_hi_z = np.abs(10**(-1 * (total_smf_tom2014_hi_z + total_lower_err_tom2014_hi_z)) - 10**(-1 * (total_smf_tom2014_hi_z)))
total_smf_tom2014_hi_z = 10**(-1 * total_smf_tom2014_hi_z)
total_err_tom2014_hi_z = np.array([total_lower_err_tom2014_hi_z, total_upper_err_tom2014_hi_z])
#
#
#
## from Baldry et al. 2012: (total pop only, no split between SF/Q)
#
# x-values (mass bin midpoints)
x_baldry2012_data = np.arange(7.9,12.1,0.2)
# TOTAL population
total_smf_baldry2012 = np.array([27.3,28.3,23.5,19.2,18.0,14.3,10.2,9.59,7.42,6.21,5.71,5.51,5.48,5.12,3.55,2.41,1.27,0.338,0.042,0.021,0.042])
total_err_baldry2012 = np.array([4.2,2.8,3.0,1.2,2.6,1.7,0.6,0.55,0.41,0.37,0.35,0.34,0.34,0.33,0.27,0.23,0.16,0.085,0.030,0.21,0.030])
# units of 10^-3
total_smf_baldry2012 = 1e-3*total_smf_baldry2012
total_err_baldry2012 = 1e-3*total_err_baldry2012
#
#
#
## from van der Burg 2018:
#
# x-values (mass bins midpoints)
x_vdb2018_data = np.arange(9.55,11.85,0.1)
# SF population
SF_smf_vdb2018 = np.array([497.4,545.4,423.0,417.6,348.0,304.1,283.8,233.8,213.6,143.9,138.5,119.6,71.63,43.93,34.14,21.63,7.43,2.03,2.03,0.68,0,0,0])
SF_upper_err_vdb2018 = np.array([17.6,16.9,19.6,18.2,14.2,14.9,12.8,12.2,14.2,10.1,8.8,8.1,6.76,5.41,6.08,3.38,1.35,0.68,0.68,0.68,0,0,0])
SF_lower_err_vdb2018 = np.array([13.5,18.9,17.6,19.6,14.9,10.8,14.2,12.8,10.1,7.4,10.1,10.8,6.76,4.05,4.736,2.70,2.03,1.35,0.68,0.68,0,0,0])
SF_err_vdb2018 = np.array([SF_lower_err_vdb2018, SF_upper_err_vdb2018])
# Q population
Q_smf_vdb2018 = np.array([79.74,93.94,83.8,108.1,109.5,110.8,150.0,140.6,148.37,138.5,150.0,152.1,134.5,123.0,103.4,77.04,59.47,31.76,22.3,14.87,6.76,4.73,0])
Q_upper_err_vdb2018 = np.array([7.43,6.76,7.43,8.1,8.8,8.1,8.1,9.5,10.1,8.8,8.8,12.2,7.4,10.1,9.5,6.76,6.08,6.08,4.73,2.70,2.03,2.03,0])
Q_lower_err_vdb2018 = np.array([5.41,8.11,7.43,8.8,10.1,7.4,8.8,10.8,9.5,9.5,9.5,12.2,6.8,10.1,8.1,6.76,7.43,3.38,4.05,3.38,1.35,2.03,0])
Q_err_vdb2018 = np.array([Q_lower_err_vdb2018, Q_upper_err_vdb2018])
# TOTAL population
total_smf_vdb2018 = np.array([577.8,639.3,506.2,522.4,458.9,417.6,435.2,374.4,364.3,281.8,287.9,270.3,206.8,166.9,138.5,99.34,65.55,33.79,24.33,14.87,6.76,4.73,0])
total_upper_err_vdb2018 = np.array([16.9,18.9,21.6,25.0,14.2,10.1,14.9,16.9,18.9,14.9,11.5,16.2,8.8,12.8,10.8,6.76,6.76,6.08,4.73,3.38,2.03,2.03,0])
total_lower_err_vdb2018 = np.array([14.2,20.3,18.2,20.3,16.9,14.9,16.9,16.2,16.2,11.5,14.2,14.2,8.8,9.5,7.4,8.11,5.41,4.01,4.01,3.38,1.35,2.03,0])
total_err_vdb2018 = np.array([total_lower_err_vdb2018, total_upper_err_vdb2018])
# apply units of 10^-5
SF_smf_vdb2018 = 1e-5*SF_smf_vdb2018
SF_err_vdb2018 = 1e-5*SF_err_vdb2018
Q_smf_vdb2018 = 1e-5*Q_smf_vdb2018
Q_err_vdb2018 = 1e-5*Q_err_vdb2018
total_smf_vdb2018 = 1e-5*total_smf_vdb2018
total_err_vdb2018 = 1e-5*total_err_vdb2018
#
#
#
## SCETION (1.2): Visualization
#
# Plot: compare environments by population (i.e. plot Cluster vs Field for Total, SF, & Q)
if plot_flag_1 == 1:
    SMF = plt.figure()
    gs = gridspec.GridSpec(1,3, wspace=0, hspace=0, width_ratios=[1,1,1])   #make a tiled-plot
    # Total population
    ax0 = plt.subplot(gs[0])
    ax0.errorbar(x_vdb2018_data,total_smf_vdb2018,yerr=total_err_vdb2018,label='vdB 2018: 0.5 < z < 0.7',marker='o',color='magenta',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    ax0.errorbar(x_baldry2012_data,total_smf_baldry2012,yerr=total_err_baldry2012,label='baldry 2012: z < 0.06',marker='s',color='cyan',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    ax0.errorbar(x_tom2014_data,total_smf_tom2014_lo_z,yerr=total_err_tom2014_lo_z,label='tomczak 2014: 0.2 < z < 0.5',marker='<',color='orangered',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    ax0.errorbar(x_tom2014_data,total_smf_tom2014_hi_z,yerr=total_err_tom2014_hi_z,label='tomczak 2014: 0.5 < z < 0.75',marker='>',color='maroon',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    #Plot Schechter fits:
    ax0.plot(x_plot,total_model_field_mcmc_plot_double, '-k', linewidth=1.0)
    ax0.set_xscale('linear')
    ax0.set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    ax0.set_xlim(7.8,12.5)
    ax0.set_yscale('log')
    ax0.set_ylim(1e-6,0.2)
    ax0.minorticks_on()
    ax0.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=False,labelleft=True,labelbottom=True,labelsize=18)
    ax0.yaxis.set_label_position("left")
    ax0.set_ylabel('${\phi}$ [$Mpc^{-3}$ $dex^{-1}$]',fontsize=20)
    ax0.set_title('Total')
    ax0.legend(scatterpoints=1,loc='upper right', frameon=False, fontsize = 'medium')
    # ax0.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')
    #
    # SF population
    ax1 = plt.subplot(gs[1])
    ax1.errorbar(x_vdb2018_data,SF_smf_vdb2018,yerr=SF_err_vdb2018,label='vdB 2018',marker='o',color='magenta',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    ax1.errorbar(x_tom2014_data,SF_smf_tom2014_lo_z,yerr=SF_err_tom2014_lo_z,label='tomczak 2014: 0.2 < z < 0.5',marker='<',color='orangered',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    ax1.errorbar(x_tom2014_data,SF_smf_tom2014_hi_z,yerr=SF_err_tom2014_hi_z,label='tomczak 2014: 0.5 < z < 0.75',marker='>',color='maroon',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    #
    #Plot Schechter fits:
    ax1.plot(x_plot,SF_model_field_mcmc_plot, '-k', linewidth=1.0)
    ax1.set_xscale('linear')
    ax1.set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    ax1.set_xlim(7.8,12.5)
    ax1.set_yscale('log')
    ax1.set_ylim(1e-6,0.2)
    ax1.minorticks_on()
    ax1.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=False,labelleft=False,labelbottom=True,labelsize=18)
    # ax0.yaxis.set_label_position("left")
    # ax0.ylabel('???')
    ax1.set_title('Star-Forming')
    # ax0.legend(scatterpoints=1,loc='lower left', frameon=False, fontsize = 'x-small')
    # ax1.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')#
    #
    # Q population
    ax2 = plt.subplot(gs[2])
    ax2.errorbar(x_vdb2018_data,Q_smf_vdb2018,yerr=Q_err_vdb2018,label='vdB 2018',marker='o',color='magenta',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    ax2.errorbar(x_tom2014_data,Q_smf_tom2014_lo_z,yerr=Q_err_tom2014_lo_z,label='tomczak 2014: 0.2 < z < 0.5',marker = '<', color = 'orangered',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    ax2.errorbar(x_tom2014_data,Q_smf_tom2014_hi_z,yerr=Q_err_tom2014_hi_z,label='tomczak 2014: 0.5 < z < 0.75',marker = '>', color = 'maroon',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    #Plot Schechter fits:
    ax2.plot(x_plot,Q_model_field_mcmc_plot_double, '-k', linewidth=1.0)
    ax2.set_xscale('linear')
    ax2.set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    ax2.set_xlim(7.8,12.5)
    ax2.set_yscale('log')
    ax2.set_ylim(1e-6,0.2)
    ax2.minorticks_on()
    ax2.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright=True, labelleft=False,labelbottom=True,labelsize=18)
    ax2.yaxis.set_label_position("right")
    ax2.set_ylabel('${\phi}$ [$Mpc^{-3}$ $dex^{-1}$]',fontsize=20)
    ax2.set_title('Quiescent')
    # ax0.legend(scatterpoints=1,loc='lower left', frameon=False, fontsize = 'x-small')
    # ax2.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')#
    #
    plt.show()
#
#
#
## SECTION (2): plotting various Schechter functions from the literature, including 1-sigma confidence regions
#  recall: vector for plotting is called "x_plot" has range [7.8,12.5]
#
## SECTION (2.1): DATA ENTRY
#
## Muzzin et al. 2013
#
# x vectors for plotting
x_muzzin2013_lo_z = np.arange(8.37,12.01,0.01)
x_muzzin2013_hi_z = np.arange(8.92,12.01,0.01)
# SF population
# lo z
SF_M_star_muzzin_lo_z, SF_phi1_muzzin_lo_z, SF_alpha1_muzzin_lo_z = 10.81,11.35e-4,-1.34
SF_M_star_upper_err_muzzin_lo_z, SF_phi1_upper_err_muzzin_lo_z, SF_alpha1_upper_err_muzzin_lo_z = 0.03,1.60e-4,0.01
SF_M_star_lower_err_muzzin_lo_z, SF_phi1_lower_err_muzzin_lo_z, SF_alpha1_lower_err_muzzin_lo_z = 0.03,1.51e-4,0.01
# hi z
SF_M_star_muzzin_hi_z, SF_phi1_muzzin_hi_z, SF_alpha1_muzzin_hi_z = 10.78,12.71e-4,-1.26
SF_M_star_upper_err_muzzin_hi_z, SF_phi1_upper_err_muzzin_hi_z, SF_alpha1_upper_err_muzzin_hi_z = 0.01,1.07e-4,0.02
SF_M_star_lower_err_muzzin_hi_z, SF_phi1_lower_err_muzzin_hi_z, SF_alpha1_lower_err_muzzin_hi_z = 0.02,0.91e-4,0.01
# Q population
# lo z
Q_M_star_muzzin_lo_z, Q_phi1_muzzin_lo_z, Q_alpha1_muzzin_lo_z, Q_phi2_muzzin_lo_z, Q_alpha2_muzzin_lo_z = 10.92,19.68e-4,-0.38,0.58e-4,-1.52
Q_M_star_upper_err_muzzin_lo_z, Q_phi1_upper_err_muzzin_lo_z, Q_alpha1_upper_err_muzzin_lo_z, Q_phi2_upper_err_muzzin_lo_z, Q_alpha2_upper_err_muzzin_lo_z = 0.06,2.64e-4,0.06,0.27e-4,0.06
Q_M_star_lower_err_muzzin_lo_z, Q_phi1_lower_err_muzzin_lo_z, Q_alpha1_lower_err_muzzin_lo_z, Q_phi2_lower_err_muzzin_lo_z, Q_alpha2_lower_err_muzzin_lo_z = 0.02,2.96e-4,0.12,0.32e-4,0.16
# hi z
Q_M_star_muzzin_hi_z, Q_phi1_muzzin_hi_z, Q_alpha1_muzzin_hi_z, Q_phi2_muzzin_hi_z, Q_alpha2_muzzin_hi_z =  10.84,14.55e-4,-0.36,0.005e-4,-2.32
Q_M_star_upper_err_muzzin_hi_z, Q_phi1_upper_err_muzzin_hi_z, Q_alpha1_upper_err_muzzin_hi_z, Q_phi2_upper_err_muzzin_hi_z, Q_alpha2_upper_err_muzzin_hi_z = 0.03,1.21e-4,0.06,0.021e-4,0.41
Q_M_star_lower_err_muzzin_hi_z, Q_phi1_lower_err_muzzin_hi_z, Q_alpha1_lower_err_muzzin_hi_z, Q_phi2_lower_err_muzzin_hi_z, Q_alpha2_lower_err_muzzin_hi_z = 0.03,1.09e-4,0.04,0.004e-4,0.38
# TOTAL population
# lo z
T_M_star_muzzin_lo_z, T_phi1_muzzin_lo_z, T_alpha1_muzzin_lo_z, T_phi2_muzzin_lo_z, T_alpha2_muzzin_lo_z = 10.97,16.27e-4,-0.53,9.47e-4,-1.37
T_M_star_upper_err_muzzin_lo_z, T_phi1_upper_err_muzzin_lo_z, T_alpha1_upper_err_muzzin_lo_z, T_phi2_upper_err_muzzin_lo_z, T_alpha2_upper_err_muzzin_lo_z = 0.06,3.88e-4,0.16,2.32e-4,0.01
T_M_star_lower_err_muzzin_lo_z, T_phi1_lower_err_muzzin_lo_z, T_alpha1_lower_err_muzzin_lo_z, T_phi2_lower_err_muzzin_lo_z, T_alpha2_lower_err_muzzin_lo_z = 0.06,2.41e-4,0.28,3.83e-4,0.06
# hi z
T_M_star_muzzin_hi_z, T_phi1_muzzin_hi_z, T_alpha1_muzzin_hi_z = 11.00,16.25e-4,-1.17
T_M_star_upper_err_muzzin_hi_z, T_phi1_upper_err_muzzin_hi_z, T_alpha1_upper_err_muzzin_hi_z = 0.02,1.17e-4,0.01
T_M_star_lower_err_muzzin_hi_z, T_phi1_lower_err_muzzin_hi_z, T_alpha1_lower_err_muzzin_hi_z = 0.01,1.28e-4,0.01
#
#
## Tomczak et al. 2014
#
# x vector for plotting
x_tomczak2014_lo_z = np.arange(8.0,11.51,0.01)
x_tomczak2014_hi_z = np.arange(8.25,11.51,0.01)
x_tomczak2014_lo_z_Q = np.arange(8.25,11.51,0.01)
x_tomczak2014_hi_z_Q = np.arange(8.50,11.51,0.01)
# SF population
# lo z
SF_M_star_tomczak_lo_z, SF_phi1_tomczak_lo_z, SF_alpha1_tomczak_lo_z, SF_phi2_tomczak_lo_z, SF_alpha2_tomczak_lo_z = 10.59,10**(-2.67),-1.08,10**(-4.46),-2.00
SF_M_star_upper_err_tomczak_lo_z, SF_phi1_upper_err_tomczak_lo_z, SF_alpha1_upper_err_tomczak_lo_z, SF_phi2_upper_err_tomczak_lo_z, SF_alpha2_upper_err_tomczak_lo_z = 0.09,10**(-2.67+0.11),0.23,10**(-4.46+0.63),0.49
SF_M_star_lower_err_tomczak_lo_z, SF_phi1_lower_err_tomczak_lo_z, SF_alpha1_lower_err_tomczak_lo_z, SF_phi2_lower_err_tomczak_lo_z, SF_alpha2_lower_err_tomczak_lo_z = 0.09,10**(-2.67-0.11),0.23,10**(-4.46-0.63),0.49
# hi z
SF_M_star_tomczak_hi_z, SF_phi1_tomczak_hi_z, SF_alpha1_tomczak_hi_z, SF_phi2_tomczak_hi_z, SF_alpha2_tomczak_hi_z = 10.65,10**(-2.97),-0.97,10**(-3.34),-1.58
SF_M_star_upper_err_tomczak_hi_z, SF_phi1_upper_err_tomczak_hi_z, SF_alpha1_upper_err_tomczak_hi_z, SF_phi2_upper_err_tomczak_hi_z, SF_alpha2_upper_err_tomczak_hi_z = 0.23,10**(-2.97+0.28),1.32,10**(-3.34+0.67),0.54
SF_M_star_lower_err_tomczak_hi_z, SF_phi1_lower_err_tomczak_hi_z, SF_alpha1_lower_err_tomczak_hi_z, SF_phi2_lower_err_tomczak_hi_z, SF_alpha2_lower_err_tomczak_hi_z = 0.23,10**(-2.97-0.28),1.32,10**(-3.34-0.67),0.54
# Q population
# lo z
Q_M_star_tomczak_lo_z, Q_phi1_tomczak_lo_z, Q_alpha1_tomczak_lo_z, Q_phi2_tomczak_lo_z, Q_alpha2_tomczak_lo_z = 10.75,10**(-2.76),-0.47,10**(-5.21),-1.97
Q_M_star_upper_err_tomczak_lo_z, Q_phi1_upper_err_tomczak_lo_z, Q_alpha1_upper_err_tomczak_lo_z, Q_phi2_upper_err_tomczak_lo_z, Q_alpha2_upper_err_tomczak_lo_z = 0.10,10**(-2.76+0.09),0.20,10**(-5.21+0.48),0.34
Q_M_star_lower_err_tomczak_lo_z, Q_phi1_lower_err_tomczak_lo_z, Q_alpha1_lower_err_tomczak_lo_z, Q_phi2_lower_err_tomczak_lo_z, Q_alpha2_lower_err_tomczak_lo_z = 0.10,10**(-2.76-0.09),0.20,10**(-5.21-0.48),0.34
# hi z
Q_M_star_tomczak_hi_z, Q_phi1_tomczak_hi_z, Q_alpha1_tomczak_hi_z, Q_phi2_tomczak_hi_z, Q_alpha2_tomczak_hi_z = 10.68,10**(-2.67),-0.10,10**(-4.29),-1.69
Q_M_star_upper_err_tomczak_hi_z, Q_phi1_upper_err_tomczak_hi_z, Q_alpha1_upper_err_tomczak_hi_z, Q_phi2_upper_err_tomczak_hi_z, Q_alpha2_upper_err_tomczak_hi_z = 0.07,10**(-2.67+0.05),0.27,10**(-4.29+0.33),0.24
Q_M_star_lower_err_tomczak_hi_z, Q_phi1_lower_err_tomczak_hi_z, Q_alpha1_lower_err_tomczak_hi_z, Q_phi2_lower_err_tomczak_hi_z, Q_alpha2_lower_err_tomczak_hi_z = 0.07,10**(-2.67-0.05),0.27,10**(-4.29-0.33),0.24
# TOTAL population
# lo z
T_M_star_tomczak_lo_z, T_phi1_tomczak_lo_z, T_alpha1_tomczak_lo_z, T_phi2_tomczak_lo_z, T_alpha2_tomczak_lo_z = 10.78,10**(-2.54),-0.98,10**(-4.29),-1.90
T_M_star_upper_err_tomczak_lo_z, T_phi1_upper_err_tomczak_lo_z, T_alpha1_upper_err_tomczak_lo_z, T_phi2_upper_err_tomczak_lo_z, T_alpha2_upper_err_tomczak_lo_z = 0.11,10**(-2.54+0.12),0.24,10**(-4.29+0.55),0.36
T_M_star_lower_err_tomczak_lo_z, T_phi1_lower_err_tomczak_lo_z, T_alpha1_lower_err_tomczak_lo_z, T_phi2_lower_err_tomczak_lo_z, T_alpha2_lower_err_tomczak_lo_z = 0.11,10**(-2.54-0.12),0.24,10**(-4.29-0.55),0.36
# hi z
T_M_star_tomczak_hi_z, T_phi1_tomczak_hi_z, T_alpha1_tomczak_hi_z, T_phi2_tomczak_hi_z, T_alpha2_tomczak_hi_z = 10.70,10**(-2.55),-0.39,10**(-3.15),-1.53
T_M_star_upper_err_tomczak_hi_z, T_phi1_upper_err_tomczak_hi_z, T_alpha1_upper_err_tomczak_hi_z, T_phi2_upper_err_tomczak_hi_z, T_alpha2_upper_err_tomczak_hi_z = 0.10,10**(-2.55+0.09),0.50,10**(-3.15+0.23),0.12
T_M_star_lower_err_tomczak_hi_z, T_phi1_lower_err_tomczak_hi_z, T_alpha1_lower_err_tomczak_hi_z, T_phi2_lower_err_tomczak_hi_z, T_alpha2_lower_err_tomczak_hi_z = 0.10,10**(-2.55-0.09),0.50,10**(-3.15-0.23),0.12
#
## van der Burg 2018
#
# x vector for plotting
x_vdb2018 = np.arange(9.55,11.76,0.01)
# SF population
SF_M_star_vdb2018, SF_phi1_vdb2018, SF_alpha1_vdb2018 = 10.69,238.72e-5,-1.33
SF_M_star_upper_err_vdb2018, SF_phi1_upper_err_vdb2018, SF_alpha1_upper_err_vdb2018 = 0.03,3.20e-5,0.03
SF_M_star_lower_err_vdb2018, SF_phi1_lower_err_vdb2018, SF_alpha1_lower_err_vdb2018 = 0.03,3.20e-5,0.03
# Q populaQion
Q_M_star_vdb2018, Q_phi1_vdb2018, Q_alpha1_vdb2018 = 10.86,313.85e-5,-0.55
Q_M_star_upper_err_vdb2018, Q_phi1_upper_err_vdb2018, Q_alpha1_upper_err_vdb2018 = 0.02,5.83e-5,0.03
Q_M_star_lower_err_vdb2018, Q_phi1_lower_err_vdb2018, Q_alpha1_lower_err_vdb2018 = 0.02,5.83e-5,0.03
# T population
T_M_star_vdb2018, T_phi1_vdb2018, T_alpha1_vdb2018 = 10.98,320.59e-5,-1.20
T_M_star_upper_err_vdb2018, T_phi1_upper_err_vdb2018, T_alpha1_upper_err_vdb2018 = 0.02,3.49e-5,0.02
T_M_star_lower_err_vdb2018, T_phi1_lower_err_vdb2018, T_alpha1_lower_err_vdb2018 = 0.02,3.49e-5,0.02
#
## baldry+2012
#
# x vector for plotting
x_baldry2012 = np.arange(8,11.61,0.01)
x_baldry2012_Q = np.arange(8.4,11.61,0.01)
# SF population
SF_M_star_baldry2012, SF_phi1_baldry2012, SF_alpha1_baldry2012 = 10.72,0.71e-3,-1.45
SF_M_star_upper_err_baldry2012, SF_phi1_upper_err_baldry2012, SF_alpha1_upper_err_baldry2012 = 0.,0.,0.
SF_M_star_lower_err_baldry2012, SF_phi1_lower_err_baldry2012, SF_alpha1_lower_err_baldry2012 = 0.,0.,0.
# Q population
Q_M_star_baldry2012, Q_phi1_baldry2012, Q_alpha1_baldry2012, Q_phi2_baldry2012, Q_alpha2_baldry2012 = 10.72,3.25e-3,-0.45,0.08e-3,-1.45
Q_M_star_upper_err_baldry2012, Q_phi1_upper_err_baldry2012, Q_alpha1_upper_err_baldry2012, Q_phi2_upper_err_baldry2012, Q_alpha2_upper_err_baldry2012 = 0.,0.,0.,0,0
Q_M_star_lower_err_baldry2012, Q_phi1_lower_err_baldry2012, Q_alpha1_lower_err_baldry2012, Q_phi2_lower_err_baldry2012, Q_alpha2_lower_err_baldry2012 = 0.,0.,0.,0,0
# T population
T_M_star_baldry2012, T_phi1_baldry2012, T_alpha1_baldry2012, T_phi2_baldry2012, T_alpha2_baldry2012 = 10.66,3.96e-3,-0.35,0.79e-3,-1.47
T_M_star_upper_err_baldry2012, T_phi1_upper_err_baldry2012, T_alpha1_upper_err_baldry2012, T_phi2_upper_err_baldry2012, T_alpha2_upper_err_baldry2012 = 0.05,0.34e-3,0.18,0.23e-3,0.05
T_M_star_lower_err_baldry2012, T_phi1_lower_err_baldry2012, T_alpha1_lower_err_baldry2012, T_phi2_lower_err_baldry2012, T_alpha2_lower_err_baldry2012 = 0.05,0.34e-3,0.18,0.23e-3,0.05


#
#
#
## SECTION (2.2): build SCHECHTER MODELS
#
## construct upper/lower uncertainty limits for MY FIT (Sarrouh & Muzzin 2021)
#
SF_model_field_upper =  np.log(10) * np.exp(-10**(x_plot-(SFM_star_mcmc_field+SFM_star_field_mcmc_upper_err))) * ((SFphi_mcmc_field+SFphi_field_mcmc_upper_err) * (10**(x_plot-(SFM_star_mcmc_field+SFM_star_field_mcmc_upper_err)))**(1+(SFalpha_mcmc_field-SFalpha_field_mcmc_upper_err)))
SF_model_field_lower =  np.log(10) * np.exp(-10**(x_plot-(SFM_star_mcmc_field-SFM_star_field_mcmc_lower_err))) * ((SFphi_mcmc_field-SFphi_field_mcmc_lower_err) * (10**(x_plot-(SFM_star_mcmc_field-SFM_star_field_mcmc_lower_err)))**(1+(SFalpha_mcmc_field+SFalpha_field_mcmc_lower_err)))
Q_model_field_upper =  np.log(10) * np.exp(-10**(x_plot-(QM_star_mcmc_field+QM_star_field_mcmc_upper_err))) * ( ((Qphi1_mcmc+Qphi1_field_mcmc_upper_err) * (10**(x_plot-(QM_star_mcmc_field+QM_star_field_mcmc_upper_err)))**(1+(Qalpha1_mcmc-Qalpha1_field_mcmc_upper_err)))  + ((Qphi2_mcmc+Qphi2_field_mcmc_upper_err)*(10**(x_plot-(QM_star_mcmc_field+QM_star_field_mcmc_upper_err)))**(1+(Qalpha2_mcmc-Qalpha2_field_mcmc_upper_err))) )
Q_model_field_lower =  np.log(10) * np.exp(-10**(x_plot-(QM_star_mcmc_field-QM_star_field_mcmc_lower_err))) * ( ((Qphi1_mcmc-Qphi1_field_mcmc_lower_err) * (10**(x_plot-(QM_star_mcmc_field-QM_star_field_mcmc_lower_err)))**(1+(Qalpha1_mcmc+Qalpha1_field_mcmc_lower_err)))  + ((Qphi2_mcmc-Qphi2_field_mcmc_lower_err)*(10**(x_plot-(QM_star_mcmc_field-QM_star_field_mcmc_lower_err)))**(1+(Qalpha2_mcmc+Qalpha2_field_mcmc_lower_err))) )
total_model_field_upper =  np.log(10) * np.exp(-10**(x_plot-(TM_star_mcmc_field+TM_star_field_mcmc_upper_err))) * ( ((Tphi1_mcmc+Tphi1_field_mcmc_upper_err) * (10**(x_plot-(TM_star_mcmc_field+TM_star_field_mcmc_upper_err)))**(1+(Talpha1_mcmc-Talpha1_field_mcmc_upper_err)))  + ((Tphi2_mcmc+Tphi2_field_mcmc_upper_err)*(10**(x_plot-(TM_star_mcmc_field+TM_star_field_mcmc_upper_err)))**(1+(Talpha2_mcmc-Talpha2_field_mcmc_upper_err))) )
total_model_field_lower =  np.log(10) * np.exp(-10**(x_plot-(TM_star_mcmc_field-TM_star_field_mcmc_lower_err))) * ( ((Tphi1_mcmc-Tphi1_field_mcmc_lower_err) * (10**(x_plot-(TM_star_mcmc_field-TM_star_field_mcmc_lower_err)))**(1+(Talpha1_mcmc+Talpha1_field_mcmc_lower_err)))  + ((Tphi2_mcmc-Tphi2_field_mcmc_lower_err)*(10**(x_plot-(TM_star_mcmc_field-TM_star_field_mcmc_lower_err)))**(1+(Talpha2_mcmc+Talpha2_field_mcmc_lower_err))) )
#
#
## build MODELS FROM LITERATURE for SMF, SMF_upper_error_bound & SMF_lower_error_bound
#
# muzzin, lo z
# SMFs
SF_model_muzzin_lo_z = np.log(10) * SF_phi1_muzzin_lo_z * (10**((x_muzzin2013_lo_z-SF_M_star_muzzin_lo_z)*(1+SF_alpha1_muzzin_lo_z))) * np.exp(-10**(x_muzzin2013_lo_z-SF_M_star_muzzin_lo_z))#
Q_model_muzzin_lo_z = np.log(10) * np.exp(-10**(x_muzzin2013_lo_z-Q_M_star_muzzin_lo_z)) * ( (Q_phi1_muzzin_lo_z*(10**(x_muzzin2013_lo_z-Q_M_star_muzzin_lo_z))**(1+Q_alpha1_muzzin_lo_z))  + (Q_phi2_muzzin_lo_z*(10**(x_muzzin2013_lo_z-Q_M_star_muzzin_lo_z))**(1+Q_alpha2_muzzin_lo_z)) )
T_model_muzzin_lo_z = np.log(10) * np.exp(-10**(x_muzzin2013_lo_z-T_M_star_muzzin_lo_z)) * ( (T_phi1_muzzin_lo_z*(10**(x_muzzin2013_lo_z-T_M_star_muzzin_lo_z))**(1+T_alpha1_muzzin_lo_z))  + (T_phi2_muzzin_lo_z*(10**(x_muzzin2013_lo_z-T_M_star_muzzin_lo_z))**(1+T_alpha2_muzzin_lo_z)) )
# upper_error_bound
SF_model_muzzin_upper_lo_z = np.log(10) * (SF_phi1_muzzin_lo_z+SF_phi1_upper_err_muzzin_lo_z) * (10**((x_muzzin2013_lo_z-(SF_M_star_muzzin_lo_z+SF_M_star_upper_err_muzzin_lo_z))*(1+(SF_alpha1_muzzin_lo_z-SF_alpha1_upper_err_muzzin_lo_z)))) * np.exp(-10**(x_muzzin2013_lo_z-(SF_M_star_muzzin_lo_z+SF_M_star_upper_err_muzzin_lo_z)))#
Q_model_muzzin_upper_lo_z = np.log(10) * np.exp(-10**(x_muzzin2013_lo_z-(Q_M_star_muzzin_lo_z+Q_M_star_upper_err_muzzin_lo_z))) * ( ((Q_phi1_muzzin_lo_z+Q_phi1_upper_err_muzzin_lo_z)*(10**(x_muzzin2013_lo_z-(Q_M_star_muzzin_lo_z-Q_M_star_upper_err_muzzin_lo_z)))**(1+(Q_alpha1_muzzin_lo_z-Q_alpha1_upper_err_muzzin_lo_z)))  + ((Q_phi2_muzzin_lo_z+Q_phi2_upper_err_muzzin_lo_z)*(10**(x_muzzin2013_lo_z-(Q_M_star_muzzin_lo_z+Q_M_star_upper_err_muzzin_lo_z)))**(1+(Q_alpha2_muzzin_lo_z-Q_alpha2_upper_err_muzzin_lo_z))) )
T_model_muzzin_upper_lo_z = np.log(10) * np.exp(-10**(x_muzzin2013_lo_z-(T_M_star_muzzin_lo_z+T_M_star_upper_err_muzzin_lo_z))) * ( ((T_phi1_muzzin_lo_z+T_phi1_upper_err_muzzin_lo_z)*(10**(x_muzzin2013_lo_z-(T_M_star_muzzin_lo_z-T_M_star_upper_err_muzzin_lo_z)))**(1+(T_alpha1_muzzin_lo_z-T_alpha1_upper_err_muzzin_lo_z)))  + ((T_phi2_muzzin_lo_z+T_phi2_upper_err_muzzin_lo_z)*(10**(x_muzzin2013_lo_z-(T_M_star_muzzin_lo_z+T_M_star_upper_err_muzzin_lo_z)))**(1+(T_alpha2_muzzin_lo_z-T_alpha2_upper_err_muzzin_lo_z))) )
# lower_error_bound
SF_model_muzzin_lower_lo_z = np.log(10) * (SF_phi1_muzzin_lo_z-SF_phi1_lower_err_muzzin_lo_z) * (10**((x_muzzin2013_lo_z-(SF_M_star_muzzin_lo_z-SF_M_star_lower_err_muzzin_lo_z))*(1+(SF_alpha1_muzzin_lo_z+SF_alpha1_lower_err_muzzin_lo_z)))) * np.exp(-10**(x_muzzin2013_lo_z-(SF_M_star_muzzin_lo_z-SF_M_star_lower_err_muzzin_lo_z)))#
Q_model_muzzin_lower_lo_z = np.log(10) * np.exp(-10**(x_muzzin2013_lo_z-(Q_M_star_muzzin_lo_z-Q_M_star_lower_err_muzzin_lo_z))) * ( ((Q_phi1_muzzin_lo_z-Q_phi1_lower_err_muzzin_lo_z)*(10**(x_muzzin2013_lo_z-(Q_M_star_muzzin_lo_z+Q_M_star_lower_err_muzzin_lo_z)))**(1+(Q_alpha1_muzzin_lo_z+Q_alpha1_lower_err_muzzin_lo_z)))  + ((Q_phi2_muzzin_lo_z-Q_phi2_lower_err_muzzin_lo_z)*(10**(x_muzzin2013_lo_z-(Q_M_star_muzzin_lo_z-Q_M_star_lower_err_muzzin_lo_z)))**(1+(Q_alpha2_muzzin_lo_z+Q_alpha2_lower_err_muzzin_lo_z))) )
T_model_muzzin_lower_lo_z = np.log(10) * np.exp(-10**(x_muzzin2013_lo_z-(T_M_star_muzzin_lo_z-T_M_star_lower_err_muzzin_lo_z))) * ( ((T_phi1_muzzin_lo_z-T_phi1_lower_err_muzzin_lo_z)*(10**(x_muzzin2013_lo_z-(T_M_star_muzzin_lo_z+T_M_star_lower_err_muzzin_lo_z)))**(1+(T_alpha1_muzzin_lo_z+T_alpha1_lower_err_muzzin_lo_z)))  + ((T_phi2_muzzin_lo_z-T_phi2_lower_err_muzzin_lo_z)*(10**(x_muzzin2013_lo_z-(T_M_star_muzzin_lo_z-T_M_star_lower_err_muzzin_lo_z)))**(1+(T_alpha2_muzzin_lo_z+T_alpha2_lower_err_muzzin_lo_z))) )
# muzzin, hi z
# SMFs
SF_model_muzzin_hi_z = np.log(10) * SF_phi1_muzzin_hi_z * (10**((x_muzzin2013_hi_z-SF_M_star_muzzin_hi_z)*(1+SF_alpha1_muzzin_hi_z))) * np.exp(-10**(x_muzzin2013_hi_z-SF_M_star_muzzin_hi_z))#
Q_model_muzzin_hi_z = np.log(10) * np.exp(-10**(x_muzzin2013_hi_z-Q_M_star_muzzin_hi_z)) * ( (Q_phi1_muzzin_hi_z*(10**(x_muzzin2013_hi_z-Q_M_star_muzzin_hi_z))**(1+Q_alpha1_muzzin_hi_z))  + (Q_phi2_muzzin_hi_z*(10**(x_muzzin2013_hi_z-Q_M_star_muzzin_hi_z))**(1+Q_alpha2_muzzin_hi_z)) )
T_model_muzzin_hi_z = np.log(10) * T_phi1_muzzin_hi_z * (10**((x_muzzin2013_hi_z-T_M_star_muzzin_hi_z)*(1+T_alpha1_muzzin_hi_z))) * np.exp(-10**(x_muzzin2013_hi_z-T_M_star_muzzin_hi_z))#
# upper_error_bound
SF_model_muzzin_upper_hi_z = np.log(10) * (SF_phi1_muzzin_hi_z+SF_phi1_upper_err_muzzin_hi_z) * (10**((x_muzzin2013_hi_z-(SF_M_star_muzzin_hi_z+SF_M_star_upper_err_muzzin_hi_z))*(1+(SF_alpha1_muzzin_hi_z-SF_alpha1_upper_err_muzzin_hi_z)))) * np.exp(-10**(x_muzzin2013_hi_z-(SF_M_star_muzzin_hi_z+SF_M_star_upper_err_muzzin_hi_z)))#
Q_model_muzzin_upper_hi_z = np.log(10) * np.exp(-10**(x_muzzin2013_hi_z-(Q_M_star_muzzin_hi_z+Q_M_star_upper_err_muzzin_hi_z))) * ( ((Q_phi1_muzzin_hi_z+Q_phi1_upper_err_muzzin_hi_z)*(10**(x_muzzin2013_hi_z-(Q_M_star_muzzin_hi_z+Q_M_star_upper_err_muzzin_hi_z)))**(1+(Q_alpha1_muzzin_hi_z-Q_alpha1_upper_err_muzzin_hi_z)))  + ((Q_phi2_muzzin_hi_z+Q_phi2_upper_err_muzzin_hi_z)*(10**(x_muzzin2013_hi_z-(Q_M_star_muzzin_hi_z+Q_M_star_upper_err_muzzin_hi_z)))**(1+(Q_alpha2_muzzin_hi_z-Q_alpha2_upper_err_muzzin_hi_z))) )
T_model_muzzin_upper_hi_z = np.log(10) * (T_phi1_muzzin_hi_z+T_phi1_upper_err_muzzin_hi_z) * (10**((x_muzzin2013_hi_z-(T_M_star_muzzin_hi_z+T_M_star_upper_err_muzzin_hi_z))*(1+(T_alpha1_muzzin_hi_z-T_alpha1_upper_err_muzzin_hi_z)))) * np.exp(-10**(x_muzzin2013_hi_z-(T_M_star_muzzin_hi_z+T_M_star_upper_err_muzzin_hi_z)))#
# lower_error_bound
SF_model_muzzin_lower_hi_z = np.log(10) * (SF_phi1_muzzin_hi_z-SF_phi1_lower_err_muzzin_hi_z) * (10**((x_muzzin2013_hi_z-(SF_M_star_muzzin_hi_z-SF_M_star_lower_err_muzzin_hi_z))*(1+(SF_alpha1_muzzin_hi_z+SF_alpha1_lower_err_muzzin_hi_z)))) * np.exp(-10**(x_muzzin2013_hi_z-(SF_M_star_muzzin_hi_z-SF_M_star_lower_err_muzzin_hi_z)))#
Q_model_muzzin_lower_hi_z = np.log(10) * np.exp(-10**(x_muzzin2013_hi_z-(Q_M_star_muzzin_hi_z-Q_M_star_lower_err_muzzin_hi_z))) * ( ((Q_phi1_muzzin_hi_z-Q_phi1_lower_err_muzzin_hi_z)*(10**(x_muzzin2013_hi_z-(Q_M_star_muzzin_hi_z-Q_M_star_lower_err_muzzin_hi_z)))**(1+(Q_alpha1_muzzin_hi_z+Q_alpha1_lower_err_muzzin_hi_z)))  + ((Q_phi2_muzzin_hi_z-Q_phi2_lower_err_muzzin_hi_z)*(10**(x_muzzin2013_hi_z-(Q_M_star_muzzin_hi_z-Q_M_star_lower_err_muzzin_hi_z)))**(1+(Q_alpha2_muzzin_hi_z+Q_alpha2_lower_err_muzzin_hi_z))) )
T_model_muzzin_lower_hi_z = np.log(10) * (T_phi1_muzzin_hi_z-T_phi1_lower_err_muzzin_hi_z) * (10**((x_muzzin2013_hi_z-(T_M_star_muzzin_hi_z-T_M_star_lower_err_muzzin_hi_z))*(1+(T_alpha1_muzzin_hi_z+T_alpha1_lower_err_muzzin_hi_z)))) * np.exp(-10**(x_muzzin2013_hi_z-(T_M_star_muzzin_hi_z-T_M_star_lower_err_muzzin_hi_z)))#
#
#
# tomczak+2014, lo z
# SMFs
SF_model_tomczak_lo_z = np.log(10) * np.exp(-10**(x_tomczak2014_lo_z-SF_M_star_tomczak_lo_z)) * ( (SF_phi1_tomczak_lo_z*(10**(x_tomczak2014_lo_z-SF_M_star_tomczak_lo_z))**(1+SF_alpha1_tomczak_lo_z))  + (SF_phi2_tomczak_lo_z*(10**(x_tomczak2014_lo_z-SF_M_star_tomczak_lo_z))**(1+SF_alpha2_tomczak_lo_z)) )
Q_model_tomczak_lo_z = np.log(10) * np.exp(-10**(x_tomczak2014_lo_z_Q-Q_M_star_tomczak_lo_z)) * ( (Q_phi1_tomczak_lo_z*(10**(x_tomczak2014_lo_z_Q-Q_M_star_tomczak_lo_z))**(1+Q_alpha1_tomczak_lo_z))  + (Q_phi2_tomczak_lo_z*(10**(x_tomczak2014_lo_z_Q-Q_M_star_tomczak_lo_z))**(1+Q_alpha2_tomczak_lo_z)) )
T_model_tomczak_lo_z = np.log(10) * np.exp(-10**(x_tomczak2014_lo_z-T_M_star_tomczak_lo_z)) * ( (T_phi1_tomczak_lo_z*(10**(x_tomczak2014_lo_z-T_M_star_tomczak_lo_z))**(1+T_alpha1_tomczak_lo_z))  + (T_phi2_tomczak_lo_z*(10**(x_tomczak2014_lo_z-T_M_star_tomczak_lo_z))**(1+T_alpha2_tomczak_lo_z)) )
# upper_error_bound
SF_model_tomczak_upper_lo_z = np.log(10) * np.exp(-10**(x_tomczak2014_lo_z-(SF_M_star_tomczak_lo_z+SF_M_star_upper_err_tomczak_lo_z))) * ( (SF_phi1_upper_err_tomczak_lo_z*(10**(x_tomczak2014_lo_z-(SF_M_star_tomczak_lo_z+SF_M_star_upper_err_tomczak_lo_z)))**(1+(SF_alpha1_tomczak_lo_z-SF_alpha1_upper_err_tomczak_lo_z)))  + (SF_phi2_upper_err_tomczak_lo_z*(10**(x_tomczak2014_lo_z-(SF_M_star_tomczak_lo_z+SF_M_star_upper_err_tomczak_lo_z)))**(1+(SF_alpha2_tomczak_lo_z-SF_alpha1_upper_err_tomczak_lo_z))) )
Q_model_tomczak_upper_lo_z = np.log(10) * np.exp(-10**(x_tomczak2014_lo_z_Q-(Q_M_star_tomczak_lo_z+Q_M_star_upper_err_tomczak_lo_z))) * ( (Q_phi1_upper_err_tomczak_lo_z*(10**(x_tomczak2014_lo_z_Q-(Q_M_star_tomczak_lo_z+Q_M_star_upper_err_tomczak_lo_z)))**(1+(Q_alpha1_tomczak_lo_z-Q_alpha1_upper_err_tomczak_lo_z)))  + (Q_phi2_upper_err_tomczak_lo_z*(10**(x_tomczak2014_lo_z_Q-(Q_M_star_tomczak_lo_z+Q_M_star_upper_err_tomczak_lo_z)))**(1+(Q_alpha2_tomczak_lo_z-Q_alpha1_upper_err_tomczak_lo_z))) )
T_model_tomczak_upper_lo_z = np.log(10) * np.exp(-10**(x_tomczak2014_lo_z-(T_M_star_tomczak_lo_z+T_M_star_upper_err_tomczak_lo_z))) * ( (T_phi1_upper_err_tomczak_lo_z*(10**(x_tomczak2014_lo_z-(T_M_star_tomczak_lo_z+T_M_star_upper_err_tomczak_lo_z)))**(1+(T_alpha1_tomczak_lo_z-T_alpha1_upper_err_tomczak_lo_z)))  + (T_phi2_upper_err_tomczak_lo_z*(10**(x_tomczak2014_lo_z-(T_M_star_tomczak_lo_z+T_M_star_upper_err_tomczak_lo_z)))**(1+(T_alpha2_tomczak_lo_z-T_alpha1_upper_err_tomczak_lo_z))) )
# lower_error_bound
SF_model_tomczak_lower_lo_z = np.log(10) * np.exp(-10**(x_tomczak2014_lo_z-(SF_M_star_tomczak_lo_z-SF_M_star_lower_err_tomczak_lo_z))) * ( (SF_phi1_lower_err_tomczak_lo_z*(10**(x_tomczak2014_lo_z-(SF_M_star_tomczak_lo_z-SF_M_star_lower_err_tomczak_lo_z)))**(1+(SF_alpha1_tomczak_lo_z+SF_alpha1_lower_err_tomczak_lo_z)))  + (SF_phi2_lower_err_tomczak_lo_z*(10**(x_tomczak2014_lo_z-(SF_M_star_tomczak_lo_z-SF_M_star_lower_err_tomczak_lo_z)))**(1+(SF_alpha2_tomczak_lo_z+SF_alpha1_lower_err_tomczak_lo_z))) )
Q_model_tomczak_lower_lo_z = np.log(10) * np.exp(-10**(x_tomczak2014_lo_z_Q-(Q_M_star_tomczak_lo_z-Q_M_star_lower_err_tomczak_lo_z))) * ( (Q_phi1_lower_err_tomczak_lo_z*(10**(x_tomczak2014_lo_z_Q-(Q_M_star_tomczak_lo_z-Q_M_star_lower_err_tomczak_lo_z)))**(1+(Q_alpha1_tomczak_lo_z+Q_alpha1_lower_err_tomczak_lo_z)))  + (Q_phi2_lower_err_tomczak_lo_z*(10**(x_tomczak2014_lo_z_Q-(Q_M_star_tomczak_lo_z-Q_M_star_lower_err_tomczak_lo_z)))**(1+(Q_alpha2_tomczak_lo_z+Q_alpha1_lower_err_tomczak_lo_z))) )
T_model_tomczak_lower_lo_z = np.log(10) * np.exp(-10**(x_tomczak2014_lo_z-(T_M_star_tomczak_lo_z-T_M_star_lower_err_tomczak_lo_z))) * ( (T_phi1_lower_err_tomczak_lo_z*(10**(x_tomczak2014_lo_z-(T_M_star_tomczak_lo_z-T_M_star_lower_err_tomczak_lo_z)))**(1+(T_alpha1_tomczak_lo_z+T_alpha1_lower_err_tomczak_lo_z)))  + (T_phi2_lower_err_tomczak_lo_z*(10**(x_tomczak2014_lo_z-(T_M_star_tomczak_lo_z-T_M_star_lower_err_tomczak_lo_z)))**(1+(T_alpha2_tomczak_lo_z+T_alpha1_lower_err_tomczak_lo_z))) )
# tomczak+2014, hi z
# SMFs
SF_model_tomczak_hi_z = np.log(10) * np.exp(-10**(x_tomczak2014_hi_z-SF_M_star_tomczak_hi_z)) * ( (SF_phi1_tomczak_hi_z*(10**(x_tomczak2014_hi_z-SF_M_star_tomczak_hi_z))**(1+SF_alpha1_tomczak_hi_z))  + (SF_phi2_tomczak_hi_z*(10**(x_tomczak2014_hi_z-SF_M_star_tomczak_hi_z))**(1+SF_alpha2_tomczak_hi_z)) )
Q_model_tomczak_hi_z = np.log(10) * np.exp(-10**(x_tomczak2014_hi_z_Q-Q_M_star_tomczak_hi_z)) * ( (Q_phi1_tomczak_hi_z*(10**(x_tomczak2014_hi_z_Q-Q_M_star_tomczak_hi_z))**(1+Q_alpha1_tomczak_hi_z))  + (Q_phi2_tomczak_hi_z*(10**(x_tomczak2014_hi_z_Q-Q_M_star_tomczak_hi_z))**(1+Q_alpha2_tomczak_hi_z)) )
T_model_tomczak_hi_z = np.log(10) * np.exp(-10**(x_tomczak2014_hi_z-T_M_star_tomczak_hi_z)) * ( (T_phi1_tomczak_hi_z*(10**(x_tomczak2014_hi_z-T_M_star_tomczak_hi_z))**(1+T_alpha1_tomczak_hi_z))  + (T_phi2_tomczak_hi_z*(10**(x_tomczak2014_hi_z-T_M_star_tomczak_hi_z))**(1+T_alpha2_tomczak_hi_z)) )
# upper_error_bound
SF_model_tomczak_upper_hi_z = np.log(10) * np.exp(-10**(x_tomczak2014_hi_z-(SF_M_star_tomczak_hi_z+SF_M_star_upper_err_tomczak_hi_z))) * ( (SF_phi1_upper_err_tomczak_hi_z*(10**(x_tomczak2014_hi_z-(SF_M_star_tomczak_hi_z+SF_M_star_upper_err_tomczak_hi_z)))**(1+(SF_alpha1_tomczak_hi_z-SF_alpha1_upper_err_tomczak_hi_z)))  + (SF_phi2_upper_err_tomczak_hi_z*(10**(x_tomczak2014_hi_z-(SF_M_star_tomczak_hi_z+SF_M_star_upper_err_tomczak_hi_z)))**(1+(SF_alpha2_tomczak_hi_z-SF_alpha1_upper_err_tomczak_hi_z))) )
Q_model_tomczak_upper_hi_z = np.log(10) * np.exp(-10**(x_tomczak2014_hi_z_Q-(Q_M_star_tomczak_hi_z+Q_M_star_upper_err_tomczak_hi_z))) * ( (Q_phi1_upper_err_tomczak_hi_z*(10**(x_tomczak2014_hi_z_Q-(Q_M_star_tomczak_hi_z+Q_M_star_upper_err_tomczak_hi_z)))**(1+(Q_alpha1_tomczak_hi_z-Q_alpha1_upper_err_tomczak_hi_z)))  + (Q_phi2_upper_err_tomczak_hi_z*(10**(x_tomczak2014_hi_z_Q-(Q_M_star_tomczak_hi_z+Q_M_star_upper_err_tomczak_hi_z)))**(1+(Q_alpha2_tomczak_hi_z-Q_alpha1_upper_err_tomczak_hi_z))) )
T_model_tomczak_upper_hi_z = np.log(10) * np.exp(-10**(x_tomczak2014_hi_z-(T_M_star_tomczak_hi_z+T_M_star_upper_err_tomczak_hi_z))) * ( (T_phi1_upper_err_tomczak_hi_z*(10**(x_tomczak2014_hi_z-(T_M_star_tomczak_hi_z+T_M_star_upper_err_tomczak_hi_z)))**(1+(T_alpha1_tomczak_hi_z-T_alpha1_upper_err_tomczak_hi_z)))  + (T_phi2_upper_err_tomczak_hi_z*(10**(x_tomczak2014_hi_z-(T_M_star_tomczak_hi_z+T_M_star_upper_err_tomczak_hi_z)))**(1+(T_alpha2_tomczak_hi_z-T_alpha1_upper_err_tomczak_hi_z))) )
# lower_error_bound
SF_model_tomczak_lower_hi_z = np.log(10) * np.exp(-10**(x_tomczak2014_hi_z-(SF_M_star_tomczak_hi_z-SF_M_star_lower_err_tomczak_hi_z))) * ( (SF_phi1_lower_err_tomczak_hi_z*(10**(x_tomczak2014_hi_z-(SF_M_star_tomczak_hi_z-SF_M_star_lower_err_tomczak_hi_z)))**(1+(SF_alpha1_tomczak_hi_z+SF_alpha1_lower_err_tomczak_hi_z)))  + (SF_phi2_lower_err_tomczak_hi_z*(10**(x_tomczak2014_hi_z-(SF_M_star_tomczak_hi_z-SF_M_star_lower_err_tomczak_hi_z)))**(1+(SF_alpha2_tomczak_hi_z+SF_alpha1_lower_err_tomczak_hi_z))) )
Q_model_tomczak_lower_hi_z = np.log(10) * np.exp(-10**(x_tomczak2014_hi_z_Q-(Q_M_star_tomczak_hi_z-Q_M_star_lower_err_tomczak_hi_z))) * ( (Q_phi1_lower_err_tomczak_hi_z*(10**(x_tomczak2014_hi_z_Q-(Q_M_star_tomczak_hi_z-Q_M_star_lower_err_tomczak_hi_z)))**(1+(Q_alpha1_tomczak_hi_z+Q_alpha1_lower_err_tomczak_hi_z)))  + (Q_phi2_lower_err_tomczak_hi_z*(10**(x_tomczak2014_hi_z_Q-(Q_M_star_tomczak_hi_z-Q_M_star_lower_err_tomczak_hi_z)))**(1+(Q_alpha2_tomczak_hi_z+Q_alpha1_lower_err_tomczak_hi_z))) )
T_model_tomczak_lower_hi_z = np.log(10) * np.exp(-10**(x_tomczak2014_hi_z-(T_M_star_tomczak_hi_z-T_M_star_lower_err_tomczak_hi_z))) * ( (T_phi1_lower_err_tomczak_hi_z*(10**(x_tomczak2014_hi_z-(T_M_star_tomczak_hi_z-T_M_star_lower_err_tomczak_hi_z)))**(1+(T_alpha1_tomczak_hi_z+T_alpha1_lower_err_tomczak_hi_z)))  + (T_phi2_lower_err_tomczak_hi_z*(10**(x_tomczak2014_hi_z-(T_M_star_tomczak_hi_z-T_M_star_lower_err_tomczak_hi_z)))**(1+(T_alpha2_tomczak_hi_z+T_alpha1_lower_err_tomczak_hi_z))) )
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
# baldry+2012
# SMFs
SF_model_baldry2012 = np.log(10) * SF_phi1_baldry2012 * (10**((x_baldry2012-SF_M_star_baldry2012)*(1+SF_alpha1_baldry2012))) * np.exp(-10**(x_baldry2012-SF_M_star_baldry2012))
Q_model_baldry2012 = np.log(10) * np.exp(-10**(x_baldry2012_Q-Q_M_star_baldry2012)) * ( (Q_phi1_baldry2012*(10**(x_baldry2012_Q-Q_M_star_baldry2012))**(1+Q_alpha1_baldry2012))  + (Q_phi2_baldry2012*(10**(x_baldry2012_Q-Q_M_star_baldry2012))**(1+Q_alpha2_baldry2012)) )
T_model_baldry2012 = np.log(10) * np.exp(-10**(x_baldry2012-T_M_star_baldry2012)) * ( (T_phi1_baldry2012*(10**(x_baldry2012-T_M_star_baldry2012))**(1+T_alpha1_baldry2012))  + (T_phi2_baldry2012*(10**(x_baldry2012-T_M_star_baldry2012))**(1+T_alpha2_baldry2012)) )
# upper_error_bound
SF_model_upper_err_baldry2012 = np.log(10) * (SF_phi1_baldry2012+SF_phi1_upper_err_baldry2012) * (10**((x_baldry2012-(SF_M_star_baldry2012+SF_M_star_upper_err_baldry2012))*(1+(SF_alpha1_baldry2012-SF_alpha1_upper_err_baldry2012)))) * np.exp(-10**(x_baldry2012-(SF_M_star_baldry2012+SF_M_star_upper_err_baldry2012)))#
Q_model_upper_err_baldry2012 = np.log(10) * np.exp(-10**(x_baldry2012_Q-(Q_M_star_baldry2012+Q_M_star_upper_err_baldry2012))) * ( ((Q_phi1_baldry2012+Q_phi1_upper_err_baldry2012)*(10**(x_baldry2012_Q-(Q_M_star_baldry2012+Q_M_star_upper_err_baldry2012)))**(1+(Q_alpha1_baldry2012-Q_alpha1_upper_err_baldry2012)))  + ((Q_phi2_baldry2012+Q_phi2_upper_err_baldry2012)*(10**(x_baldry2012_Q-(Q_M_star_baldry2012+Q_M_star_upper_err_baldry2012)))**(1+(Q_alpha2_baldry2012-Q_alpha2_upper_err_baldry2012))) )
T_model_upper_err_baldry2012 = np.log(10) * np.exp(-10**(x_baldry2012-(T_M_star_baldry2012+T_M_star_upper_err_baldry2012))) * ( ((T_phi1_baldry2012+T_phi1_upper_err_baldry2012)*(10**(x_baldry2012-(T_M_star_baldry2012+T_M_star_upper_err_baldry2012)))**(1+(T_alpha1_baldry2012-T_alpha1_upper_err_baldry2012)))  + ((T_phi2_baldry2012+T_phi2_upper_err_baldry2012)*(10**(x_baldry2012-(T_M_star_baldry2012+T_M_star_upper_err_baldry2012)))**(1+(T_alpha2_baldry2012-T_alpha2_upper_err_baldry2012))) )
# lower_error_bound
SF_model_lower_err_baldry2012 = np.log(10) * (SF_phi1_baldry2012-SF_phi1_lower_err_baldry2012) * (10**((x_baldry2012-(SF_M_star_baldry2012-SF_M_star_lower_err_baldry2012))*(1+(SF_alpha1_baldry2012+SF_alpha1_lower_err_baldry2012)))) * np.exp(-10**(x_baldry2012-(SF_M_star_baldry2012-SF_M_star_lower_err_baldry2012)))
Q_model_lower_err_baldry2012 = np.log(10) * np.exp(-10**(x_baldry2012_Q-(Q_M_star_baldry2012-Q_M_star_lower_err_baldry2012))) * ( ((Q_phi1_baldry2012-Q_phi1_lower_err_baldry2012)*(10**(x_baldry2012_Q-(Q_M_star_baldry2012-Q_M_star_lower_err_baldry2012)))**(1+(Q_alpha1_baldry2012+Q_alpha1_lower_err_baldry2012)))  + ((Q_phi2_baldry2012-Q_phi2_lower_err_baldry2012)*(10**(x_baldry2012_Q-(Q_M_star_baldry2012-Q_M_star_lower_err_baldry2012)))**(1+(Q_alpha2_baldry2012+Q_alpha2_lower_err_baldry2012))) )
T_model_lower_err_baldry2012 = np.log(10) * np.exp(-10**(x_baldry2012-(T_M_star_baldry2012-T_M_star_lower_err_baldry2012))) * ( ((T_phi1_baldry2012-T_phi1_lower_err_baldry2012)*(10**(x_baldry2012-(T_M_star_baldry2012-T_M_star_lower_err_baldry2012)))**(1+(T_alpha1_baldry2012+T_alpha1_lower_err_baldry2012)))  + ((T_phi2_baldry2012-T_phi2_lower_err_baldry2012)*(10**(x_baldry2012-(T_M_star_baldry2012-T_M_star_lower_err_baldry2012)))**(1+(T_alpha2_baldry2012+T_alpha2_lower_err_baldry2012))) )
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
    ax0.plot(x_plot,total_model_field_mcmc_plot_double,color='black',linestyle='-',label='this work', linewidth=1.0)
    ax0.plot(x_plot,total_model_field_upper, color='grey', linestyle='--',linewidth=1.0,alpha=0.3)
    ax0.plot(x_plot,total_model_field_lower, color='grey', linestyle='--',linewidth=1.0,alpha=0.3)
    ax0.fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_upper,facecolor='grey',interpolate=True,alpha=0.1)
    ax0.fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_lower,facecolor='grey',interpolate=True,alpha=0.1)
    # muzzin+2013
    ax0.plot(x_muzzin2013_lo_z,T_model_muzzin_lo_z, color='cornflowerblue', linestyle='-',label='Muzzin+13: 0.2 < z < 0.5', linewidth=1.0)
    ax0.plot(x_muzzin2013_lo_z,T_model_muzzin_upper_lo_z, color='cornflowerblue', linestyle='--',linewidth=1.0,alpha=0.3)
    ax0.plot(x_muzzin2013_lo_z,T_model_muzzin_lower_lo_z, color='cornflowerblue', linestyle='--', linewidth=1.0,alpha=0.3)
    ax0.fill_between(x_muzzin2013_lo_z,T_model_muzzin_lo_z,T_model_muzzin_upper_lo_z,facecolor='cornflowerblue',interpolate=True,alpha=0.1)
    ax0.fill_between(x_muzzin2013_lo_z,T_model_muzzin_lo_z,T_model_muzzin_lower_lo_z,facecolor='cornflowerblue',interpolate=True,alpha=0.1)
    ax0.plot(x_muzzin2013_hi_z,T_model_muzzin_hi_z, color='navy', linestyle='-',label='Muzzin+13: 0.5 < z < 1.0', linewidth=1.0)
    ax0.plot(x_muzzin2013_hi_z,T_model_muzzin_upper_hi_z, color='navy', linestyle='--', linewidth=1.0,alpha=0.3)
    ax0.plot(x_muzzin2013_hi_z,T_model_muzzin_lower_hi_z, color='navy', linestyle='--', linewidth=1.0,alpha=0.3)
    ax0.fill_between(x_muzzin2013_hi_z,T_model_muzzin_hi_z,T_model_muzzin_upper_hi_z,facecolor='navy',interpolate=True,alpha=0.1)
    ax0.fill_between(x_muzzin2013_hi_z,T_model_muzzin_hi_z,T_model_muzzin_lower_hi_z,facecolor='navy',interpolate=True,alpha=0.1)
    #tomczak+2014
    ax0.plot(x_tomczak2014_lo_z,T_model_tomczak_lo_z, color='orangered', linestyle='-',label='tomczak+13: 0.2 < z < 0.5', linewidth=1.0)
    ax0.plot(x_tomczak2014_lo_z,T_model_tomczak_upper_lo_z, color='orangered', linestyle='--',linewidth=1.0,alpha=0.3)
    ax0.plot(x_tomczak2014_lo_z,T_model_tomczak_lower_lo_z, color='orangered', linestyle='--', linewidth=1.0,alpha=0.3)
    ax0.fill_between(x_tomczak2014_lo_z,T_model_tomczak_lo_z,T_model_tomczak_upper_lo_z,facecolor='orangered',interpolate=True,alpha=0.1)
    ax0.fill_between(x_tomczak2014_lo_z,T_model_tomczak_lo_z,T_model_tomczak_lower_lo_z,facecolor='orangered',interpolate=True,alpha=0.1)
    ax0.plot(x_tomczak2014_hi_z,T_model_tomczak_hi_z, color='maroon', linestyle='-',label='tomczak+13: 0.5 < z < 1.0', linewidth=1.0)
    ax0.plot(x_tomczak2014_hi_z,T_model_tomczak_upper_hi_z, color='maroon', linestyle='--', linewidth=1.0,alpha=0.3)
    ax0.plot(x_tomczak2014_hi_z,T_model_tomczak_lower_hi_z, color='maroon', linestyle='--', linewidth=1.0,alpha=0.3)
    ax0.fill_between(x_tomczak2014_hi_z,T_model_tomczak_hi_z,T_model_tomczak_upper_hi_z,facecolor='maroon',interpolate=True,alpha=0.1)
    ax0.fill_between(x_tomczak2014_hi_z,T_model_tomczak_hi_z,T_model_tomczak_lower_hi_z,facecolor='maroon',interpolate=True,alpha=0.1)
    # vdb+2018
    ax0.plot(x_vdb2018,T_model_vdb2018, color='magenta', linestyle='-',label='vdB+2018: 0.5 < z < 0.7', linewidth=1.0)
    ax0.plot(x_vdb2018,T_model_upper_err_vdb2018, color='magenta', linestyle='--',linewidth=1.0,alpha=0.3)
    ax0.plot(x_vdb2018,T_model_lower_err_vdb2018, color='magenta', linestyle='--', linewidth=1.0,alpha=0.3)
    ax0.fill_between(x_vdb2018,T_model_vdb2018,T_model_upper_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
    ax0.fill_between(x_vdb2018,T_model_vdb2018,T_model_lower_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
    # baldry+2012
    ax0.plot(x_baldry2012,T_model_baldry2012, color='cyan', linestyle='-',label='baldry+2012: z < 0.06', linewidth=1.0)
    ax0.plot(x_baldry2012,T_model_upper_err_baldry2012, color='cyan', linestyle='--',linewidth=1.0,alpha=0.3)
    ax0.plot(x_baldry2012,T_model_lower_err_baldry2012, color='cyan', linestyle='--', linewidth=1.0,alpha=0.3)
    ax0.fill_between(x_baldry2012,T_model_baldry2012,T_model_upper_err_baldry2012,facecolor='cyan',interpolate=True,alpha=0.1)
    ax0.fill_between(x_baldry2012,T_model_baldry2012,T_model_lower_err_baldry2012,facecolor='cyan',interpolate=True,alpha=0.1)
    #
    ax0.set_xscale('linear')
    ax0.set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    ax0.set_xlim(7.8,12.5)
    ax0.set_yscale('log')
    ax0.set_ylim(1e-6,0.2)
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
    ax1.plot(x_plot,SF_model_field_mcmc_plot,color='black',linestyle='-',linewidth=1.0)
    ax1.plot(x_plot,SF_model_field_upper,color='grey',linestyle='--',linewidth=1.0)
    ax1.plot(x_plot,SF_model_field_lower,color='grey',linestyle='--',linewidth=1.0)
    ax1.fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_upper,facecolor='grey',interpolate=True,alpha=0.1)
    ax1.fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_lower,facecolor='grey',interpolate=True,alpha=0.1)
    # muzzin+2013
    ax1.plot(x_muzzin2013_lo_z,SF_model_muzzin_lo_z, color='cornflowerblue', linestyle='-',label='Muzzin+13: 0.2 < z < 0.5', linewidth=1.0)
    ax1.plot(x_muzzin2013_lo_z,SF_model_muzzin_upper_lo_z, color='cornflowerblue', linestyle='--',linewidth=1.0,alpha=0.3)
    ax1.plot(x_muzzin2013_lo_z,SF_model_muzzin_lower_lo_z, color='cornflowerblue', linestyle='--', linewidth=1.0,alpha=0.3)
    ax1.fill_between(x_muzzin2013_lo_z,SF_model_muzzin_lo_z,SF_model_muzzin_upper_lo_z,facecolor='cornflowerblue',interpolate=True,alpha=0.1)
    ax1.fill_between(x_muzzin2013_lo_z,SF_model_muzzin_lo_z,SF_model_muzzin_lower_lo_z,facecolor='cornflowerblue',interpolate=True,alpha=0.1)
    ax1.plot(x_muzzin2013_hi_z,SF_model_muzzin_hi_z, color='navy', linestyle='-',label='Muzzin+13: 0.5 < z < 1.0', linewidth=1.0)
    ax1.plot(x_muzzin2013_hi_z,SF_model_muzzin_upper_hi_z, color='navy', linestyle='--', linewidth=1.0,alpha=0.3)
    ax1.plot(x_muzzin2013_hi_z,SF_model_muzzin_lower_hi_z, color='navy', linestyle='--', linewidth=1.0,alpha=0.3)
    ax1.fill_between(x_muzzin2013_hi_z,SF_model_muzzin_hi_z,SF_model_muzzin_upper_hi_z,facecolor='navy',interpolate=True,alpha=0.1)
    ax1.fill_between(x_muzzin2013_hi_z,SF_model_muzzin_hi_z,SF_model_muzzin_lower_hi_z,facecolor='navy',interpolate=True,alpha=0.1)
    # tomczak+2014
    ax1.plot(x_tomczak2014_lo_z,SF_model_tomczak_lo_z, color='orangered', linestyle='-',label='tomczak+13: 0.2 < z < 0.5', linewidth=1.0)
    ax1.plot(x_tomczak2014_lo_z,SF_model_tomczak_upper_lo_z, color='orangered', linestyle='--',linewidth=1.0,alpha=0.3)
    ax1.plot(x_tomczak2014_lo_z,SF_model_tomczak_lower_lo_z, color='orangered', linestyle='--', linewidth=1.0,alpha=0.3)
    ax1.fill_between(x_tomczak2014_lo_z,SF_model_tomczak_lo_z,SF_model_tomczak_upper_lo_z,facecolor='orangered',interpolate=True,alpha=0.1)
    ax1.fill_between(x_tomczak2014_lo_z,SF_model_tomczak_lo_z,SF_model_tomczak_lower_lo_z,facecolor='orangered',interpolate=True,alpha=0.1)
    ax1.plot(x_tomczak2014_hi_z,SF_model_tomczak_hi_z, color='maroon', linestyle='-',label='tomczak+13: 0.5 < z < 1.0', linewidth=1.0)
    ax1.plot(x_tomczak2014_hi_z,SF_model_tomczak_upper_hi_z, color='maroon', linestyle='--', linewidth=1.0,alpha=0.3)
    ax1.plot(x_tomczak2014_hi_z,SF_model_tomczak_lower_hi_z, color='maroon', linestyle='--', linewidth=1.0,alpha=0.3)
    ax1.fill_between(x_tomczak2014_hi_z,SF_model_tomczak_hi_z,SF_model_tomczak_upper_hi_z,facecolor='maroon',interpolate=True,alpha=0.1)
    ax1.fill_between(x_tomczak2014_hi_z,SF_model_tomczak_hi_z,SF_model_tomczak_lower_hi_z,facecolor='maroon',interpolate=True,alpha=0.1)
    # vdb+2018
    ax1.plot(x_vdb2018,SF_model_vdb2018, color='magenta', linestyle='-',label='vdB+2018: 0.5 < z < 0.7', linewidth=1.0)
    ax1.plot(x_vdb2018,SF_model_upper_err_vdb2018, color='magenta', linestyle='--',linewidth=1.0,alpha=0.3)
    ax1.plot(x_vdb2018,SF_model_lower_err_vdb2018, color='magenta', linestyle='--', linewidth=1.0,alpha=0.3)
    ax1.fill_between(x_vdb2018,SF_model_vdb2018,SF_model_upper_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
    ax1.fill_between(x_vdb2018,SF_model_vdb2018,SF_model_lower_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
    # baldry+2012
    ax1.plot(x_baldry2012,SF_model_baldry2012, color='cyan', linestyle='-',label='baldry+2012: z < 0.06', linewidth=1.0)
    ax1.plot(x_baldry2012,SF_model_upper_err_baldry2012, color='cyan', linestyle='--',linewidth=1.0,alpha=0.3)
    ax1.plot(x_baldry2012,SF_model_lower_err_baldry2012, color='cyan', linestyle='--', linewidth=1.0,alpha=0.3)
    ax1.fill_between(x_baldry2012,SF_model_baldry2012,SF_model_upper_err_baldry2012,facecolor='cyan',interpolate=True,alpha=0.1)
    ax1.fill_between(x_baldry2012,SF_model_baldry2012,SF_model_lower_err_baldry2012,facecolor='cyan',interpolate=True,alpha=0.1)
    #
    ax1.set_xscale('linear')
    ax1.set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    ax1.set_xlim(7.8,12.5)
    ax1.set_yscale('log')
    ax1.set_ylim(1e-6,0.2)
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
    ax2.plot(x_plot,Q_model_field_mcmc_plot_double,color='black',linestyle='-',label='this work', linewidth=1.0)
    ax2.plot(x_plot,Q_model_field_upper,color='grey',linestyle='--', linewidth=1.0)
    ax2.plot(x_plot,Q_model_field_lower,color='grey',linestyle='--',linewidth=1.0)
    ax2.fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_upper,facecolor='grey',interpolate=True,alpha=0.1)
    ax2.fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_lower,facecolor='grey',interpolate=True,alpha=0.1)
    # muzzin+2013
    ax2.plot(x_muzzin2013_lo_z,Q_model_muzzin_lo_z, color='cornflowerblue', linestyle='-',label='Muzzin+2013: 0.2 < z < 0.5', linewidth=1.0)
    ax2.plot(x_muzzin2013_lo_z,Q_model_muzzin_upper_lo_z, color='cornflowerblue', linestyle='--',linewidth=1.0,alpha=0.3)
    ax2.plot(x_muzzin2013_lo_z,Q_model_muzzin_lower_lo_z, color='cornflowerblue', linestyle='--', linewidth=1.0,alpha=0.3)
    ax2.fill_between(x_muzzin2013_lo_z,Q_model_muzzin_lo_z,Q_model_muzzin_upper_lo_z,facecolor='cornflowerblue',interpolate=True,alpha=0.1)
    ax2.fill_between(x_muzzin2013_lo_z,Q_model_muzzin_lo_z,Q_model_muzzin_lower_lo_z,facecolor='cornflowerblue',interpolate=True,alpha=0.1)
    ax2.plot(x_muzzin2013_hi_z,Q_model_muzzin_hi_z, color='navy', linestyle='-',label='Muzzin+2013: 0.5 < z < 1.0', linewidth=1.0)
    ax2.plot(x_muzzin2013_hi_z,Q_model_muzzin_upper_hi_z, color='navy', linestyle='--', linewidth=1.0,alpha=0.3)
    ax2.plot(x_muzzin2013_hi_z,Q_model_muzzin_lower_hi_z, color='navy', linestyle='--', linewidth=1.0,alpha=0.3)
    ax2.fill_between(x_muzzin2013_hi_z,Q_model_muzzin_hi_z,Q_model_muzzin_upper_hi_z,facecolor='navy',interpolate=True,alpha=0.1)
    ax2.fill_between(x_muzzin2013_hi_z,Q_model_muzzin_hi_z,Q_model_muzzin_lower_hi_z,facecolor='navy',interpolate=True,alpha=0.1)
    # tomczak+2014
    ax2.plot(x_tomczak2014_lo_z_Q,Q_model_tomczak_lo_z, color='orangered', linestyle='-',label='tomczak+2014: 0.2 < z < 0.5', linewidth=1.0)
    ax2.plot(x_tomczak2014_lo_z_Q,Q_model_tomczak_upper_lo_z, color='orangered', linestyle='--',linewidth=1.0,alpha=0.3)
    ax2.plot(x_tomczak2014_lo_z_Q,Q_model_tomczak_lower_lo_z, color='orangered', linestyle='--', linewidth=1.0,alpha=0.3)
    ax2.fill_between(x_tomczak2014_lo_z_Q,Q_model_tomczak_lo_z,Q_model_tomczak_upper_lo_z,facecolor='orangered',interpolate=True,alpha=0.1)
    ax2.fill_between(x_tomczak2014_lo_z_Q,Q_model_tomczak_lo_z,Q_model_tomczak_lower_lo_z,facecolor='orangered',interpolate=True,alpha=0.1)
    ax2.plot(x_tomczak2014_hi_z_Q,Q_model_tomczak_hi_z, color='maroon', linestyle='-',label='tomczak+2014: 0.5 < z < 0.75', linewidth=1.0)
    ax2.plot(x_tomczak2014_hi_z_Q,Q_model_tomczak_upper_hi_z, color='maroon', linestyle='--', linewidth=1.0,alpha=0.3)
    ax2.plot(x_tomczak2014_hi_z_Q,Q_model_tomczak_lower_hi_z, color='maroon', linestyle='--', linewidth=1.0,alpha=0.3)
    ax2.fill_between(x_tomczak2014_hi_z_Q,Q_model_tomczak_hi_z,Q_model_tomczak_upper_hi_z,facecolor='maroon',interpolate=True,alpha=0.1)
    ax2.fill_between(x_tomczak2014_hi_z_Q,Q_model_tomczak_hi_z,Q_model_tomczak_lower_hi_z,facecolor='maroon',interpolate=True,alpha=0.1)
    # vdb+2018
    ax2.plot(x_vdb2018,Q_model_vdb2018, color='magenta', linestyle='-',label='vdB+2018: 0.5 < z < 0.7', linewidth=1.0)
    ax2.plot(x_vdb2018,Q_model_upper_err_vdb2018, color='magenta', linestyle='--',linewidth=1.0,alpha=0.3)
    ax2.plot(x_vdb2018,Q_model_lower_err_vdb2018, color='magenta', linestyle='--', linewidth=1.0,alpha=0.3)
    ax2.fill_between(x_vdb2018,Q_model_vdb2018,Q_model_upper_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
    ax2.fill_between(x_vdb2018,Q_model_vdb2018,Q_model_lower_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
    # baldry+2012
    ax2.plot(x_baldry2012_Q,Q_model_baldry2012, color='cyan', linestyle='-',label='baldry+2012: z < 0.06', linewidth=1.0)
    ax2.plot(x_baldry2012_Q,Q_model_upper_err_baldry2012, color='cyan', linestyle='--',linewidth=1.0,alpha=0.3)
    ax2.plot(x_baldry2012_Q,Q_model_lower_err_baldry2012, color='cyan', linestyle='--', linewidth=1.0,alpha=0.3)
    ax2.fill_between(x_baldry2012_Q,Q_model_baldry2012,Q_model_upper_err_baldry2012,facecolor='cyan',interpolate=True,alpha=0.1)
    ax2.fill_between(x_baldry2012_Q,Q_model_baldry2012,Q_model_lower_err_baldry2012,facecolor='cyan',interpolate=True,alpha=0.1)
    #
    ax2.set_xscale('linear')
    ax2.set_xlabel('log($M_{*}$/$M_{\odot}$)',fontsize=20)
    ax2.set_xlim(7.8,12.5)
    ax2.set_yscale('log')
    ax2.set_ylim(1e-6,0.2)
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
if plot_flag_3 == 1:
    #
    ## TOTAL population
    fig, axs = plt.subplots(2,2, sharex=True,sharey=True,gridspec_kw={'hspace': 0,'wspace': 0})
    #
    fig.suptitle('Field: Total',fontsize=30)
    # Panel 1: Baldry+2012
    # my fits
    axs[0,0].plot(x_plot,total_model_field_mcmc_plot_double,color='black',linestyle='-',label='this work', linewidth=1.0)
    axs[0,0].plot(x_plot,total_model_field_upper, color='grey', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].plot(x_plot,total_model_field_lower, color='grey', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_upper,facecolor='grey',interpolate=True,alpha=0.1)
    axs[0,0].fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_lower,facecolor='grey',interpolate=True,alpha=0.1)
    # baldry+2012 data pts + schechter fits
    axs[0,0].errorbar(x_baldry2012_data,total_smf_baldry2012,yerr=total_err_baldry2012,marker='s',color='cyan',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[0,0].plot(x_baldry2012,T_model_baldry2012, color='cyan', linestyle='-',label='baldry+2012: z < 0.06\n(GAMA)', linewidth=1.0)
    axs[0,0].plot(x_baldry2012,T_model_upper_err_baldry2012, color='cyan', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].plot(x_baldry2012,T_model_lower_err_baldry2012, color='cyan', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,0].fill_between(x_baldry2012,T_model_baldry2012,T_model_upper_err_baldry2012,facecolor='cyan',interpolate=True,alpha=0.1)
    axs[0,0].fill_between(x_baldry2012,T_model_baldry2012,T_model_lower_err_baldry2012,facecolor='cyan',interpolate=True,alpha=0.1)
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
    axs[0,1].plot(x_plot,total_model_field_mcmc_plot_double,color='black',linestyle='-',label='this work', linewidth=1.0)
    axs[0,1].plot(x_plot,total_model_field_upper, color='grey', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_plot,total_model_field_lower, color='grey', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_upper,facecolor='grey',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_lower,facecolor='grey',interpolate=True,alpha=0.1)
    # muzzin+2013 data + fits
    axs[0,1].plot(x_muzzin2013_lo_z,T_model_muzzin_lo_z, color='cornflowerblue', linestyle='-',label='muzzin+2013: 0.2 < z < 0.5', linewidth=1.0)
    axs[0,1].plot(x_muzzin2013_lo_z,T_model_muzzin_upper_lo_z, color='cornflowerblue', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_muzzin2013_lo_z,T_model_muzzin_lower_lo_z, color='cornflowerblue', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_muzzin2013_lo_z,T_model_muzzin_lo_z,T_model_muzzin_upper_lo_z,facecolor='cornflowerblue',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_muzzin2013_lo_z,T_model_muzzin_lo_z,T_model_muzzin_lower_lo_z,facecolor='cornflowerblue',interpolate=True,alpha=0.1)
    axs[0,1].plot(x_muzzin2013_hi_z,T_model_muzzin_hi_z, color='navy', linestyle='-',label='muzzin+2013: 0.5 < z < 1.0\n(COSMOS/UltraVISTA)', linewidth=1.0)
    axs[0,1].plot(x_muzzin2013_hi_z,T_model_muzzin_upper_hi_z, color='navy', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_muzzin2013_hi_z,T_model_muzzin_lower_hi_z, color='navy', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_muzzin2013_hi_z,T_model_muzzin_hi_z,T_model_muzzin_upper_hi_z,facecolor='navy',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_muzzin2013_hi_z,T_model_muzzin_hi_z,T_model_muzzin_lower_hi_z,facecolor='navy',interpolate=True,alpha=0.1)
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
    axs[1,0].plot(x_plot,total_model_field_mcmc_plot_double,color='black',linestyle='-',label='this work', linewidth=1.0)
    axs[1,0].plot(x_plot,total_model_field_upper, color='grey', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_plot,total_model_field_lower, color='grey', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_upper,facecolor='grey',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_lower,facecolor='grey',interpolate=True,alpha=0.1)
    # tomczak+2014 data + fits
    axs[1,0].errorbar(x_tom2014_data,total_smf_tom2014_lo_z,yerr=total_err_tom2014_lo_z,marker='<',color='orangered',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[1,0].errorbar(x_tom2014_data,total_smf_tom2014_hi_z,yerr=total_err_tom2014_hi_z,marker='>',color='maroon',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[1,0].plot(x_tomczak2014_lo_z,T_model_tomczak_lo_z, color='orangered', linestyle='-',label='tomczak+2014: 0.2 < z < 0.5', linewidth=1.0)
    axs[1,0].plot(x_tomczak2014_lo_z,T_model_tomczak_upper_lo_z, color='orangered', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_tomczak2014_lo_z,T_model_tomczak_lower_lo_z, color='orangered', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_tomczak2014_lo_z,T_model_tomczak_lo_z,T_model_tomczak_upper_lo_z,facecolor='orangered',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_tomczak2014_lo_z,T_model_tomczak_lo_z,T_model_tomczak_lower_lo_z,facecolor='orangered',interpolate=True,alpha=0.1)
    axs[1,0].plot(x_tomczak2014_hi_z,T_model_tomczak_hi_z, color='maroon', linestyle='-',label='tomczak+2014: 0.5 < z < 0.75\n(ZFOURGE/CANDELS)', linewidth=1.0)
    axs[1,0].plot(x_tomczak2014_hi_z,T_model_tomczak_upper_hi_z, color='maroon', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_tomczak2014_hi_z,T_model_tomczak_lower_hi_z, color='maroon', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_tomczak2014_hi_z,T_model_tomczak_hi_z,T_model_tomczak_upper_hi_z,facecolor='maroon',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_tomczak2014_hi_z,T_model_tomczak_hi_z,T_model_tomczak_lower_hi_z,facecolor='maroon',interpolate=True,alpha=0.1)
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
    # Panel 4: vanderBurg+2018
    # my fits
    axs[1,1].plot(x_plot,total_model_field_mcmc_plot_double,color='black',linestyle='-',label='this work', linewidth=1.0)
    axs[1,1].plot(x_plot,total_model_field_upper, color='grey', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,1].plot(x_plot,total_model_field_lower, color='grey', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,1].fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_upper,facecolor='grey',interpolate=True,alpha=0.1)
    axs[1,1].fill_between(x_plot,total_model_field_mcmc_plot_double,total_model_field_lower,facecolor='grey',interpolate=True,alpha=0.1)
    # baldry+2012 data pts + schechter fits
    axs[1,1].errorbar(x_vdb2018_data,total_smf_vdb2018,yerr=total_err_vdb2018,marker='s',color='magenta',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[1,1].plot(x_vdb2018,T_model_vdb2018, color='magenta', linestyle='-',label='vdB+2018: 0.5 < z < 0.7\n(COSMOS/UltraVISTA)', linewidth=1.0)
    axs[1,1].plot(x_vdb2018,T_model_upper_err_vdb2018, color='magenta', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,1].plot(x_vdb2018,T_model_lower_err_vdb2018, color='magenta', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,1].fill_between(x_vdb2018,T_model_vdb2018,T_model_upper_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
    axs[1,1].fill_between(x_vdb2018,T_model_vdb2018,T_model_lower_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
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
    axs[0,0].plot(x_plot,SF_model_field_mcmc_plot,color='b',linestyle='-',label='this work', linewidth=1.0)
    axs[0,0].plot(x_plot,SF_model_field_upper, color='b', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].plot(x_plot,SF_model_field_lower, color='b', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_upper,facecolor='grey',interpolate=True,alpha=0.1)
    axs[0,0].fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_lower,facecolor='grey',interpolate=True,alpha=0.1)
    # baldry+2012 data pts + schechter fits
    # axs[0,0].errorbar(x_baldry2012_data,total_smf_baldry2012,yerr=total_err_baldry2012,marker='s',color='cyan',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[0,0].plot(x_baldry2012,SF_model_baldry2012, color='cyan', linestyle='-',label='baldry+2012: z < 0.06\n(GAMA)', linewidth=1.0)
    axs[0,0].plot(x_baldry2012,SF_model_upper_err_baldry2012, color='cyan', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].plot(x_baldry2012,SF_model_lower_err_baldry2012, color='cyan', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,0].fill_between(x_baldry2012,SF_model_baldry2012,SF_model_upper_err_baldry2012,facecolor='cyan',interpolate=True,alpha=0.1)
    axs[0,0].fill_between(x_baldry2012,SF_model_baldry2012,SF_model_lower_err_baldry2012,facecolor='cyan',interpolate=True,alpha=0.1)
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
    axs[0,1].plot(x_plot,SF_model_field_mcmc_plot,color='b',linestyle='-',label='this work', linewidth=1.0)
    axs[0,1].plot(x_plot,SF_model_field_upper, color='b', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_plot,SF_model_field_lower, color='b', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_upper,facecolor='grey',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_lower,facecolor='grey',interpolate=True,alpha=0.1)
    # muzzin+2013 data + fits
    axs[0,1].plot(x_muzzin2013_lo_z,SF_model_muzzin_lo_z, color='cornflowerblue', linestyle='-',label='muzzin+2013: 0.2 < z < 0.5', linewidth=1.0)
    axs[0,1].plot(x_muzzin2013_lo_z,SF_model_muzzin_upper_lo_z, color='cornflowerblue', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_muzzin2013_lo_z,SF_model_muzzin_lower_lo_z, color='cornflowerblue', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_muzzin2013_lo_z,SF_model_muzzin_lo_z,SF_model_muzzin_upper_lo_z,facecolor='cornflowerblue',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_muzzin2013_lo_z,SF_model_muzzin_lo_z,SF_model_muzzin_lower_lo_z,facecolor='cornflowerblue',interpolate=True,alpha=0.1)
    axs[0,1].plot(x_muzzin2013_hi_z,SF_model_muzzin_hi_z, color='navy', linestyle='-',label='muzzin+2013: 0.5 < z < 1.0\n(COSMOS/UltraVISTA)', linewidth=1.0)
    axs[0,1].plot(x_muzzin2013_hi_z,SF_model_muzzin_upper_hi_z, color='navy', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_muzzin2013_hi_z,SF_model_muzzin_lower_hi_z, color='navy', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_muzzin2013_hi_z,SF_model_muzzin_hi_z,SF_model_muzzin_upper_hi_z,facecolor='navy',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_muzzin2013_hi_z,SF_model_muzzin_hi_z,SF_model_muzzin_lower_hi_z,facecolor='navy',interpolate=True,alpha=0.1)
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
    axs[1,0].plot(x_plot,SF_model_field_mcmc_plot,color='b',linestyle='-',label='this work', linewidth=1.0)
    axs[1,0].plot(x_plot,SF_model_field_upper, color='b', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_plot,SF_model_field_lower, color='b', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_upper,facecolor='grey',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_lower,facecolor='grey',interpolate=True,alpha=0.1)
    # tomczak+2014 data + fits
    axs[1,0].errorbar(x_tom2014_data,SF_smf_tom2014_lo_z,yerr=SF_err_tom2014_lo_z,marker='<',color='orangered',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[1,0].errorbar(x_tom2014_data,SF_smf_tom2014_hi_z,yerr=SF_err_tom2014_hi_z,marker='>',color='maroon',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[1,0].plot(x_tomczak2014_lo_z,SF_model_tomczak_lo_z, color='orangered', linestyle='-',label='tomczak+2014: 0.2 < z < 0.5', linewidth=1.0)
    axs[1,0].plot(x_tomczak2014_lo_z,SF_model_tomczak_upper_lo_z, color='orangered', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_tomczak2014_lo_z,SF_model_tomczak_lower_lo_z, color='orangered', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_tomczak2014_lo_z,SF_model_tomczak_lo_z,SF_model_tomczak_upper_lo_z,facecolor='orangered',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_tomczak2014_lo_z,SF_model_tomczak_lo_z,SF_model_tomczak_lower_lo_z,facecolor='orangered',interpolate=True,alpha=0.1)
    axs[1,0].plot(x_tomczak2014_hi_z,SF_model_tomczak_hi_z, color='maroon', linestyle='-',label='tomczak+2014: 0.5 < z < 0.75\n(ZFOURGE/CANDELS)', linewidth=1.0)
    axs[1,0].plot(x_tomczak2014_hi_z,SF_model_tomczak_upper_hi_z, color='maroon', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_tomczak2014_hi_z,SF_model_tomczak_lower_hi_z, color='maroon', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_tomczak2014_hi_z,SF_model_tomczak_hi_z,SF_model_tomczak_upper_hi_z,facecolor='maroon',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_tomczak2014_hi_z,SF_model_tomczak_hi_z,SF_model_tomczak_lower_hi_z,facecolor='maroon',interpolate=True,alpha=0.1)
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
    # Panel 4: vanderBurg+2018
    # my fits
    axs[1,1].plot(x_plot,SF_model_field_mcmc_plot,color='b',linestyle='-',label='this work', linewidth=1.0)
    axs[1,1].plot(x_plot,SF_model_field_upper, color='b', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,1].plot(x_plot,SF_model_field_lower, color='b', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,1].fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_upper,facecolor='grey',interpolate=True,alpha=0.1)
    axs[1,1].fill_between(x_plot,SF_model_field_mcmc_plot,SF_model_field_lower,facecolor='grey',interpolate=True,alpha=0.1)
    # vdb+2018 data pts + schechter fits
    axs[1,1].errorbar(x_vdb2018_data,SF_smf_vdb2018,yerr=total_err_vdb2018,marker='s',color='magenta',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[1,1].plot(x_vdb2018,SF_model_vdb2018, color='magenta', linestyle='-',label='vdB+2018: 0.5 < z < 0.7\n(COSMOS/UltraVISTA)', linewidth=1.0)
    axs[1,1].plot(x_vdb2018,SF_model_upper_err_vdb2018, color='magenta', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,1].plot(x_vdb2018,SF_model_lower_err_vdb2018, color='magenta', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,1].fill_between(x_vdb2018,SF_model_vdb2018,SF_model_upper_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
    axs[1,1].fill_between(x_vdb2018,SF_model_vdb2018,SF_model_lower_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
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
    axs[0,0].plot(x_plot,Q_model_field_mcmc_plot_double,color='r',linestyle='-',label='this work', linewidth=1.0)
    axs[0,0].plot(x_plot,Q_model_field_upper, color='r', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].plot(x_plot,Q_model_field_lower, color='r', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_upper,facecolor='grey',interpolate=True,alpha=0.1)
    axs[0,0].fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_lower,facecolor='grey',interpolate=True,alpha=0.1)
    # baldry+2012 data pts + schechter fits
    # axs[0,0].errorbar(x_baldry2012_data,total_smf_baldry2012,yerr=total_err_baldry2012,marker='s',color='cyan',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[0,0].plot(x_baldry2012_Q,Q_model_baldry2012, color='cyan', linestyle='-',label='baldry+2012: z < 0.06\n(GAMA)', linewidth=1.0)
    axs[0,0].plot(x_baldry2012_Q,Q_model_upper_err_baldry2012, color='cyan', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,0].plot(x_baldry2012_Q,Q_model_lower_err_baldry2012, color='cyan', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,0].fill_between(x_baldry2012_Q,Q_model_baldry2012,Q_model_upper_err_baldry2012,facecolor='cyan',interpolate=True,alpha=0.1)
    axs[0,0].fill_between(x_baldry2012_Q,Q_model_baldry2012,Q_model_lower_err_baldry2012,facecolor='cyan',interpolate=True,alpha=0.1)
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
    axs[0,1].plot(x_plot,Q_model_field_mcmc_plot_double,color='r',linestyle='-',label='this work', linewidth=1.0)
    axs[0,1].plot(x_plot,Q_model_field_upper, color='r', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_plot,Q_model_field_lower, color='r', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_upper,facecolor='grey',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_lower,facecolor='grey',interpolate=True,alpha=0.1)
    # muzzin+2013 data + fits
    axs[0,1].plot(x_muzzin2013_lo_z,Q_model_muzzin_lo_z, color='cornflowerblue', linestyle='-',label='muzzin+2013: 0.2 < z < 0.5', linewidth=1.0)
    axs[0,1].plot(x_muzzin2013_lo_z,Q_model_muzzin_upper_lo_z, color='cornflowerblue', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_muzzin2013_lo_z,Q_model_muzzin_lower_lo_z, color='cornflowerblue', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_muzzin2013_lo_z,Q_model_muzzin_lo_z,Q_model_muzzin_upper_lo_z,facecolor='cornflowerblue',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_muzzin2013_lo_z,Q_model_muzzin_lo_z,Q_model_muzzin_lower_lo_z,facecolor='cornflowerblue',interpolate=True,alpha=0.1)
    axs[0,1].plot(x_muzzin2013_hi_z,Q_model_muzzin_hi_z, color='navy', linestyle='-',label='muzzin+2013: 0.5 < z < 1.0\n(COSMOS/UltraVISTA)', linewidth=1.0)
    axs[0,1].plot(x_muzzin2013_hi_z,Q_model_muzzin_upper_hi_z, color='navy', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].plot(x_muzzin2013_hi_z,Q_model_muzzin_lower_hi_z, color='navy', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[0,1].fill_between(x_muzzin2013_hi_z,Q_model_muzzin_hi_z,Q_model_muzzin_upper_hi_z,facecolor='navy',interpolate=True,alpha=0.1)
    axs[0,1].fill_between(x_muzzin2013_hi_z,Q_model_muzzin_hi_z,Q_model_muzzin_lower_hi_z,facecolor='navy',interpolate=True,alpha=0.1)
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
    axs[1,0].plot(x_plot,Q_model_field_mcmc_plot_double,color='r',linestyle='-',label='this work', linewidth=1.0)
    axs[1,0].plot(x_plot,Q_model_field_upper, color='r', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_plot,Q_model_field_lower, color='r', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_upper,facecolor='grey',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_lower,facecolor='grey',interpolate=True,alpha=0.1)
    # tomczak+2014 data + fits
    axs[1,0].errorbar(x_tom2014_data,Q_smf_tom2014_lo_z,yerr=Q_err_tom2014_lo_z,marker='<',color='orangered',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[1,0].errorbar(x_tom2014_data,Q_smf_tom2014_hi_z,yerr=Q_err_tom2014_hi_z,marker='>',color='maroon',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[1,0].plot(x_tomczak2014_lo_z_Q,Q_model_tomczak_lo_z, color='orangered', linestyle='-',label='tomczak+2014: 0.2 < z < 0.5', linewidth=1.0)
    axs[1,0].plot(x_tomczak2014_lo_z_Q,Q_model_tomczak_upper_lo_z, color='orangered', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_tomczak2014_lo_z_Q,Q_model_tomczak_lower_lo_z, color='orangered', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_tomczak2014_lo_z_Q,Q_model_tomczak_lo_z,Q_model_tomczak_upper_lo_z,facecolor='orangered',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_tomczak2014_lo_z_Q,Q_model_tomczak_lo_z,Q_model_tomczak_lower_lo_z,facecolor='orangered',interpolate=True,alpha=0.1)
    axs[1,0].plot(x_tomczak2014_hi_z_Q,Q_model_tomczak_hi_z, color='maroon', linestyle='-',label='tomczak+2014: 0.5 < z < 0.75\n(ZFOURGE/CANDELS)', linewidth=1.0)
    axs[1,0].plot(x_tomczak2014_hi_z_Q,Q_model_tomczak_upper_hi_z, color='maroon', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].plot(x_tomczak2014_hi_z_Q,Q_model_tomczak_lower_hi_z, color='maroon', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,0].fill_between(x_tomczak2014_hi_z_Q,Q_model_tomczak_hi_z,Q_model_tomczak_upper_hi_z,facecolor='maroon',interpolate=True,alpha=0.1)
    axs[1,0].fill_between(x_tomczak2014_hi_z_Q,Q_model_tomczak_hi_z,Q_model_tomczak_lower_hi_z,facecolor='maroon',interpolate=True,alpha=0.1)
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
    # Panel 4: vanderBurg+2018
    # my fits
    axs[1,1].plot(x_plot,Q_model_field_mcmc_plot_double,color='r',linestyle='-',label='this work', linewidth=1.0)
    axs[1,1].plot(x_plot,Q_model_field_upper, color='r', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,1].plot(x_plot,Q_model_field_lower, color='r', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,1].fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_upper,facecolor='grey',interpolate=True,alpha=0.1)
    axs[1,1].fill_between(x_plot,Q_model_field_mcmc_plot_double,Q_model_field_lower,facecolor='grey',interpolate=True,alpha=0.1)
    # vdb+2018 data pts + schechter fits
    axs[1,1].errorbar(x_vdb2018_data,Q_smf_vdb2018,yerr=total_err_vdb2018,marker='s',color='magenta',fillstyle='none',lolims=False, uplims=False, linewidth=0.0, elinewidth=0.5)
    axs[1,1].plot(x_vdb2018,Q_model_vdb2018, color='magenta', linestyle='-',label='vdB+2018: 0.5 < z < 0.7\n(COSMOS/UltraVISTA)', linewidth=1.0)
    axs[1,1].plot(x_vdb2018,Q_model_upper_err_vdb2018, color='magenta', linestyle='--',linewidth=1.0,alpha=0.3)
    axs[1,1].plot(x_vdb2018,Q_model_lower_err_vdb2018, color='magenta', linestyle='--', linewidth=1.0,alpha=0.3)
    axs[1,1].fill_between(x_vdb2018,Q_model_vdb2018,Q_model_upper_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
    axs[1,1].fill_between(x_vdb2018,Q_model_vdb2018,Q_model_lower_err_vdb2018,facecolor='magenta',interpolate=True,alpha=0.1)
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
print('\n\n"smf_comparison_field.py"  terminated successfully.\n')
#
#
#
###################     PROGRAM END     ###################
