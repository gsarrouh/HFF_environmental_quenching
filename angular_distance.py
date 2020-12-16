# Created on Fri Nov 20 14:17:33 2020
#
### angular_distance.py ###
#
### WHAT THIS PROGRAM DOES:
### This script computes the angular distance on the sky from the centre of each cluster to the position of each respective galaxy. This is done to refine selection of the field sample to exclude anything within ~0.5' of the cluster centre
#
#
## SECTION (0): modules, flags, definitions
#
## MODULES
import numpy as np
import astropy
from astropy.table import Table
from astropy.table import Column
import time
#
#
## FLAGS & USER-SET INPUTS
#
#
diag_flag_1 = 1         # various diagnostic outputs
diag_flag_2 = 0         # this is a plot flag, to produce a histogram of the angular distance distribution
#
## DEFINITIONS
#
## define a function to convert RA hours into degrees
def hrs_to_deg(hrs):
    deg = 15 * hrs
    return deg
#
## define a function to convert RA min into degrees
def min_to_deg(min):
    deg = (15/60) * min
    return deg
#
## define a function to convert RA sec into degrees
def sec_to_deg(sec):
    deg = (15/3600) * sec
    return deg
#
## define a function to convert DEC arcmin into degrees
def arcmin_to_deg(arcmin):
    deg = (1/60) * arcmin
    return deg
#
## define a function to convert DEC arcmin into degrees
def arcsec_to_deg(arcsec):
    deg = (1/3600) * arcsec
    return deg
#
## define a function to convert degrees to arcmin
def deg_to_arcmin(deg):
    arcmin = 60 * deg
    return arcmin
#
#
#
## SECTION (1): convert cluster coordinates to degree-decimal point format
## NOTE: all RA/DEC values retrieved from Shipley+18 Table 1
#
#
# cluster_names = ['M0416','M1149','M0717','A370','A1063','A2744']
RA_cluster_h = np.array([4,11,7,2,22,0])                          # RA: in hours, minutes, seconds, stored separately to facilitate unit conversions
RA_cluster_m = np.array([16,49,17,39,48,14])
RA_cluster_s = np.array([8.38,35.43,34.00,52.80,44.30,21.20])
DEC_cluster_d = np.array([-24,22,37,-1,-44,-30])                    # DEC: in degrees, arcmin, arcsec
DEC_cluster_m = np.array([-4,23,44,-34,-31,-23])
DEC_cluster_s = np.array([-20.80,44.63,49.00,-36.00,-48.40,-50.10])
#
## convert all RA values into degrees
RA_cluster_h = hrs_to_deg(RA_cluster_h)
RA_cluster_m = min_to_deg(RA_cluster_m)
RA_cluster_s = sec_to_deg(RA_cluster_s)
#
## convert all DEC values to degrees
DEC_cluster_m = arcmin_to_deg(DEC_cluster_m)
DEC_cluster_s = arcsec_to_deg(DEC_cluster_s)
#
## sum totals to get decimal-degree format for cluster centre RA & DEC
RA_cluster = np.sum([RA_cluster_h,RA_cluster_m,RA_cluster_s],axis=0)
DEC_cluster = np.sum([DEC_cluster_d,DEC_cluster_m,DEC_cluster_s],axis=0)
#
## update the user on the result of the calculation
if diag_flag_1:
    print('RA of clusters (in degrees): %s'%RA_cluster)
    print('DEC of clusters (in degrees): %s'%DEC_cluster)
#
#
#
## SECTION (2): compute angular distance
#
master_cat.add_column(Column([-99]*len(master_cat),name='ang_dist', dtype=np.float64))        # angular distance from the cluster centre, in arcmin
#
for counter in range(len(master_cat)):                                  # loop through catalogue one line at a time
    for cluster in range(len(RA_cluster)):                          # identify cluster
        if master_cat['cluster'][counter] == (cluster + 1):
            master_cat['ang_dist'][counter] = deg_to_arcmin( np.sqrt( (RA_cluster[cluster] - master_cat['ra'][counter])**2 + (DEC_cluster[cluster] - master_cat['dec'][counter])**2  ) )
#
#
#
#
## SECTION (3): optional diversion - plot a histogram of the angular distance distribution
#
if diag_flag_2 == 1:
    bins_ang_dist = np.arange(0.0,3.6,0.1)
    ang_dist_list =  [[] for x in range(len(z_cluster))]
    #
    for counter in range(len(master_cat)):
        for cluster in range(len(z_cluster)):
            if master_cat['cluster'][counter] == (cluster+1):
                if master_cat['z_peak'][counter] > z_field_bounds[0] and master_cat['z_peak'][counter] < z_field_bounds[1]:
                    ang_dist_list[cluster].append(master_cat['ang_dist'][counter])
    #
    nrows=2
    ncols=3
    #plot_titles = ['SF','Q']
    #
    fig, axs = plt.subplots(nrows,ncols,sharex=True,sharey=True,tight_layout=True)
    #fig.suptitle('abs(DEL_Z PHOT)',fontsize=12)
    #
    counter = 0
    #
    for jj in range(nrows):
        for kk in range(ncols):
            if counter < len(z_cluster):
            #
                ax = axs[jj][kk]
                #
                ax.hist(np.abs(ang_dist_list[counter]),bins=bins_ang_dist,color='deepskyblue',edgecolor='steelblue',alpha=0.7)
                ax.grid(axis='y', alpha=0.75)
                # ax.set_xlim([min(bins_phot),max(bins_phot)])
                # ax.title.set_text(cluster_names[kk])
                #ax.text(0.05,2500,cluster_names[counter],fontsize=9)
                counter+=1
                if jj == 0:
                    ax.title.set_text(cluster_names[kk])
                elif jj == 1:      # put xlabels on 1st row
                    ax.set_xlabel('Angular Distance (arcmin)',fontsize=10)
                    ax.title.set_text(cluster_names[kk+3])
                if kk == 0:
                    ax.set_ylabel('# count',fontsize=10)
                fig.suptitle('CLU Angular distance from cluster center (%s < '%np.round(z_field_bounds[0],decimals=2)+' z < %s )'%np.round(z_field_bounds[1],decimals=2))
    #####
    #
    #
    plt.show()
    #
##
print('\nProgram "angular_distance.py" terminated successfully.')
##
#
### PROGRAM END
