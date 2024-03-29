#Created on Thu Aug 20 20:25:44 2020
#
#####################  field_normalization.py  #####################
#
###  This program calculates the volume subtended by the UltraVISTA/COSMOS (UVC) survey catalogue, (1.62 degrees sq., see catalogue Muzzin et al. 2013a), and the volume subtended by the 6 HFF parallel fields, and the 6 HFF cluster fields (less the volume subtended by the cluster in each frame). It then scales the volume of the UltraVISTA survey to the volume of the HFF, and applies this scaling factor to the UVC SMF, effectively scaling it to the same "size" as the HFF field SMF
#
#
### Section summary:
#
### PROGRAM START
#
### (0)    import MODULES, set FLAGS, define FUNCTIONS;
### (1)    record volume of ENTIRE SKY within redshift range of clusters
### (2)    determine volume of UlatraVISTA catalogue (UVC)
### (3)    determine volume of HFF Parallel Fields
### (4)    determine volume of field sample from HFF Cluster Fields
### (5)    determine the SCALING FACTOR to apply to UVC s.t. it matches total vol of HFF
#
### PROGRAM END
#
#
###################     PROGRAM START     ###################
#

## SECTION (0): modules, flags, constants, definitions
#
## Modules & Time flag
#
import numpy as np
import time
#
#
time_flag = 1     # 0= all timers off;   1=on, time entire program;    2=off, allow individual sections to be timed
#
## FLAGS
#
diag_flag_1 = 1                        # print proper distances and volumes subtended by each survey
#
## Constants & conversion factors
#
# areas of survey images
area_sky = 41252.96125                          # in units of degrees squared; REF: WIKI (https://en.wikipedia.org/wiki/Square_degree)
area_uvc = 1.62                                  # in units of degrees squared
area_cluster = np.array([14.1,12.5,15.4,15.1,14.6,18.2])  # [M0416,M1149,M0717,A370,A1063,A2744]
area_parallel = np.array([11.9,14.3,13.0,11.9,12.2,11.9])  # in units of arcmin squared
#
# constants/conversion factors
##
#H_0 = 70                      # units [kms^-1Mpc^-1]; Freedman et al 2019, APJ
#G = 6.67430e-11*((1/3.0857e+16)**3)*((1.9885e30/1)) # gravitational constant; units [m^3*k^-1g*s^-2]*[Mpc/m]^3*[kg/Msol] = [Mpc^3*(Msol*s^2)^-1]; physics.nist.gov
#c = 2.99792458e+8               # speed of light, units [m/s]; from physics.nist.gov
#omega_m = 0.3                   # energy density of matter (***AT THE PRESENT EPOCH???***)
#omega_lambda = 0.7              # energy density of Lambda (both of these omegas are dimensionless)
#
#D_H = (c/1e3) / H_0             # Hubble distance, in Mpc; see Hogg 2000

#
## DEFINITIONS
#
## define a function to convert sqaure arcmin to square degrees; conversion factor (3600 sq arcmin = 1 sq degree)
def arcmin_sq_to_deg_sq(arcmin_sq):
    deg_sq = arcmin_sq / 3600
    return deg_sq
#
## define a function to convert cubic Gpc to cubic Mpc (Gps^3 to Mpc^3)
def Gpc_to_Mpc_cube(Gpc3):
    Mpc3 = 1e9 * Gpc3
    return Mpc3
#
#
#
#
#
#
## compute redshift range of galaxies in cluster sample, i.e. the "height" of the volume enclosed (whose geometry is that of a pyramid, with the observer at the vertex).
# ## for z_cutoff = [0.012,0.055]
# if z_field_bounds_flag == 1:
#     lower_bound = (min(z_cluster) - z_cutoff[1]) / (1 + z_cutoff[1])  # == ~0.240
#     upper_bound = (max(z_cluster) + z_cutoff_field[1]) / (1 - z_cutoff_field[1])  # ==
#     z_field_bounds = [lower_bound, upper_bound]
# #
## define the lower_bound & upper_bound for each cluster
lower_bound_cluster = np.array([0.]*6)
upper_bound_cluster = np.array([0.]*6)
for cluster in range(len(z_cluster)):
    lower_bound_cluster[cluster] = (z_cluster[cluster] - z_cutoff_field[1]) / (1 + z_cutoff_field[1])
    upper_bound_cluster[cluster] = (z_cluster[cluster] + z_cutoff_field[1]) / (1 - z_cutoff_field[1])
    #
print('lower_bound_cluster: %s'%lower_bound_cluster )
print('upper_bound_cluster: %s'%upper_bound_cluster )
#
## RESULT
## z_cutoff_field[2] == 0.14
# lower_bound_cluster: [0.2245614  0.35350877 0.35526316 0.20614035 0.18245614 0.14736842]
# upper_bound_cluster: [0.62325581 0.79418605 0.79651163 0.59883721 0.56744186 0.52093023]
#
## z_cutoff_field[1] == 0.15
# lower_bound_cluster: [0.21391304 0.34173913 0.34347826 0.19565217 0.17217391 0.1373913 ]
# upper_bound_cluster: [0.64235294 0.81529412 0.81764706 0.61764706 0.58588235 0.53882353]

## SECTION (1): determine volume subtended by ENTIRE SKY over redshift range [lower_bound, upper_bound]
## NOTE: comoving volumes calculated using Ned Wright's Java Script Cosmological Calculatror (http://www.astro.ucla.edu/~wright/CosmoCalc.html); *** cite  Wright (2006, PASP, 118, 1711) ***
#
# Parameters: omega_matter = 0.3, omega_lambda = 0.7, H_0 = 70 km/s/Mpc
## COMOVING VOLUMES
if z_field_bounds_flag == 0:
    vol_sky_lower = 4.292            # volume subtended of entire sky from z=0  to z=lower_bound = 0.25; UNITS = [Gpc^3]
    vol_sky_upper = 77.829           # volume subtended of entire sky from z=0  to z_upper_bound = 0.75
elif z_field_bounds_flag == 1 and z_cutoff_field[1] == 0.125:
    vol_sky_lower = 4.292            # volume subtended of entire sky from z=0  to z=lower_bound; UNITS = [Gpc^3]    <-- z_lower_bound= 0.240
    vol_sky_upper = 85.793
elif z_field_bounds_flag == 1 and z_cutoff_field[1] == 0.14:
    vol_sky_lower = 3.826            # volume subtended of entire sky from z=0  to z=lower_bound; UNITS = [Gpc^3]    <-- z_lower_bound= 0.240
    vol_sky_upper = 89.793           # volume subtended of entire sky from z=0  to z=upper_bound; UNITS = [Gpc^3]    <-- for z_cutoff_field[1] == 0.14
elif z_field_bounds_flag == 1 and z_cutoff_field[1] == 0.15:
    vol_sky_lower = 3.826            # volume subtended of entire sky from z=0  to z=lower_bound; UNITS = [Gpc^3]    <-- z_lower_bound= 0.240
    vol_sky_upper = 95.494           # volume subtended of entire sky from z=0  to z=upper_bound; UNITS = [Gpc^3]    <-- for z_cutoff_field[1] == 0.15

#
vol_sky = vol_sky_upper - vol_sky_lower      # volume subtended of entire sky over redshift range of clusters; UNITS = [Gpc^3]
vol_sky = Gpc_to_Mpc_cube(vol_sky)      # UNITS = [Mpc^3]
#
## do the same but for the volume occupied by each cluster at redshift
## NOTE: the lower bound is fixedby the lower redshift bound of the clusters, and so is the same for all field cuts; it is also the floor redshift cut for all low-redshift clusters (i.e. the the lower_bound based on del_z calculation is below z ~.240, volume corresponding to z=.240 is adopted)
if z_field_bounds_flag == 0 and z_cutoff_field[1] == 0.125:
    vol_cluster_lower = np.array([4.292,12.200,12.347,4.292,4.292,4.292])        # volume subtended of entire sky from z=0  to z=lower_bound_cluster; UNITS = [Gpc^3]
    vol_cluster_upper = np.array([47.480,77.829,77.829,42.763,37.063,29.376])       # volume subtended of entire sky from z=0  to z=upper_bound_cluster; UNITS = [Gpc^3
elif z_field_bounds_flag == 0 and z_cutoff_field[1] == 0.14:
    vol_cluster_lower = np.array([4.292,11.200,11.347,4.292,4.292,4.292])        # volume subtended of entire sky from z=0  to z=lower_bound_cluster; UNITS = [Gpc^3]
    vol_cluster_upper = np.array([49.480,77.829,77.829,44.763,39.063,31.376])       # volume subtended of entire sky from z=0  to z=upper_bound_cluster; UNITS = [Gpc^3]
elif z_field_bounds_flag == 0 and z_cutoff_field[1] == 0.15:
    vol_cluster_lower = np.array([4.292,10.212,10.355,4.292,4.292,4.292])        # volume subtended of entire sky from z=0  to z=lower_bound_cluster; UNITS = [Gpc^3]
    vol_cluster_upper = np.array([53.340,77.829,77.829,48.376,42.362,34.225])       # volume subtended of entire sky from z=0  to z=upper_bound_cluster; UNITS = [Gpc^3]
elif z_field_bounds_flag == 1 and z_cutoff_field[1] == 0.14:
    vol_cluster_lower = np.array([3.826,11.200,11.347,3.826,3.826,3.826])        # volume subtended of entire sky from z=0  to z=lower_bound_cluster; UNITS = [Gpc^3]
    vol_cluster_upper = np.array([49.480,89.174,89.793,44.763,39.063,31.376])       # volume subtended of entire sky from z=0  to z=upper_bound_cluster; UNITS = [Gpc^3]
elif z_field_bounds_flag == 1 and z_cutoff_field[1] == 0.15:
    vol_cluster_lower = np.array([3.826,10.212,10.355,3.826,3.826,3.826])        # volume subtended of entire sky from z=0  to z=lower_bound_cluster; UNITS = [Gpc^3]
    vol_cluster_upper = np.array([53.340,94.851,95.494,48.376,42.362,34.225])       # volume subtended of entire sky from z=0  to z=upper_bound_cluster; UNITS = [Gpc^3]
#
#

vol_cluster_sky = vol_cluster_upper - vol_cluster_lower     # volume subtended by all 6 clusters, in units of [Gpc^3]
vol_cluster_sky = Gpc_to_Mpc_cube(vol_cluster_sky)                  # units [Mpc^3]
#
#
#
if (diag_flag_1 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        print('\nVolume of ENTIRE SKY from %s'%np.round(lower_bound,decimals=2)+' < z < %s'%np.round(upper_bound,decimals=2)+': %s'%np.round(vol_sky,decimals=4)+' (Mpc^3)')
        print('Volume CLUSTERS (ENTIRE SKY) from %s'%np.round(lower_bound_cluster,decimals=2)+' < z < %s'%np.round(upper_bound_cluster,decimals=2)+': %s'%np.round(vol_cluster_sky,decimals=4)+' (Mpc^3)')
#
#
#
#
## SECTION (2): determine volume of UlatraVISTA catalogue (UVC)
# determine what fraction of the entire sky was observed by UVC. apply this fraction to the volume subtended by the entire sky to determine the volume subtended by the UVC survey
#
frac_sky_uvc = area_uvc / area_sky         # dimensionless
#
vol_UVC = vol_sky * frac_sky_uvc #* 1e9          # UNITS = [Mpc^3]
#
#
if (diag_flag_1 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        print('\nVolume of UVC from %s'%np.round(lower_bound,decimals=3)+' < z < %s'%np.round(upper_bound,decimals=3)+': %s'%np.round(vol_UVC,decimals=3)+' (Mpc^3)')
#
#
## SECTION (3): determine volume of HFF PARALLEL FIELDS
# determine what fraction of the entire sky was observed by the PARALLEL fields. apply this fraction to the volume subtended by the entire sky to determine the volume subtended by the HFF Parallel fields.
#
area_parallel_total = area_parallel                              # total area of parallel fields
area_parallel_total = arcmin_sq_to_deg_sq(area_parallel_total)           # convert to sq. degrees
frac_sky_parallel = area_parallel_total / area_sky                       # fraction of entire sky covered by parallel fields
#
vol_parallel_field = np.array([0.]*6)
vol_parallel_cluster = np.array([0.]*6)
for cluster in range(len(frac_sky_parallel)):
    vol_parallel_field[cluster] = vol_sky * frac_sky_parallel[cluster]         # volume subtended by the entire field images, UNITS = [Mpc^3]
    vol_parallel_cluster[cluster] = vol_cluster_sky[cluster] * frac_sky_parallel[cluster]          # volume subtended by the clusters in the field images, UNITS = [Mpc^3]
#
vol_parallel = np.sum(vol_parallel_field) - np.sum(vol_parallel_cluster)            # UNITS = [Mpc^3]
#
if (diag_flag_1 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        print('\nVolume of HFF Parallel Fields from %s'%np.round(lower_bound,decimals=3)+' < z < %s'%np.round(upper_bound,decimals=3)+': %s'%np.round(vol_parallel,decimals=3)+' (Mpc^3)')
        print('Vol. Parallel IMAGES: %s'%np.round(np.sum(vol_parallel_field),decimals=3)+' [Mpc^3]')
        print('Vol. CLUTSERS in Parallel images: %s'%np.round(np.sum(vol_parallel_cluster),decimals=3)+' [Mpc^3]')
        print('Note: vol_parallel = vol_parallel_field - vol_parallel_cluster')
#
#
#
#
## SECTION (4): determine volume of HFF CLUSTER FIELDS
# determine what fraction of the entire sky was observed by the CLUSTER fields. apply this fraction to the volume subtended by the entire sky to determine the volume subtended by the HFF Parallel fields.
#
area_cluster_image_total = area_cluster                                    # total area of cluster fields, including volume occupied by cluster
area_cluster_image_total = arcmin_sq_to_deg_sq(area_cluster_image_total)           # convert to sq. degrees
#
frac_sky_cluster_image = area_cluster_image_total / area_sky                       # fraction of entire sky covered by cluster images
#
vol_cluster_image = np.array([0.]*6)
vol_cluster = np.array([0.]*6)
#
for cluster in range(len(frac_sky_cluster_image)):
    vol_cluster_image[cluster] = vol_sky * frac_sky_cluster_image[cluster]# * 1e9          # UNITS = [Mpc^3];
    #
    # now do the same to determine the volume occupied in each image by the clusters themselves.
    vol_cluster[cluster] = (vol_cluster_sky[cluster] * frac_sky_cluster_image[cluster]) #* 1e9          # UNITS = [Mpc^3];
#
#
if (diag_flag_1 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        print('\nVolume of HFF Cluster IMAGES (total) from %s'%np.round(lower_bound,decimals=3)+' < z < %s'%np.round(upper_bound,decimals=3)+': %s'%np.round(np.sum(vol_cluster_image),decimals=3)+' (Mpc^3)')
        print('Volume of HFF CLUSTERS (cluster volume inside "cluster" images) from %s'%np.round(lower_bound_cluster,decimals=3)+' < z < %s'%np.round(upper_bound_cluster,decimals=3)+': %s'%np.round(np.sum(vol_cluster),decimals=3)+' (Mpc^3)')
        print('\nVolume of HFF Cluster FIELD (total volume of image less volume of clusters): %s'%np.round(np.sum(vol_cluster_image),decimals=3)+' - %s'%np.round(np.sum(vol_cluster),decimals=3)+' = %s'%np.round((np.sum(vol_cluster_image) - np.sum(vol_cluster)),decimals=3)+' (Mpc^3)')
#
#
#
#
## SECTION (5): determine SCALING FACTOR which will scale UVC volume to match the TOTAL HFF volume (i.e. vol_parallel + vol_cluster_image - vol_cluster)
if cluster_field_inclusion_flag == 0:
    vol_HFF = vol_parallel
elif cluster_field_inclusion_flag == 1:
    vol_HFF = vol_parallel + (np.sum(vol_cluster_image) - np.sum(vol_cluster))           # units: [Mpc^3]
#
## SCALING FACTOR
scale_factor_uvc_HFF = vol_HFF / vol_UVC
#
#
if (diag_flag_1 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        print('Volume of HFF TOTAL FIELD SURVEY: %s'%np.round(vol_HFF,decimals=3)+' (Mpc^3)')
        print('\nSCALE FACTOR: %s'%scale_factor_uvc_HFF)
#
#
#
#
## TIME_FLAG END
#
if time_flag == 1:
    if project_time_flag == 1:
        pass
    else:
        print('Program "field_normalization.py" took: %s seconds to run.\n\n' % (time.time() - start_time))
#
#
#
#
print('\n\n"field_normalization.py"  terminated successfully.\n')
#
#
###################     PROGRAM END     ###################
