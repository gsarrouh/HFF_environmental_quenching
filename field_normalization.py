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
area_cluster = [14.1,12.5,15.4,15.1,14.6,18.2]  # [M0416,M1149,M0717,A370,A1063,A2744]
area_parallel = [11.9,14.3,13.0,11.9,12.2,11.9]  # in units of arcmin squared
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
#
#
#
#
#
#
## compute redshift range of galaxies in cluster sample, i.e. the "height" of the volume enclosed (whose geometry is that of a pyramid, with the observer at the vertex).
lower_bound = (min(z_cluster) - z_cutoff[1]) / (1 + z_cutoff[1])  # == ~0.234
upper_bound = (max(z_cluster) + z_cutoff[1]) / (1 - z_cutoff[1])  # == ~0.644
z_field_bounds = [lower_bound, upper_bound]
#
## define the lower_bound & upper_bound for a cluster at z=0.4, z_cutoff=0.06
lower_bound_cluster = (0.4 - 0.06) / (1 + 0.06) # == ~0.321
upper_bound_cluster = (0.4 + 0.06) / (1 - 0.06) # == ~0.489
#
#
## SECTION (1): determine volume subtended by ENTIRE SKY over redshift range [lower_bound, upper_bound]
## NOTE: comoving volumes calculated using Ned Wright's Java Script Cosmological Calculatror (http://www.astro.ucla.edu/~wright/CosmoCalc.html); *** cite  Wright (2006, PASP, 118, 1711) ***
#
# Parameters: omega_matter = 0.3, omega_lambda = 0.7, H_0 = 70 km/s/Mpc
## COMOVING VOLUMES
vol_sky_lower = 3.561            # volume subtended of entire sky from z=0  to z=lower_bound; UNITS = [Mpc^3]
vol_sky_upper = 53.598           # volume subtended of entire sky from z=0  to z=upper_bound; UNITS = [Mpc^3]
#
vol_sky = vol_sky_upper - vol_sky_lower      # volume subtended of entire sky over redshift range of clusters
#
## do the same but for the volume occupied by a cluster at redshift z=0.4, del_z = 0.06
vol_cluster_lower = 8.583        # volume subtended of entire sky from z=0  to z=lower_bound_clustser; UNITS = [Mpc^3]
vol_cluster_upper = 26.677       # volume subtended of entire sky from z=0  to z=upper_bound_clustser; UNITS = [Mpc^3]
#
vol_cluster_sky = vol_cluster_upper - vol_cluster_lower
#
#
if (diag_flag_1 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        print('Volume of ENTIRE SKY from %s'%np.round(lower_bound,decimals=2)+' < z < %s'%np.round(upper_bound,decimals=2)+': %s'%np.round(vol_sky,decimals=2)+' (Mpc^3)')
        print('Volume a cluster (ENTIRE SKY) from %s'%np.round(lower_bound_cluster,decimals=2)+' < z < %s'%np.round(upper_bound_cluster,decimals=2)+': %s'%np.round(vol_cluster_sky,decimals=2)+' (Mpc^3)')
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
        print('\nVolume of UVC from %s'%np.round(lower_bound,decimals=3)+' < z < %s'%np.round(upper_bound,decimals=3)+': %s'%np.round(vol_UVC,decimals=3)+' (kpc^3)')
#
#
## SECTION (3): determine volume of HFF PARALLEL FIELDS
# determine what fraction of the entire sky was observed by the PARALLEL fields. apply this fraction to the volume subtended by the entire sky to determine the volume subtended by the HFF Parallel fields.
#
area_parallel_total = np.sum(area_parallel)                              # total area of parallel fields
area_parallel_total = arcmin_sq_to_deg_sq(area_parallel_total)           # convert to sq. degrees
frac_sky_parallel = area_parallel_total / area_sky                       # fraction of entire sky covered by parallel fields [Mpc^3]
#
vol_parallel = vol_sky * frac_sky_parallel# * 1e9          # UNITS = [kpc^3]
#
#
if (diag_flag_1 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        print('\nVolume of HFF Parallel Fields from %s'%np.round(lower_bound,decimals=3)+' < z < %s'%np.round(upper_bound,decimals=3)+': %s'%np.round(vol_parallel,decimals=3)+' (kpc^3)')
#
#
#
#
## SECTION (4): determine volume of HFF PARALLEL FIELDS
# determine what fraction of the entire sky was observed by the PARALLEL fields. apply this fraction to the volume subtended by the entire sky to determine the volume subtended by the HFF Parallel fields.
#
area_cluster_image_total = np.sum(area_cluster)                                    # total area of cluster fields, including volume occupied by cluster
area_cluster_image_total = arcmin_sq_to_deg_sq(area_cluster_image_total)           # convert to sq. degrees
#
frac_sky_cluster_image = area_cluster_image_total / area_sky                       # fraction of entire sky covered by cluster images [Mpc^3]
#
vol_cluster_image = vol_sky * frac_sky_cluster_image# * 1e9          # UNITS = [kpc^3];
#
## now do the same to determine the volume occupied in each image by the clusters themselves. to do this, assume all 6 clusters were at the redshift z=0.4, with cutoff +/- 0.06
vol_cluster = (vol_cluster_sky * frac_sky_cluster_image) #* 1e9          # UNITS = [kpc^3];
#
#
if (diag_flag_1 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        print('\nVolume of HFF Cluster IMAGES from %s'%np.round(lower_bound,decimals=3)+' < z < %s'%np.round(upper_bound,decimals=3)+': %s'%vol_cluster_image+' (Mpc^3)')
        print('Volume of HFF CLUSTERS from %s'%np.round(lower_bound_cluster,decimals=3)+' < z < %s'%np.round(upper_bound_cluster,decimals=3)+': %s'%vol_cluster+' (Mpc^3)')
        print('Volume of HFF Cluster FIELD: %s'%np.round(vol_cluster_image,decimals=3)+' - %s'%np.round(vol_cluster,decimals=3)+' = %s'%(vol_cluster_image - vol_cluster)+' (Mpc^3)')
#
#
#
#
## SECTION (5): determine SCALING FACTOR which will scale UVC volume to match the TOTAL HFF volume (i.e. vol_parallel + vol_cluster_image - vol_cluster)
vol_HFF = vol_parallel + vol_cluster_image #- vol_cluster           # units: [Mpc^3]
#
## SCALING FACTOR
scale_factor_uvc_HFF = vol_HFF / vol_UVC
#
#
if (diag_flag_1 == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        print('Volume of HFF TOTAL FIELD SURVEY: %s'%vol_HFF+' (Mpc^3)')
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
