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
### (1)    determine volume of UltraVISTA catalogue 
### (2)    determine volume of HFF Parallel Fields
### (3)    determine volume of field sample from HFF Cluster Fields
### (4)    
#
### PROGRAM END
#
#
### NOTE: if DIAG_FLAG_1 IS NOT TURNED ON, THIS PROGRAM WILL NOT COMPUTE LIMITING MASSES. This is a pretty big design flaw on my part. 
#
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
#
## Constants & conversion factors
#
# areas of survey images
area_uvc = 1.62                                  # in units of degrees squared
area_clusters = [14.1,12.5,15.4,15.1,14.6,18.2]  # [M0416,M1149,M0717,A370,A1063,A2744]
area_parallel = [11.9,14.3,13.0,11.9,12.2,11.9]  # in units of arcmin squared
#
# constants/conversion factors
#
H_0 = 67.4                      # units [kms^-1Mpc^-1]; Freedman et al 2019, APJ
G = 6.67430e-11*((1/3.0857e+16)**3)*((1.9885e30/1)) # gravitational constant; units [m^3*k^-1g*s^-2]*[Mpc/m]^3*[kg/Msol] = [Mpc^3*(Msol*s^2)^-1]; physics.nist.gov
c = 2.99792458e+8               # speed of light, units [m/s]; from physics.nist.gov
omega_m = 0.3                   # energy density of matter (***AT THE PRESENT EPOCH???***)
omega_lambda = 0.7              # energy density of Lambda (both of these omegas are dimensionless)
#
D_H = (c/1e3) / H_0             # Hubble distance, in Mpc; see Hogg 2000

#
## DEFINITIONS
#
## define a function to convert redshift "z" to physical distance in Mpc "D"
def z_to_Mpc(z):
    # insert function D(z) = ???
    return D
    
#
## define a function for the small angle approximation; wh define a function ere theta is in radians and D is in Mpc
def small_angle_approx(theta,D):
    d = theta * D
    return d
#
## define a function 
    
#
## define a function to compute the volume of a right square pyramid, where all quantities are in units of Mpc; [vol] = Mpc^3
def vol_of_pyramid(base_side_length,height):          # length of side of square base of pyramid
    vol = (base_side_length**2)*height/3
    return vol
#
#
#
#
## compute redshift range of galaxies in cluster sample, i.e. the "height" of the volume enclosed (whose geometry is that of a pyramid, with the observer at the vertex). 
lower_bound = (min(z_cluster) - z_cutoff[1]) / (1 + z_cutoff[1])
upper_bound = (max(z_cluster) + z_cutoff[1]) / (1 - z_cutoff[1])
z_field_bounds = [lower_bound, upper_bound]
#
#
## SECTION (1): determine volume subtended by UVC survey 
#
# The volume subtended by a survey from redshift range [0,z], for a square field of view, takes on the geometry of a pyramid with a square base. The angular area of the base is known, but needs to be converted into a physical distance using the small angle approximation (i.e. sin(theta) = theta = d/D, where d is (half) the side length of the square base, and D is the redshift converted into a physical distance (i.e. the "height" of the pyramid)). 
# To compute the volume subtended between two redshifts [z_a, z_b], calculate the volume of [0,z_b] and [0,z_a] and subtract the latter from the former
#
## Calculate the angular dimensions of the (square) field of view: area = theta^2
#
theta_uvc = np.sqrt(area_uvc)*(np.pi/180)             # in radians
#
## First, convert redshift z to Mpc. Next, determine side length of base in Mpc. Volume calculated. 
#
## height of pyramid
lower_bound_height = z_to_Mpc(lower_bound)
upper_bound_height = z_to_Mpc(upper_bound)
## side length of pyramid
lower_bound_side_length = small_angle_approx(theta_uvc,lower_bound)*2
upper_bound_side_length = small_angle_approx(theta_uvc,upper_bound)*2
#
## Compute volume of a pyramid of height "lower bound"; the volume of a pyramid of height "upper bound"; then subtract the two to get the volume subtended by the survey
#
lower_bound_volume = vol_of_pyramid(lower_bound_side_length, lower_bound_height)
upper_bound_volume = vol_of_pyramid(upper_bound_side_length, upper_bound_height)
#
volume_uvc = upper_bound_volume - lower_bound_volume
#
#
#
## SECTION (2): determine volume subtended by parallel fields
#
# Do a similar thing as per above, but do it in a loop for each of the cluster fields, since the area of each field differs. You're after the total volume, in Mpc^3







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