#Created on Sun Jul 12 17:31:48 2020
#
#####################  velocity_dispersion.py  #####################
#
###  This program is designed to measure the size of a galaxy cluster. From the observable redshift (z), it calculates the velocity dispersion and, applying the Virial Theorem, calculates the total Kinetic Energy of the cluster. Equating this to the Potential Energy of the cluster, we can define and measure a radius (and a mass contained within a sphere of this radius), where the average density is 200x greater than the critical density for a spatially flat, matter-lambda dominated universe. The key results are this scale radius, r_200, and the mass contained within a sphere of this radius, M_200.
#
#
#
### NOTE: To find a section, search "SECTION (*)" in all caps, where * is the section number you want. To find fields which require user input, search "MAY NEED TO EDIT" in all caps. Some of these fields may direct you to manually enter input into a sub-program. 
#
### Section summary:
#
### PROGRAM START
#
### (0)    import modules
### (0.1)  define FUNCTIONS; 
### (0.2)  define CONSTANTS & USER INPUTS
### (1)    convert REDSHIFT TO RECESSIONAL VELOCITY
### (1.1)   prepare SUMMARY TABLES to check done for all galaxies in MEMBERSHIP sample
### (2)    compute VELOCITY DISPERSION
### (3)    apply the VIRIAL THEOREM
### (3.1)   LHS: total KINETIC ENERGY of cluster
### (3.2)   RHS: total POTENTIAL ENERGY of cluster
### (3.3)   summary table of calculated MASS_200, RADIUS_200 for all 6 clusters
#
### PROGRAM END
#
#
#
### NOTE: if DIAG_FLAG_1 IS NOT TURNED ON, THIS PROGRAM WILL NOT COMPUTE LIMITING MASSES. This is a pretty big design flaw on my part. 
#
#
#
###################     PROGRAM START     ###################
#
#
time_flag = 1     # 0= all timers off;   1=on, time entire program;    2=off, allow individual sections to be timed
#
if time_flag == 1:
    if project_time_flag == 1:
        pass
    else:
        start_time = time.time()
#
#
## SECTION (0): Modules, FLAGS
#
## import modules & definitions
#
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
import time
#
## FLAGS
#
diag_flag = 1
#
## SECTION (0.1): Modules, Definitions, FLAGS
#
## DEFINITIONS 
#
## define a function to convert a redshift to a recessional velocity
def recessional_velocity(z):
    v = c * ( ((1+z)**2 - 1) / ((1+z)**2 + 1) )
    return v
#
## define a function to compute the "Hubble Parameter" by solving the Friedmann Equation for a universe containing only matter and Lambda-dark energy.
def friedmann_equation(z_clu):
    H2 = (H_0**2) * ( (omega_m  * (1+z_clu)**3) + omega_lambda )
    return H2
#
## define a function to convert metres to parsecs
def m_to_kpc(r_m):
    r_kpc = r_m / 3.0857e19                # conversion factor from wikipedia 
    return r_kpc
#
## define a function to convert kilograms to solar masses
def kg_to_log10_Msol(m_kg):
    m_Msol = m_kg / 1988500e24             # conversion factor from NASA sun fact sheet 
    log10_Msol = np.log10(m_Msol)
    return log10_Msol
#
## define a function to convert kilograms to solar masses
def log10_Msol_to_kg(m_log10_Msol):
    m_kg = np.float128
    m_kg = 10**m_log10_Msol              # get mass in solar units, not log_10()
    m_kg = m_kg * 1988500e24             # conversion factor from NASA sun fact sheet: 1 M_sol = 1,988,500e24 kg 
    return m_kg
#
#
## SECTION (0.2): Modules, Definitions, FLAGS
#
## CONSTANTS & USER INPUTS
## MAY NEED TO EDIT
#
H_0 = 70                          # units [s^-1]; NO REFERENCE - THIS WILL NEED TO BE UPDATED: MAY NEED TO EDIT
G = 6.67430e-11                   # gravitational constant; units [m^3*(kg*s^2)^-1]; physics.nist.gov
c = 2.99792458e+8                 # speed of light, units [m/s]; from physics.nist.gov
omega_m = 0.3                     # energy density of matter (***AT THE PRESENT EPOCH???***)
omega_lambda = 0.7                # energy density of Lambda (both of these omegas are dimensionless)
#
##
#
## SECTION (1): translate REDSHIFT to RECESSIONAL VELOCITY
#
## start by adding a column to master_cat store recessional velocities, if it doesn't already exist (add try/except for repeatablility)
#
try: 
    F = Column([-99]*len(master_cat),name='v_rec',dtype='float128')   # add a column in which to store magnitudes
    master_cat.add_column(F)
    print("This is the first time you run this program today, isn't it ya ghassan?")
except:
    pass
#
## convert redshifts to recessional velocities for all cluster members
#
counting_array = np.array([0]*3)
cluster_counting_array = np.array([[0]*6]*2)      # row1=spec;  row2=phot
cluster_velocities = [ [], [], [], [], [], [] ]
cluster_mass = [ [], [], [], [], [], [] ]
lost_below_lim_mass = np.array([[0]*6]*2)         # row1=spec;  row2=phot
#
for counter in range(len(master_cat)):
    if master_cat['member'][counter] == 0:     # select only cluster members to perform this calculation
        counting_array[0]+=1                   # member counter
        for ii in range(len(cluster_velocities)):# track velocity dispersion by cluster
            if master_cat['cluster'][counter] == (ii+1):
                if master_cat['lmass'][counter] >= limiting_mass[ii]:
                    if master_cat['sub'][counter] == 1:    # use z_spec where available
                        counting_array[1]+=1    # spec-only counter
                        master_cat['v_rec'][counter] = recessional_velocity(master_cat['z_spec'][counter])
                        cluster_counting_array[0][ii]+=1   # cluster counter
                        cluster_velocities[ii].append(master_cat['v_rec'][counter])
                        cluster_mass[ii].append(master_cat['lmass'][counter])
                    elif master_cat['sub'][counter] == 2:  # otherwise use z_phot
                        counting_array[2]+=1    # phot-only counter
                        #master_cat['v_rec'][counter] = recessional_velocity(master_cat['z_peak'][counter])
                        #cluster_velocities[ii].append(master_cat['v_rec'][counter])
                        #cluster_mass[ii].append(master_cat['lmass'][counter])  
                        #cluster_counting_array[1][ii]+=1   # cluster counter
                elif master_cat['lmass'][counter] < limiting_mass[ii]:
                    if master_cat['sub'][counter] == 1:
                        lost_below_lim_mass[0][ii]+=1
                    elif master_cat['sub'][counter] == 2:
                        lost_below_lim_mass[1][ii]+=1
                    
#
#
#
if (diag_flag == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        ## Summarize initial data stats in table
        vel_names = Column(['FULL PHOT (Parent)','Phot only (Parent)','Spec+Phot (Parent)','Phot only','Spec+Phot','Lost below lim. mass','SUM'],name='Property')
        col_names = cluster_names
        vel0 = Column([np.sum([mem_spec[0],mem_spec[1],mem_phot[0],mem_phot[1]]),np.sum([mem_phot[0],mem_phot[1]]),np.sum([mem_spec[0],mem_spec[1]]),counting_array[2],counting_array[1],np.sum(lost_below_lim_mass),np.sum([counting_array[1],counting_array[2],np.sum(lost_below_lim_mass)])],name='Total')  # total column
        vel_stats = Table([vel_names,vel0])
        for ii in range(len(spec_only)):
            vel_col = Column([np.sum([mem_spec[0][ii],mem_spec[1][ii],mem_phot[0][ii],mem_phot[1][ii]]),np.sum([mem_phot[0][ii],mem_phot[1][ii]]),np.sum([mem_spec[0][ii],mem_spec[1][ii]]),cluster_counting_array[1][ii],cluster_counting_array[0][ii],np.sum([lost_below_lim_mass[0][ii],lost_below_lim_mass[1][ii]]),np.sum([cluster_counting_array[1][ii],cluster_counting_array[0][ii],lost_below_lim_mass[0][ii],lost_below_lim_mass[1][ii]])],name=col_names[ii])  # add columns to table one cluster at a time
            vel_stats.add_column(vel_col)
        #
        print('\nSummary Table - Redshift-to-velocity conversion:\n%s'%vel_stats,'\nNOTE: Redshift were converted to recessional velocities for the SPEC SUBSAMPLE ONLY.')
        #
#
#
#
#
## SECTION (2): compute VELOCITY DISPERSION
#
## calculate the velocity dispersion in the z-direction (sigma_z), which is the standard deviation of the recessional velocities (i.e. the velocities in the z-direction); ASSUMPTION: velocity is ISOTROPIC, so velocity dispersion squared is sum of squared components (x,y,z). we calculate the z-component from the redshift (recessional velocity), square it, and multiply by 3 to get the squared total velocity dispersion.
#
vel_disp_mean = np.array([0]*6,dtype='float128')
vel_disp_z = np.array([0]*6,dtype='float128')
vel_disp_sq = np.array([0]*6,dtype='float128')
for ii in range(len(vel_disp)):
    vel_disp_mean[ii] = np.mean(cluster_velocities[ii])
    vel_disp_z[ii] = np.std(cluster_velocities[ii])
    vel_disp_sq[ii] = 3*(np.std(cluster_velocities[ii])**2)
#
if (diag_flag == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        print('\nMean recessional velocity (km/s): \n%s'%(vel_disp_mean/1e3),'\nRecessional velocity dispersion (z-direction) by cluster (km/s): \n%s'%(vel_disp_z/1e3),'\nTOTAL velocity dispersion by cluster (km/s): \n%s'%(np.sqrt(vel_disp_sq)/1e3),'\nVelocity dispersion squared (total) by cluster (km/s): \n%s'%(vel_disp_sq/1e3))
#
#
#
#
## SECTION (3): apply the VIRIAL THEOREM for a closed system w/ constant angular momentum: <KE> = - (1/2)*<PE>; where <KE> = 0.5 * M_star * velocity_dispersion**2 is the LHS; M_star is the total stellar mass of all objects in the spectroscopic subsample; the velocity dispersion is the standard deviation of the velocities of the spec subsample; the RHS is <PE> = G*M_200/r_200, where M_200 = (4/3)*pi*(r_200**3)*rho_200; rho_200 = 200 * ( (3 * H**2(t)) / (8*pi*G) ); H**2(t) = H_o**2 { Omega_m*(1+z)**3 + Omega_lambda } for a flat universe with matter + lambda, no radiation. 
#
## SECTION (3.1): LHS: compute the Kinetic Energy of each cluster
#
## first calculate the total stellar mass in each cluster; then apply KE = 0.5 * (M_star total) * vel_disp_sq
#
cluster_stellar_mass_total_Msol = np.array([0]*6,dtype='float128')
#
for ii in range(len(cluster_mass)): 
    cluster_stellar_mass_total_Msol[ii] = np.sum(cluster_mass[ii])
#
## convert masses from [log10(solar mass units)] to [kg]
cluster_stellar_mass_total_kg = np.empty_like(cluster_stellar_mass_total_Msol,dtype='float128')
cluster_stellar_mass_total_kg = log10_Msol_to_kg(cluster_stellar_mass_total_Msol)
#
#
kinetic_energy = 0
for ii in range(len(cluster_stellar_mass_total)):
    kinetic_energy+=(0.5*cluster_stellar_mass_total_kg[ii]*vel_disp_sq[ii])
#
#
if (diag_flag == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        print('\nTotal STELLAR Mass by cluster (units: [M_sol]:\n%s'%cluster_stellar_mass_total_Msol)
#
## SECTION (3.2): RHS: compute the POTENTIAL ENERGY --> really, compute r_200 & M_200: the radius within which the average density is 200 times the critical density of a spatially flat univervse, and the mass enclosed within a sphere of that radius
#
## r_200: for a derivation, see Notes; first calculate the Hubble Parameter ("H2") for each cluster 
#
H2 = np.array([0]*6,dtype='float128')
for ii in range(len(z_cluster)):
    H2[ii] = friedmann_equation(z_cluster[ii])
#
r_200_m = np.array([0]*6,dtype='float128')
for ii in range(len(r_200_m)):
    r_200_m[ii] = np.sqrt((1/100)*(1/H2[ii])*cluster_stellar_mass_total_kg[ii]*vel_disp_sq[ii])    # dimensionally correct: units of 
#
## convert r_200 from metres to kiloparsecs
r_200_kpc = np.array([0]*6,dtype='float128')
r_200_kpc = m_to_kpc(r_200_m)
## M_200: for a derivation, see Notes; 
#
M_200_kg = np.array([0]*6,dtype='float128')
for ii in range(len(M_200_kg)):
    M_200_kg[ii] = (100)*(H2[ii])*(r_200_m[ii]**3)/G        # dimensionally correct
#
## convert M_200 from kg to log10(M_sol)
M_200_Msol = np.array([0]*6,dtype='float128')
M_200_Msol = kg_to_log10_Msol(M_200_kg)
#
if (diag_flag == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        print('\nHubble factor by cluser:\n%s'%H2,'\nr_200 by cluster (kpc):\n%s'%r_200_kpc,'\nM_200 by cluster (log10(M_sol)):\n%s'%M_200_Msol)
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
        print('Program "data_mass_completeness*.py" took: %s seconds to run.\n\n' % (time.time() - start_time))
#
#
#
#
print('\n\n"velocity_dispersion.py"  terminated successfully.\n')
#
#                        
###################     PROGRAM END     ###################