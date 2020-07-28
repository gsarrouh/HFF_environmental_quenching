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
from astropy.stats import sigma_clip
#
## FLAGS
#
diag_flag = 1
#
## SECTION (0.1): Modules, Definitions, FLAGS
#
## DEFINITIONS 
#
## define a function to convert a redshift to a recessional velocity in UNITS of m/s
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
def m_to_Mpc(r_m):
    r_Mpc = r_m / 3.0857e22                # conversion factor from wikipedia 
    return r_Mpc
#
## define a function to convert kilograms to solar masses
def kg_to_log10_Msol(m_kg):             
    m_log10_Msol = np.log10(m_kg / 1.9885e30) # conversion factor from NASA sun fact sheet 
    return m_log10_Msol
#
## define a function to convert log10(Msol) to log10(log10(Msol)), to deal with doing math on huge numbers
def log10_mSol_to_log10_log10_Msol(m_log10_Msol):
    m_log10_log10_Msol = np.log10(m_log10_Msol)        
    return m_log10_log10_Msol
#
## define a function to convert log10(log10(Msol)) back to log10(Msol)
def log10_log10_Msol_to_log10_mSol(log10_Msol):
    m_Msol = 10**log10_Msol     
    return m_Msol
#
## define a function to convert kilograms to solar masses
def log10_Msol_to_kg(m_log10_Msol):
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
H_0 = 67.4*(1/3.0857e19)        # units [kms^-1Mpc^-1]*[Mpc per km]; Freedman et al 2019, APJ
G = 6.67430e-11*((1/3.0857e+16)**3)*((1.9885e30/1)) # gravitational constant; units [m^3*k^-1g*s^-2]*[Mpc/m]^3*[kg/Msol] = [Mpc^3*(Msol*s^2)^-1]; physics.nist.gov
c = 2.99792458e+8               # speed of light, units [m/s]; from physics.nist.gov
omega_m = 0.3                   # energy density of matter (***AT THE PRESENT EPOCH???***)
omega_lambda = 0.7              # energy density of Lambda (both of these omegas are dimensionless)
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
vel_disp_z_m = np.array([0]*6,dtype='float128')
vel_disp_sq_m = np.array([0]*6,dtype='float128')
for ii in range(len(vel_disp_mean)):
    vel_disp_mean[ii] = np.mean(cluster_velocities[ii])
    vel_disp_z_m[ii] = = sigma_clip(cluster_velocities[ii], sigma=3, maxiters=5,cenfunc=mean)
    vel_disp_sq_m[ii] = 3*(vel_disp_z_m[ii]**2)
#
## convert squared velocity dispersion to units of [Mpc/s}^2
vel_disp_sq_Mpc = np.array([0]*6,dtype='float128')
#
for ii in range(len(vel_disp_sq_Mpc)):
    vel_disp_sq_Mpc[ii] = m_to_Mpc(vel_disp_z_m[ii])
vel_disp_sq_Mpc = 3*(vel_disp_sq_Mpc**2)
#
if (diag_flag == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        print('\nMean recessional velocity (km/s): \n%s'%(vel_disp_mean/1e3),'\nRecessional velocity dispersion (z-direction) by cluster (km/s): \n%s'%(vel_disp_z_m/1e3),'\nTOTAL velocity dispersion by cluster (km/s): \n%s'%(np.sqrt(vel_disp_sq_m)/1e3),'\nTOTAL velocity dispersion (squared) by cluster (km/s): \n%s'%(vel_disp_sq_m/1e3))
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
cluster_stellar_mass_total_log10_Msol = np.array([0]*6,dtype='float128')
#
for ii in range(len(cluster_mass)): 
    cluster_stellar_mass_total_log10_Msol[ii] = np.sum(cluster_mass[ii])
#
####### convert masses from [log10(solar mass units)] to [kg]
#####cluster_stellar_mass_total_kg = np.empty_like(cluster_stellar_mass_total_log10_Msol,dtype='float128')
#####cluster_stellar_mass_total_kg = log10_Msol_to_kg(cluster_stellar_mass_total_log10_Msol)
#
#
#####kinetic_energy = 0
#####for ii in range(len(cluster_stellar_mass_total)):
#####    kinetic_energy+=(0.5*cluster_stellar_mass_total_kg[ii]*vel_disp_sq[ii])
#
#
if (diag_flag == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        print('\nTotal STELLAR Mass by cluster (units: [log_10(M_sol)]:\n%s'%cluster_stellar_mass_total_log10_Msol)
#
## SECTION (3.2): RHS: compute the POTENTIAL ENERGY --> really, compute r_200 & M_200: the radius within which the average density is 200 times the critical density of a spatially flat univervse, and the mass enclosed within a sphere of that radius
#
## r_200: for a derivation, see Notes; first calculate the Hubble Parameter ("H2") for each cluster 
#
H2 = np.array([0]*6,dtype='float128')
for ii in range(len(z_cluster)):
    H2[ii] = friedmann_equation(z_cluster[ii])
#
## r_200 in UNITS [Mpc]
r_200_Mpc = np.array([0]*6,dtype='float128')
for ii in range(len(r_200_m)):
    #r_200_m[ii] = np.sqrt((1/100)*(1/H2[ii])*vel_disp_sq_m[ii])
    #r_200_Mpc = m_to_Mpc(r_200_m)
    r_200_Mpc[ii] = np.sqrt((1/100)*(1/H2[ii])*vel_disp_sq_Mpc[ii])    # dimensionally correct: units of 
#
## M_200: for a derivation, see Notes; in UNITS [Msol]
#
M_200_log10_Msol = np.array([0]*6,dtype='float128')
for ii in range(len(M_200_kg)):
    M_200_log10_Msol[ii] = np.log10((100)*(H2[ii])*(r_200_Mpc[ii]**3)/G)        # dimensionally correct
#
## convert M_200 from kg to log10(M_sol)
#M_200_Msol = np.array([0]*6,dtype='float128')
#M_200_Msol = kg_to_log10_Msol(M_200_kg)
#
if (diag_flag == 1 and project_diagnostic_flag == 2) or project_diagnostic_flag == 1:
    if project_diagnostic_flag == 0:
        pass
    else:
        print('\nHubble factor by cluser:\n%s'%H2,'\n\nr_200 by cluster (Mpc):\n%s'%r_200_Mpc,'\n\nM_200 by cluster (log10(M_sol)):\n%s'%M_200_log10_Msol)
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