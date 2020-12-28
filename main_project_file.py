# Created on Mon Jun 29 22:01:47 2020
#
################## main_project_file.py ##################
#
#
### WHAT THIS PROGRAM IS:
### This is the main wrapper program for the project "HFF_environmental_quenching", which is a study of the Hubble Frontier Fields Stellar Mass Function (SMF). It is mostly comments, with little code. It calls other programs one at a time to execute section by section the analysis detailed below.
#
### HOW THIS PROGRAM WORKS:
### Each section calls a program file (i.e. *.py script file) which accomplishes a task. Some programs call other programs. They are all listed in the Section summary which follows. This program provides the main organizational framework for the project.
#
#
#
### Section summary:
#
### PROGRAM START
#
### (0)    import modules & definitions;
#
### (0.1)  FLAGS
#
### (1)    main DATA preparation file; imports raw data, sorts it, and
###        determines cluster membership; it also performs a variational
###        analysis (VAR) to determine the optimal memebrship cut which yields
###        the spec. completeness correction factors closest to 1 in each
###        mass bin
###        MAIN program: master_data_7_final.py;
###        SUB-programs: spec_completeness_binning.py, phot_completeness_binning.py,
###        spec_correction_factors.py (VAR), data_mass_completeness*.py (VAR),
###        delz_hist_diagnostic.py, bcg_hist_diagnostic.py,
###        spec_asymmetric_binning.py (VAR),;
#
### (2):   PARALLEL FIELDS: by and large a copy of master_data*.py, imports and anlyzes
###        the field sample from the Hubble Parallel Fields ;
###        MAIN program: masterparallel_2.py;
#
### (3):   determine LIMITING MASS - self-explanatory; MAIN program: ******.py
#
### (4):   import COSMOS/UltraVISTA catalogue and select FIELD SAMPLE
###        MAIN program: UVC_master_data.py;
#
### (5):   REDSHIFT plots; executes a file which produces figures
###        ANALYZING the QUALITY of the data set;
###        MAIN program: master_zplots_2_final.py;
#
### (6):   UVJ Figures; diagnostic & publication quality
###        MAIN program: UVJ_plots.py;
#
### (7):   determine SCALE RADIUS and MASS "r_200" & "M_200" of each cluster
###        MAIN program: velocity_dispersion.py
#
### (8):   produce Stellar Mass Function (SMF); ;
###        MAIN program: master_smfz_9_final.py;
#
#
### PROGRAM END
#
#
#
### NOTE: To find a section, search "SECTION (*)" in all caps, where * is the section number you want. To find fields which require user input, search "MAY NEED TO EDIT" in all caps. Some of these fields may direct you to manually enter input into a sub-program.
#
#
#
###################     PROGRAM START     ###################
#
#
## SECTION (0): modules & definitions
#
## MODULES
#
import sys
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
#
## PROJECT TIME_FLAG: START
## MAY NEED TO EDIT
## Superior time_flag which measures execution time of the entire project
project_time_flag = 1     # 0= off;   1=on, time program, not section-by-section by in its entirety. to time individual sections within programs, go to that file and turn on individual time flags by searching "time_flag";
#
## PROJECT TIME_FLAG START
#
if project_time_flag == 1:
    start_time = time.time()
#
#
## DEFINITIONS & USER INPUTS
#
## USER INPUTS
## MAY NEED TO EDIT: hard code the cluster membership definition cuts if not running Variational Analysis
#
range2_flag = 2            # 0==off, range of SMF set by lowest limiting mass / highest-mass object; 1==on, limits defined by user below
#
if range2_flag == 1:
    range2 = [7.8,12.2]
elif range2_flag == 2:
    range2 = [8.0,12.4]
#
# store cluster redshifts; obtained from https://archive.stsci.edu/prepds/frontier/
z_cluster = [0.396,0.543,0.545,0.375,0.348,0.308]
# cluster_names = ['M0416','M1149','M0717','A370','A1063','A2744']
#
z_cutoff_field = [0.04,0.15]    # same but for definition of "field" galaxies
bin_width = 0.4  # of the SMF, in dex
#
ang_dist_TOL = 0.5                  # angular distance TOLerance (in units of arcmin) for field sample (from cluster images only)
#
#
bins_exempt_from_membership_correction = 0      # the 1st 'n' bins will be exempt from the membership correction due to spectroscopic incompleteness
#
if range2_flag == 1:
    if bins_exempt_from_membership_correction != 4:             # for bins_exempt_from_membership_correction == 0,1,2,3
        membership_correction_lower_mass = range2[0] + (bin_width * bins_exempt_from_membership_correction)      # mass below which the membership correction (i.e. spec. completeness correction) will NOT be applied
    elif bins_exempt_from_membership_correction == 4:
        membership_correction_lower_mass = 12.2
elif range2_flag == 2:
    if bins_exempt_from_membership_correction == 0:
        membership_correction_lower_mass = range2[0]
    elif bins_exempt_from_membership_correction == 1:
        membership_correction_lower_mass = range2[0] + (bin_width * bins_exempt_from_membership_correction)      # mass below which the membership correction (i.e. spec. completeness correction) will NOT be applied
#
#
## choose your binning method for calculating membership correction (i.e. false pos/neg ratios)
## 1 = symmetric binning;   2 = asymmetric binning by # count;   3 = asymmetric binning by mass;   4 = best of all methods; 5 = final values for range2 = [7.8,12.2]
membership_correction_binning_flag = 5
#
if membership_correction_binning_flag == 1:    # symmetric
    z_cutoff = [0.012,0.060]
    num_bins_SF_pos_neg = [7.36,9.007,10.653,12.3]  #2nd try
    num_bins_Q_pos_neg = [7.36,8.348,9.336,10.324,11.312,12.3]
elif membership_correction_binning_flag == 2:     # asymmetric by # count
    z_cutoff = [0.012,0.055]
    num_bins_SF_pos_neg = [7.36,8.77,9.5,10.18,12.3]   #1st try
    num_bins_Q_pos_neg = [7.36,9.41,10.41,12.3]
elif membership_correction_binning_flag == 3:     # asymmetric by mass
    z_cutoff = [0.012,0.055]
    num_bins_SF_pos_neg = [7.36,8.93,10.16,12.3]   #1st try
    num_bins_Q_pos_neg = [7.36,9.43,10.11,10.47,12.3]
elif membership_correction_binning_flag == 4:     # best of all 3 methods
    z_cutoff = [0.012,0.055]
    num_bins_SF_pos_neg = [7.36,8.93,10.16,12.3]   #1st try
    num_bins_Q_pos_neg = [7.36,9.41,10.41,12.3]
elif membership_correction_binning_flag == 5:     # FINAL - SMF range fixed at
    if range2_flag == 1:                # range2 = [7.8,12.2]
        z_cutoff = [0.013,0.055]
        z_cutoff_lo_mass = 0.1
        if bins_exempt_from_membership_correction == 1:
            num_bins_SF_pos_neg = [8.2,9.64,10.11,12.2]   #binning method = 2
            num_bins_Q_pos_neg = [8.2,9.6,10.42,12.2]   #binning method = 2
        elif bins_exempt_from_membership_correction == 2:
            num_bins_SF_pos_neg = [8.6,9.5,10.18,12.2]   #binning method = 2
            num_bins_Q_pos_neg = [8.6,9.47,10.41,12.2]   #binning method = 2
        elif bins_exempt_from_membership_correction == 3:
            num_bins_SF_pos_neg = [8.6,9.5,10.18,12.2]   # PLACEHOLDER
            num_bins_Q_pos_neg = [8.6,9.47,10.41,12.2]   #
        elif bins_exempt_from_membership_correction == 4:
            num_bins_SF_pos_neg = [8.6,9.5,10.18,12.2]   # PLACEHOLDER
            num_bins_Q_pos_neg = [8.6,9.47,10.41,12.2]   #
            #
    elif range2_flag == 2:          # range2 = [8.0,12.4]
        z_cutoff = [0.012,0.055]
        z_cutoff_lo_mass = z_cutoff[1] #0.1
        if bins_exempt_from_membership_correction == 0:
            num_bins_SF_pos_neg = [8.0,8.93,10.16,12.4]   #binning method = 3
            num_bins_Q_pos_neg = [8.0,9.46,10.41,12.4]
        elif bins_exempt_from_membership_correction == 1:
            num_bins_SF_pos_neg = [8.4,9.4,10.18,12.4]   #binning method = 2
            num_bins_Q_pos_neg = [8.4,9.46,10.41,12.4]
#
## define redshift bounds for the FIELD SAMPLE
z_field_bounds_flag = 0         # 0 = user-defined, below;   1 = automatically defined based on definition of z_cutoff_field
if z_field_bounds_flag == 0:
    lower_bound = 0.25               # lower bound of field sample
    upper_bound = 0.75               # upper bound of field sample, in redshift space
    z_field_bounds = [lower_bound,upper_bound]

#
#
#
## SECTION (0.1): FLAGS
## MAY NEED TO EDIT: the following flags allow you to choose which sections (i.e. which progams) to execute (e.g. for testing purposes, you may want to run "master_data*.py" once to load the data, but not again on subsequent executions of the program for debugging)
#
## ADAM's flag
adams_flag = 1                     # diagnostic summary table to account for every object, after EACH calculation, given a redshift cut
#                                  #  0=off;  1=on
#
## PROJECT DIAGNOSTIC FLAG
project_diagnostic_flag = 2        # 0=off, turn OFF ALL diagnostic flags in all files;    1=on, turn ON ALL diagnostics for every file that will be executed; or go into each individual program and turn on specific diagnostic flags by searching "diag_flag";      2=on, allow INDIVIDUAL diagnostic flags within separate programs to be turned on/off
#
#
## PROJECT VARIATIONAL ANALYSIS FLAG
project_master_variational_flag = 0        # 0=off, don't perform variational analaysis;  1=on, do it. the analysis may also be turned on/off within "master_data*.py"
diagnostic_round_flag = 2                  # variational analysis performed in 2 rounds: 1st round (flag==1) DEPRECATED; 2nd round (flag==2), compare/sort results as you like
#
#
## PROJECT PLOT FLAG
project_plot_flag = 2        # 0=off, make no figures;  1=on, make all figures;  2=on, allow individual figures
#
#
##       0=off, skip section;     1=on, execute section
#
#
section_1_flag = 1                 # HFF cluster catalogue data prep
section_2_flag = 1                 # HFF parallel field catalogue data prep
section_3_flag = 0                 # limiting mass # DEPRECATED: now called in "master_data*.py"; s/b set to == 0
section_4_flag = 1                 # import UltraVISTA catalogue for Field SMF
section_5_flag = 0                 # z-plots
section_6_flag = 0                 # UVJ diagram
section_61_flag = 0                # UVJ diagram for UltraVISTA catalogue
section_7_flag = 0                 # velocity dispersion, r_200, M_200 calculation - DEPRECATED; should always be set to ==0
section_8_flag = 1                 # SMF
#
## Update the user on what this program will run
#
## MAY NEED TO EDIT: choose the filter in which to determine limiting mass
limiting_mass_flag = 1             #   1 = F160W+F814W;   2 = F814W; DEPRECATED: SHOULD ALWAYS BE SET TO ==1
lim_mass_offset_flag = 1           #   0 = off, do not apply offset; 1 = on, apply (F160W - F180W) offset
#
#
#
## MCMC FLAGS
## MAY NEED TO EDIT: choose whether to enter the MCMC simulation
mcmc_flag = 1             #   0 = off - skip sim;   1 = on - perform MCMC sim & Exit program;
mcmc_field_flag = 1         # 0 = off - fit cluster SMFs;   1 = on - fit field SMFs
## Update the user on what this program will run
#
## MAY NEED TO EDIT: choose if you want to include in the field sample galaxies drawn from the HFF cluster images (outside the redshift range of the image's respective cluster)
cluster_field_inclusion_flag = 1            # 0 == off (NO), do NOT include galaxies from cluster images; 1 == on (YES), include them
field_smf_schechter_flag = 2                # fitting the SF FIELD SMF with a schechter function:  1==single;   2==double
#
## MAY NEED TO EDIT: choose which populations to fit
SF_flag = 1
Q_flag = 1
T_flag = 1
#
#
#
## Update the user on which programs will be run
#
print('"main_project_file.py" will run the following:')
#
if section_1_flag == 1:
    print('Section 1: Import & prepare CLUSTER DATA ("master_data*.py")')
#
if section_2_flag == 1:
    print('Section 2: Import & prepare PARALLEL FIELD DATA ("master_parallel*.py")')
#
if section_3_flag == 1:
    print('Section 3: Determine LIMITING MASS ("data_mass_completeness*.py")')
    if limiting_mass_flag == 1:
        print('Limiting mass calculated in F160W.')
    elif limiting_mass_flag == 2:
        print('Limiting mass calculated in F814W.')
#
if section_4_flag == 1:
    print('Section 4: Import & prepare Field SMF data from UltraVISTA catalogue ("UVC_master_data.py")')
#
#
if section_5_flag == 1:
    print('Section 5: Prepare redshift FIGURES ("master_zplots*.py")')
#
#
if section_6_flag == 1:
    print('Section 6: Prepare UVJ diagram(s) - HFF ("UVJ_plots.py")')
#
#
if section_61_flag == 1:
    print('Section 6.1: Prepare UVJ diagram(s) - UVC ("UVJ_plots_UVC.py")')
#
#
if section_7_flag == 1:
    print('Section 7: Calculate Velocity Dispersion, r_200, M_200 ("velocity_dispersion.py")')
#
#
if section_8_flag == 1:
    print('Section 8: Produce SMF ("master_smf*.py")')
#
print('\nNotes:')
if z_field_bounds_flag == 0:
    print('z_field_bounds fixed at: %s'%np.round(z_field_bounds,decimals=3))
elif z_field_bounds_flag == 1:
    print('z_field_bounds set by definition of z_cutoff_field.\nz_cutoff_field = %s'%np.round(z_cutoff_field,decimlas=3))
#
print('Cluster cutoff (del_z) [spec,phot]: %s'%z_cutoff)
print('Lo-mass cluster cutoff (del_z) [phot]: %s'%z_cutoff_lo_mass)
print('Field cutoff (del_z) [spec,phot]: %s'%z_cutoff_field+'      NOTE: no spec cutoff applied to field')
#
if range2_flag == 1:
    print('SMF range fixed at: %s'%range2)
elif range2_flag == 0:
    print('SMF range: auto')
#
print('Membership correction not applied below: %.2f'%membership_correction_lower_mass+' (affects the first %i'%bins_exempt_from_membership_correction+' bin(s).)')
#
if lim_mass_offset_flag == 0:
    print('(F160W - F814W) offset: OFF')
elif lim_mass_offset_flag == 1:
    print('(F160W - F814W) offset: ON')
#
if cluster_field_inclusion_flag == 0:
    print('Field SMF: PAR ONLY')
elif cluster_field_inclusion_flag == 1:
    print('Field SMF: CLU + PAR')
#
if project_master_variational_flag == 1:
    print('Membership variational analysis: ON')
elif project_master_variational_flag == 0:
    print('Membership variational analysis: OFF')
#
if mcmc_flag == 0:
    print('MCMC simulation: OFF')
elif mcmc_flag == 1:
    print('MCMC simulation: ON')
    if mcmc_field_flag == 0:
        print('MCMC simulation will fit the CLUSTER SMFs.')
    elif mcmc_field_flag == 1:
        print('MCMC simulation will fit the FIELD SMFs.')
        if field_smf_schechter_flag == 1:
            print('FIELD SMFs: SF: single-schechter; Q/Total: double-schechter')
        elif field_smf_schechter_flag == 2:
            print('FIELD SMFs: SF/Q/Total: double-schechter')
#
## SECTION (1): main DATA preparation file; imports raw data, sorts it by data type (i.e. photometric/spectroscopic data, stars, etc...), and classifies all galaxies with good photometric redshift estimates as Star-Forming (SF) or Quiescent (Q); it then runs a VARIATIONAL ANALYSIS to determine optimal redshift definitions (redshift "cuts") for cluster membership, based on which cuts yield an equal number of false positives/negatives in each mass bin; executes redshift cuts, classifies SF/Q galaxies as either cluster members, false pos/neg, or field galaxy; and finally, checks the catalogue for Brightest Cluster Galaxies (bCGs); MAIN program: master_data_7_final.py; SUB-programs: spec_completeness_binning.py, correction_factors.py;

# by data type, and classifies them as SF/Q; runs a VARIATIONAL ANALYSIS to determine optimal definition of cluster membership; executes redshift cuts, classifies SF/Q galaxies as cluster members; finally, checks the catalogue for Brightest Cluster Galaxies (bCGs);
#
## Call and execute the "master_data*.py" file, to import and classify data
#
#
#
if section_1_flag == 1:
    print('\nBeginning "master_data*.py"')
    if project_master_variational_flag == 1:
        print('"master_data*.py" will perform the Variational Analysis, then Exit.')
    exec(open('master_data_7_final.py').read())      #opens and executes the script
#
if project_master_variational_flag == 1:
    #
    if project_time_flag == 1:
        print('Program "main_project_file.py" took: %s seconds to run Variational Analysis.\n\n' % (time.time() - start_time))
    #
    sys.exit()
#
#
#
#
## SECTION (2): PARALLEL FIELD
#
#
#
if section_2_flag == 1:
    print('\nBeginning "master_parallel*.py"')
    #
    exec(open('master_parallel_2.py').read())      #opens and executes the script
#
#
#
#
## SECTION (3): LIMITING MASS CALCULATION
#
#
#
if section_3_flag == 1:
    print('\nBeginning "data_mass_completeness*.py"')
    #
    if limiting_mass_flag == 1:
        exec(open('data_mass_completeness_F160W.py').read())      #opens and executes the script
    elif limiting_mass_flag == 2:
        exec(open('data_mass_completeness_F814W.py').read())      #opens and executes the script
#
#
#
#
## SECTION (4): import UltraVISTA catalogue and sort SF/Q, and select galaxies for Field SMF above 10^9
#
## Call and execute the "master_zplots*.py" file, to create plots assessing the quality of the data and visualizes galaxy classification (i.e. SF/Q).  Fig. 1: z_phot v z_spec;  Fig. 2: cluster members/field/false pos/false neg;  Fig. 3: UVJ diagram
#
#
#
if section_4_flag == 1:
    print('\nBeginning "UVC_master_data.py"')
    exec(open('UVC_master_data.py').read())      #opens and executes the script
#
#
#
#
#
## SECTION (5): SF/Q classification & spec. subsample selection FIGURES
#
## Call and execute the "master_zplots*.py" file, to create plots assessing the quality of the data and visualizes galaxy classification (i.e. SF/Q).  Fig. 1: z_phot v z_spec;  Fig. 2: cluster members/field/false pos/false neg;  Fig. 3: UVJ diagram
#
#
#
if section_5_flag == 1:
    print('\nBeginning "master_zplots*.py"')
    exec(open('master_zplots_2_final.py').read())      #opens and executes the script
#
#
#
#
## SECTION (6): UVJ Diagram(s) - HFF
#
#
#
if section_6_flag == 1:
    print('\nBeginning "UVJ_plots.py" (HFF)')
    exec(open('UVJ_plots.py').read())      #opens and executes the script
#
#
## SECTION (6.1): UVJ Diagram(s) - UVC
#
#
#
if section_61_flag == 1:
    print('\nBeginning "UVJ_plots.py" (UVC)')
    exec(open('UVJ_plots_UVC.py').read())      #opens and executes the script
#
#
#
#
## SECTION (7): VELOCITY DISPERSION
#
#
#
if section_7_flag == 1:
    print('\nBeginning "velocity_dispersion.py"')
    exec(open('velocity_dispersion.py').read())      #opens and executes the script
#
#
#
#
## SECTION (8): SMF
#
if section_8_flag == 1:
    print('\nBeginning "master_smfz*.py"')
    exec(open('master_smfz_9_final.py').read())      #opens and executes the script
#
#
#
#
#
#
## PROJECT TIME_FLAG END
#
if project_time_flag == 1:
    print('Program "main_project_file.py" took: %s seconds to run.\n\n' % (time.time() - start_time))
#
#
#
print('\n\nProgram "main_project_file.py" terminated successfully.')
#
#
#
#
#
###################     PROGRAM END     ###################
