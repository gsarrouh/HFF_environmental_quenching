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
### (2):   by and large a copy of master_data*.py, imports and anlyzes 
###        the field sample from the Hubble Parallel Fields ; 
###        MAIN program: masterparallel_2.py;
#
### (3):   determine LIMITING MASS - self-explanatory; MAIN program: ******.py 
#
### (4):   REDSHIFT plots; executes a file which produces figures 
###        ANALYZING the QUALITY of the data set; 
###        MAIN program: master_zplots_2_final.py;
#
### (5):   determine SCALE RADIUS and MASS "r_200" & "M_200" of each cluster
###        MAIN program: velocity_dispersion.py
#
### (6):   produce Stellar Mass Function (SMF); ;
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
z_cutoff = [0.012,0.055]#1st try[0.02,0.06]     # [spec,phot] cutoffs for cluster membership
z_cutoff_field = [0.08,0.15]    # same but for definition of "field" galaxies
bin_width = 0.4  # of the SMF, in dex
num_bins_SF_pos_neg = [6.64,8.527,10.413 ,12.3]#[6.6,9.05,10.13,12.3]   # bin edges of the false pos/neg SF histogram, in units of log10(mass); METHOD = 1
num_bins_Q_pos_neg = [6.64,8.527, 10.413,12.3]#1st try -[7.05,10.11,12.3]  # bin edges of the false pos/neg Q histogram;    METHOD = 2
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
diagnostic_round_flag = 2                  # variational analysis performed in 2 rounds: 1st round (flag==1), try all possible cuts; 2nd round (flag==2), compare the top 6 (3 best SF/Q)
#
#
## PROJECT PLOT FLAG
project_plot_flag = 2        # 0=off, make no figures;  1=on, make all figures;  2=on, allow individual figures
#
#
##       0=off, skip section;     1=on, execute section
#
#
section_1_flag = 1                 # cluster catalogue data prep
section_2_flag = 1                 # parallel field catalogue data prep# z plots
section_3_flag = 1                 # limiting mass
section_4_flag = 0                 # z-plots
section_5_flag = 0                 # velocity dispersion, r_200, M_200 calculation
section_6_flag = 0                 # SMF
#    
#
## MAY NEED TO EDIT: choose the filter in which to determine limiting mass
limiting_mass_flag = 1             #   1 = F160W;   2 = F814W
## Update the user on what this program will run
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
    print('Section 4: Prepare redshift & UVJ classification FIGURES ("master_zplots*.py")')

#
if section_5_flag == 1:
    print('Section 5: calculate Velocity Dispersion, r_200, M_200 ("velocity_dispersion.py")')
#
#
if section_6_flag == 1:
    print('Section 6: Produce SMF ("master_smf*.py")')
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
## SECTION (4): SF/Q classification & spec. subsample selection FIGURES
#
## Call and execute the "master_zplots*.py" file, to create plots assessing the quality of the data and visualizes galaxy classification (i.e. SF/Q).  Fig. 1: z_phot v z_spec;  Fig. 2: cluster members/field/false pos/false neg;  Fig. 3: UVJ diagram
#
#
#
if section_4_flag == 1:
    print('\nBeginning "master_zplots*.py"')
    exec(open('master_zplots_2_final.py').read())      #opens and executes the script 
#
#
#
#
## SECTION (5): SMF
#
#
#
if section_5_flag == 1:
    print('\nBeginning "velocity_dispersion.py"')
    exec(open('velocity_dispersion.py').read())      #opens and executes the script 
#
#
#
#
## SECTION (5): SMF
#
#
#
if section_6_flag == 1:
    print('\nBeginning "master_smfz*.py"')
    exec(open('master_smfz_9_final.py').read())      #opens and executes the script 
#
#
#
## SECTION (6): 
#
#
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
print('\n\nProgram terminated successfully.')
#
#
#
#
#
###################     PROGRAM END     ###################