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
###        determine cluster membership;  
###        MAIN program: master_data_7_final.py;  
###        SUB-programs: spec_completeness_binning.py, correction_factors.py; 
#
### (2):   REDSHIFT plots; executes a file which produces figures 
###        ANALYZING the QUALITY of the data set; 
###        MAIN program: master_zplots_2_final.py;
#
### (3):   determine LIMITING MASS - self-explanatory; MAIN program: ******.py 
#
### (4):   produce Stellar Mass Function (SMF); 
###        MAIN program: master_smfz_8_final.py;
#
### (5):   ;
#
### (6):   ;
#
### (7):   ;
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
## PROJECT TIME_FLAG: START
## MAY NEED TO EDIT
## Superior time_flag which measures execution time of the entire project
project_time_flag = 1     # 0= off;   1=on, time program, not section-by-section by in its entirety. to time individual sections within programs, go to that file and turn on individual time flags by searching "time_flag";
#
#
## SECTION (0): modules & definitions
#
## MODULES
#
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
#
## DEFINITIONS
#
## MAY NEED TO EDIT: hard code the cluster membership definition cuts if not running Variational Analysis
z_cutoff = [0.01,0.03]     # [spec,phot] cutoffs for cluster membership
z_cutoff_field = [0.08,0.15]    # same but for definition of "field" galaxies
bin_width = 0.2  # of the SMF, in dex
#
#
## SECIONT (0.1): FLAGS
## MAY NEED TO EDIT: the following flags allow you to choose which sections (i.e. which progams) to execute (e.g. for testing purposes, you may want to run "master_data*.py" once to load the data, but not again on subsequent executions of the program for debugging)
#
## ADAM's flag
adams_flag = 1                     # diagnostic summary table to account for every object given a redshift cut
#                                  #  0=off;  1=on
#
## PROJECT DIAGNOSTIC FLAG
project_diagnostic_flag = 2        # 0=off, turn OFF ALL diagnostic flags in all files;    1=on, turn ON ALL diagnostics for every file that will be executed; or go into each individual program and turn on specific diagnostic flags by searching "diag_flag";      2=on, allow INDIVIDUAL diagnostic flags within separate programs to be turned on/off
#
#
## PROJECT DIAGNOSTIC FLAG
project_master_variational_flag = 0        # 0=off, don't perform variational analaysis;  1=on, do it. the analysis may also be turned on/off within "master_data*.py"
#
#
##       0=off, skip section;     1=on, execute section
#
#
section_1_flag = 1                 # data prep
section_2_flag = 0                 # z plots
section_3_flag = 0                 # limiting mass
section_4_flag = 0                 # SMF
#    
## Update the user on what this program will run
#
print('"main_project_file" will run the following:')
#
if section_1_flag == 1:
    print('Section 1: "master_data*.py"')
#
if section_2_flag == 1:
    print('Section 2: "master_zplots*.py"')
#
if section_3_flag == 1:
    print('Section 3: "data_mass_completeness*.py"')
#
if section_4_flag == 1:
    print('Section 4: "master_smf*.py"')
#
## SECTION (1): main DATA preparation file; imports raw data, sorts it by data type (i.e. photometric/spectroscopic data, stars, etc...), and classifies all galaxies with good photometric redshift estimates as Star-Forming (SF) or Quiescent (Q); it then runs a VARIATIONAL ANALYSIS to determine optimal redshift definitions (redshift "cuts") for cluster membership, based on which cuts yield an equal number of false positives/negatives in each mass bin; executes redshift cuts, classifies SF/Q galaxies as either cluster members, false pos/neg, or field galaxy; and finally, checks the catalogue for Brightest Cluster Galaxies (bCGs); MAIN program: master_data_7_final.py; SUB-programs: spec_completeness_binning.py, correction_factors.py;

# by data type, and classifies them as SF/Q; runs a VARIATIONAL ANALYSIS to determine optimal definition of cluster membership; executes redshift cuts, classifies SF/Q galaxies as cluster members; finally, checks the catalogue for Brightest Cluster Galaxies (bCGs);
#
## Call and execute the "master_data*.py" file, to import and classify data 
#
#
#
if section_1_flag == 1:
    print('\nBeginning master_data*.py')
    exec(open('master_data_7_final.py').read())      #opens and executes the script 
#
if project_master_variational_flag == 1:
    sys.exit()
#
#
## SECTION (2): 
#
## Call and execute the "master_zplots*.py" file, to create plots assessing the quality of the data and visualizes galaxy classification (i.e. SF/Q).  Fig. 1: z_phot v z_spec;  Fig. 2: cluster members/field/false pos/false neg;  Fig. 3: UVJ diagram
#
if section_2_flag == 1:
    print('\nBeginning master_zplots*.py')
    exec(open('master_zplots_2_final.py').read())      #opens and executes the script 
#
#
#
#
## SECTION (3): 
#
#
#
if section_3_flag == 1:
    print('\nBeginning data_mass_completeness*.py')
    exec(open('data_mass_completeness_5.py').read())      #opens and executes the script 
#
#
#
## SECTION (4): 
#
#
#
if section_4_flag == 1:
    print('\nBeginning master_smfz*.py')
    exec(open('master_smfz_8_final.py').read())      #opens and executes the script 
#
#
#
## SECTION (5): Prepare summary tables for ADAMS_FLAG to account for every object
#
if adams_flag == 1:
    pass
#
#
#
#
## SECTION (6): 
#
#
#
#
#
## SECTION (7): 
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