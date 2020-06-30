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
### NOTE: To find a section, search "SECTION (*)" in all caps, where * is the section number you want. To find fields which require user input, search "MAY NEED TO EDIT" in all caps. Some of these fields may direct you to manually enter input into a sub-program. 
#
#
#
### Section summary:
#
### PROGRAM START
#
### (0)    import modules & definitions;
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
### (3):   determine LIMITING MASS; MAIN program: ******.py 
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
#
###################     PROGRAM START     ###################
#
## PROJECT TIME_FLAG: START
## superior time_flag which measures execution time of the entire project
project_time_flag = 1     # 0= off;   1=on, time entire program;
#
#
## SECTION (0): modules & definitions
#
## MODULES
#
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
#
## DEFINITIONS
#
#
#
## SECTION (1): main DATA preparation file; imports raw data, sorts it by data type (i.e. photometric/spectroscopic data, stars, etc...), and classifies all galaxies with good photometric redshift estimates as Star-Forming (SF) or Quiescent (Q); it then runs a VARIATIONAL ANALYSIS to determine optimal redshift definitions (redshift "cuts") for cluster membership, based on which cuts yield an equal number of false positives/negatives in each mass bin; executes redshift cuts, classifies SF/Q galaxies as either cluster members, false pos/neg, or field galaxy; and finally, checks the catalogue for Brightest Cluster Galaxies (bCGs); MAIN program: master_data_7_final.py; SUB-programs: spec_completeness_binning.py, correction_factors.py;

# by data type, and classifies them as SF/Q; runs a VARIATIONAL ANALYSIS to determine optimal definition of cluster membership; executes redshift cuts, classifies SF/Q galaxies as cluster members; finally, checks the catalogue for Brightest Cluster Galaxies (bCGs);
#
#
#
#
#
## SECTION (2): 
#
#
#
#
#
## SECTION (3): 
#
#
#
#
#
## SECTION (4): 
#
#
#
#
#
## SECTION (5): 
#
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