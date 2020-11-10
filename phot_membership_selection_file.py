# Created on Thu Jul 02 16:20:11 2020
#
####
#
### This used to Section (5) of "master_data*.py", but I needed to use it in the Variational Analysis loop as well, so instead of copying & pasting the code I put it in this script, which "master_data*.py" calls twice.
#
## FILTER 4 - 'member': :   identifies type of data each object has
## NOTE: this designation is for objects with phot only (sub =2) and spec&phot (sub=1) only, as membership determinaion requires a good 'z_phot' estimate. as such, member=2 & =3 are for spec&phot (sub=1) subsample only, as only they can be classified as false pos/neg
## 0: secure cluster member
## 1: secure field    <-- this comprises the sample of field galaxies to be
##                        compared with cluster, at similar redshifts
## 2: false positive
## 3: false negative
## 4: FAR field       <-- this comprises galaxies well outside the allowable redshift range of our study
## 5: BCGs            <-- identified in the last section of the program, over-writing MEMBER assignment from section 4
## 6: cluster MEMBERS lost BELOW limiting mass
##
#
###################     PROGRAM START     ###################
#
#
### Quiescent sub-sample: segregate into secure member, field
# initialize arrays, format is row_1=SF, row_2=Q
#
SF_field_list = [ [], [], [], [], [], [] ]                   # THESE LISTS WILL STORE THE SMF FIELD SAMPLE MASSES
Q_field_list = [ [], [], [], [], [], [] ]
#
mem_phot = np.array([[0]*6]*2)    # for tracking phot members/field galaxies, by cluster
field_phot = np.array([[0]*6]*2)
far_field_phot = np.array([[0]*6]*2)
other_phot = 0      # objects not in sub=2 phot only subsample
n_phot_only = 0     # number of objects in sub=2 subsample
stars_outliers = 0  # number of objects in sub=0 or sub=3 subsamples
field_outliers = np.array([[0]*6]*2)  # track field objects with very large (i.e. highly discrepant) redshifts, far outside the redshift range of the clusters
n_SF = 0
n_Q = 0
lost_due_to_buffer_phot = np.array([[0]*6]*2)    # objects lost due to buffer b/w definition of cluster and field
#
for counter in range(len(master_cat)):
    if master_cat[counter]['sub'] ==2:      # sub=2 identifies phot only subsample;   this cut works fine
        n_phot_only+=1
        if master_cat[counter]['type'] == 0 or master_cat[counter]['type'] ==3: #skip stars and outliers
            stars_outliers+=1               # this cut does nothing
            pass
        elif master_cat[counter]['type'] ==1:       # type=1 identifies SF
            n_SF+=1
            if abs(master_cat[counter]['z_clusterphot']) > z_cutoff_field[1]:     # identify field galaxies
                if master_cat['z_peak'][counter] > z_field_bounds[0] and master_cat['z_peak'][counter] < z_field_bounds[1]:
                    if cluster_field_inclusion_flag == 1:
                        master_cat['member'][counter] = 1         # member=1 for FIELD
                    elif cluster_field_inclusion_flag == 0:
                        master_cat['member'][counter] = -99
                    for ii in range(len(field_phot[0])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                            SF_field_list[ii].append(master_cat['lmass'][counter])
                            field_phot[0][ii]+=1
                else:
                    master_cat['member'][counter] = 4         # member=4 for FAR FIELD
                    for ii in range(len(far_field_phot[0])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                            far_field_phot[0][ii]+=1
                            #
                #
            elif abs(master_cat[counter]['z_clusterphot']) < z_cutoff[1]:
                master_cat[counter]['member'] = 0           # member=0 is secure cluster member
                for ii in range(len(mem_phot[0])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of field galaxies by cluster
                        mem_phot[0][ii]+=1
            else:
                for ii in range(len(lost_due_to_buffer_phot[0])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of field galaxies by cluster
                        lost_due_to_buffer_phot[0][ii]+=1
        elif master_cat[counter]['type'] ==2:       #Q
            n_Q+=1
            if abs(master_cat[counter]['z_clusterphot']) > z_cutoff_field[1]:     # identify field galaxies
                if master_cat['z_peak'][counter] > z_field_bounds[0] and master_cat['z_peak'][counter] < z_field_bounds[1]:
                    if cluster_field_inclusion_flag == 1:
                        master_cat['member'][counter] = 1         # member=1 for FIELD
                    elif cluster_field_inclusion_flag == 0:
                        master_cat['member'][counter] = -99
                    for ii in range(len(field_phot[1])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                            Q_field_list[ii].append(master_cat['lmass'][counter])
                            field_phot[1][ii]+=1
                else:
                    master_cat['member'][counter] = 4         # member=4 for FAR FIELD
                    for ii in range(len(far_field_phot[1])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                            far_field_phot[1][ii]+=1
                            #
                #
            elif abs(master_cat[counter]['z_clusterphot']) < z_cutoff[1]:
                master_cat[counter]['member'] = 0           # member=0 is secure cluster member
                for ii in range(len(mem_phot[1])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of field galaxies by cluster
                        mem_phot[1][ii]+=1
            else:
                for ii in range(len(lost_due_to_buffer_phot[1])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of field galaxies by cluster
                        lost_due_to_buffer_phot[1][ii]+=1
    else:
        other_phot+=1
#
#
#
print('\n\n"phot_membership_selection_file.py"  terminated successfully.\n')
#
#
###################     PROGRAM END     ###################
