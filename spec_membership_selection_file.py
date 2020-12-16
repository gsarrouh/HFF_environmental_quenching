# Created Thu Jul 02 16:23:55 2020
#
#
###################     PROGRAM START     ###################
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
## 7: FIELD lost BELOW limiting mass
#
#
##  import modules
import numpy as np
#
#
#def spec_membership_selection(master_cat,z_cutoff):
#
#####
## INDENT HERE to make this a function...
#
mem_spec = np.array([[0]*6]*2)       # initialize arrays to track cluster members, field, false pos/neg by cluster
field_spec = np.array([[0]*6]*2)     # row_1=SF; row_2=Q; for all arrays
far_field_spec = np.array([[0]*6]*2)
pos_spec = np.array([[0]*6]*2)
neg_spec = np.array([[0]*6]*2)
other_member_spec = np.array([[0]*6]*2)                # track objects outside of (phot + spec) subsample
lost_due_to_buffer_spec = np.array([[0]*6]*2)
tctc_spec = np.array([[0]*6]*2)
#
## The following loop isolates the (spec + phot) sample, i.e. 'sub'=1, and makes the cuts defined above, assigning different classifications to the MEMBER FILTER
for counter in range(len(master_cat)):
    if master_cat['sub'][counter] == 1:                   # sub=1 identifies subsample with both spec & phot
        if master_cat['type'][counter]==1:                # type=1 identifies SF sample
            if abs(master_cat['z_clusterspec'][counter]) > z_cutoff_field[0] and abs(master_cat['z_clusterphot'][counter]) > z_cutoff_field[1]:
                if master_cat['z_peak'][counter] > z_field_bounds[0] and master_cat['z_peak'][counter] < z_field_bounds[1]:
                    if master_cat['ang_dist'][counter] > ang_dist_TOL:
                        if cluster_field_inclusion_flag == 1:
                            master_cat['member'][counter] = 1         # member=1 for FIELD
                        elif cluster_field_inclusion_flag == 0:
                            master_cat['member'][counter] = -99
                        for ii in range(len(field_spec[0])):
                            if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                                field_spec[0][ii]+=1
                    else:
                        if cluster_field_inclusion_flag == 1:
                            master_cat['member'][counter] = 8         # member=1 for FIELD
                        elif cluster_field_inclusion_flag == 0:
                            master_cat['member'][counter] = -99
                        for ii in range(len(field_spec[0])):
                            if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                                tctc_spec[0][ii]+=1
                else:
                    master_cat['member'][counter] = 4         # member=4 for FAR FIELD
                    for ii in range(len(far_field_spec[0])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                            far_field_spec[0][ii]+=1

            elif master_cat['lmass'][counter] < membership_correction_lower_mass:       # for the 1st two mass bins, do not apply spectroscopic cut at all. just phot cut, due to not applying the memberhsip correction to these bins due to insufficient spectroscopic completeness (<5%)
                if abs(master_cat[counter]['z_clusterphot']) < z_cutoff_lo_mass:
                    master_cat[counter]['member'] = 0           # member=0 is secure cluster member
                    for ii in range(len(mem_spec[0])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field galaxies by cluster
                            mem_spec[0][ii]+=1
                else:
                    for ii in range(len(lost_due_to_buffer_spec[0])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of cluster members by cluster
                            lost_due_to_buffer_spec[0][ii]+=1

            elif master_cat['lmass'][counter] >= membership_correction_lower_mass:
                if abs(master_cat['z_clusterspec'][counter]) > z_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) < z_cutoff[1]: #SF_cutoff[1]: #
                    master_cat['member'][counter] = 2         # member=2 for FALSE POSITIVE
                    for ii in range(len(pos_spec[0])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of false pos by cluster
                            pos_spec[0][ii]+=1
                elif abs(master_cat['z_clusterspec'][counter]) < z_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) > z_cutoff[1]: #SF_cutoff[1]: #
                    master_cat['member'][counter] = 3         # member=3 for FALSE NEGATIVE
                    for ii in range(len(neg_spec[0])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of false neg by cluster
                            neg_spec[0][ii]+=1
                elif abs(master_cat['z_clusterspec'][counter]) < z_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) < z_cutoff[1]: #SF_cutoff[1]: #
                    master_cat['member'][counter] = 0         # member=0 for cluster MEMBERS
                    for ii in range(len(mem_spec[0])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of cluster members by cluster
                            mem_spec[0][ii]+=1
                else:
                    for ii in range(len(lost_due_to_buffer_spec[0])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of cluster members by cluster
                            lost_due_to_buffer_spec[0][ii]+=1
        elif master_cat['type'][counter]==2:                # type=2 identifies Q sample
            if abs(master_cat['z_clusterspec'][counter]) > z_cutoff_field[0] and abs(master_cat['z_clusterphot'][counter]) > z_cutoff_field[1]:
                if master_cat['z_peak'][counter] > z_field_bounds[0] and master_cat['z_peak'][counter] < z_field_bounds[1]:
                    if master_cat['ang_dist'][counter] > ang_dist_TOL:
                        if cluster_field_inclusion_flag == 1:
                            master_cat['member'][counter] = 1         # member=1 for FIELD
                        elif cluster_field_inclusion_flag == 0:
                            master_cat['member'][counter] = -99       # ember == -99 to remove this galaxy from further consideration and all further calculations
                        for ii in range(len(field_spec[1])):
                            if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                                field_spec[1][ii]+=1
                    else:
                        if cluster_field_inclusion_flag == 1:
                            master_cat['member'][counter] = 8         # member=1 for FIELD
                        elif cluster_field_inclusion_flag == 0:
                            master_cat['member'][counter] = -99
                        for ii in range(len(field_spec[0])):
                            if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                                tctc_spec[1][ii]+=1
                else:
                    master_cat['member'][counter] = 4         # member=4 for FAR FIELD
                    for ii in range(len(far_field_spec[1])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                            far_field_spec[1][ii]+=1
                            #
                #
            elif master_cat['lmass'][counter] < membership_correction_lower_mass:       # for the 1st two mass bins, do not apply spectroscopic cut at all. just phot cut, due to not applying the memberhsip correction to these bins due to insufficient spectroscopic completeness (<5%)
                if abs(master_cat[counter]['z_clusterphot']) < z_cutoff_lo_mass:
                    master_cat[counter]['member'] = 0           # member=0 is secure cluster member
                    for ii in range(len(mem_spec[1])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field galaxies by cluster
                            mem_spec[1][ii]+=1
                else:
                    for ii in range(len(lost_due_to_buffer_spec[1])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of cluster members by cluster
                            lost_due_to_buffer_spec[1][ii]+=1
            elif master_cat['lmass'][counter] >= membership_correction_lower_mass:
                if abs(master_cat['z_clusterspec'][counter]) > z_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) < z_cutoff[1]: #Q_cutoff[1]: #
                    master_cat['member'][counter] = 2         # member=2 for FALSE POSITIVE
                    for ii in range(len(pos_spec[1])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of false pos by cluster
                            pos_spec[1][ii]+=1
                elif abs(master_cat['z_clusterspec'][counter]) < z_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) > z_cutoff[1]: #Q_cutoff[1]: #
                    master_cat['member'][counter] = 3         # member=3 for FALSE NEGATIVE
                    for ii in range(len(neg_spec[1])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of false neg by cluster
                            neg_spec[1][ii]+=1
                elif abs(master_cat['z_clusterspec'][counter]) < z_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) < z_cutoff[1]: #Q_cutoff[1]: #
                    master_cat['member'][counter] = 0         # member=0 for cluster MEMBERS
                    for ii in range(len(mem_spec[1])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of cluster members by cluster
                            mem_spec[1][ii]+=1
                else:
                    for ii in range(len(lost_due_to_buffer_spec[1])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of cluster members by cluster
                            lost_due_to_buffer_spec[1][ii]+=1
        else:
            for ii in range(len(mem_spec[1])):
                if master_cat['cluster'][counter] == (ii+1):
                    if master_cat['type'][counter] == 1:
                        other_member_spec[0][ii]+=1
                    elif master_cat['type'][counter] == 2:
                        other_member_spec[1][ii]+=1
#
#    return master_cat, mem_spec, field_spec, pos_spec, neg_spec, lost_due_to_buffer_spec
#
#
print('\n\n"spec_membership_selection_file.py"  terminated successfully.\n')
#
#
###################     PROGRAM END     ###################
