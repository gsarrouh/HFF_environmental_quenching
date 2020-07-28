# Created Thu Jul 02 16:23:55 2020
#
#
###################     PROGRAM START     ###################
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
#
## The following loop isolates the (spec + phot) sample, i.e. 'sub'=1, and makes the cuts defined above, assigning different classifications to the MEMBER FILTER 
for counter in range(len(master_cat)):
    if master_cat['sub'][counter] == 1:                   # sub=1 identifies subsample with both spec & phot
        if master_cat['type'][counter]==1:                # type=1 identifies SF sample
            if abs(master_cat['z_clusterspec'][counter]) > z_cutoff_field[0] and abs(master_cat['z_clusterphot'][counter]) > z_cutoff_field[1]: 
                if master_cat['z_peak'][counter] > z_field_bounds[0] and master_cat['z_peak'][counter] < z_field_bounds[1]:
                    master_cat['member'][counter] = 1         # member=1 for FIELD
                    for ii in range(len(field_spec[0])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                            field_spec[0][ii]+=1
                else:
                    master_cat['member'][counter] = 4         # member=4 for FAR FIELD
                    for ii in range(len(far_field_spec[0])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                            far_field_spec[0][ii]+=1
                
            elif abs(master_cat['z_clusterspec'][counter]) > z_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) < z_cutoff[1]: #SF_cutoff[1]: #
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
                    master_cat['member'][counter] = 1         # member=1 for FIELD
                    for ii in range(len(field_spec[1])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                            field_spec[1][ii]+=1
                else:
                    master_cat['member'][counter] = 4         # member=4 for FIELD
                    for ii in range(len(far_field_spec[1])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                            far_field_spec[1][ii]+=1
            elif abs(master_cat['z_clusterspec'][counter]) > z_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) < z_cutoff[1]: #Q_cutoff[1]: #
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