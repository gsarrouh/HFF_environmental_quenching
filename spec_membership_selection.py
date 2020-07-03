# Created Thu Jul 02 16:23:55 2020
#
#
###################     PROGRAM START     ###################
#
#
#
## INDENT HERE for diagnostic 
mem = np.array([[0]*6]*2)       # initialize arrays to track cluster members, field, false pos/neg by cluster
field = np.array([[0]*6]*2)     # row_1=SF; row_2=Q; for all arrays
pos = np.array([[0]*6]*2)
neg = np.array([[0]*6]*2)
other_member = 0                 # track objects outside of (phot + spec) subsample
lost_due_to_buffer = np.array([[0]*6]*2)
#
## The following loop isolates the (spec + phot) sample, i.e. 'sub'=1, and makes the cuts defined above, assigning different classifications to the MEMBER FILTER 
for counter in range(len(master_cat)):
    if master_cat['sub'][counter] == 1:                   # sub=1 identifies subsample with both spec & phot
        if master_cat['type'][counter]==1:                # type=1 identifies SF sample
            if abs(master_cat['z_clusterspec'][counter]) > 0.02 and abs(master_cat['z_clusterphot'][counter]) > 0.1: 
                master_cat['member'][counter] = 1         # member=1 for FIELD
                for ii in range(len(field[0])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                        field[0][ii]+=1
            elif abs(master_cat['z_clusterspec'][counter]) > z_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) < z_cutoff[1]: #SF_cutoff[1]: #
                master_cat['member'][counter] = 2         # member=2 for FALSE POSITIVE
                for ii in range(len(pos[0])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of false pos by cluster
                        pos[0][ii]+=1
            elif abs(master_cat['z_clusterspec'][counter]) < z_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) > z_cutoff[1]: #SF_cutoff[1]: #
                master_cat['member'][counter] = 3         # member=3 for FALSE NEGATIVE
                for ii in range(len(neg[0])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of false neg by cluster
                        neg[0][ii]+=1
            elif abs(master_cat['z_clusterspec'][counter]) < z_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) < z_cutoff[1]: #SF_cutoff[1]: #
                master_cat['member'][counter] = 0         # member=0 for cluster MEMBERS
                for ii in range(len(mem[0])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of cluster members by cluster
                        mem[0][ii]+=1
            else: 
                for ii in range(len(lost_due_to_buffer[0])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of cluster members by cluster
                        lost_due_to_buffer[0][ii]+=1
        elif master_cat['type'][counter]==2:                # type=2 identifies Q sample
            if abs(master_cat['z_clusterspec'][counter]) > 0.02 and abs(master_cat['z_clusterphot'][counter]) > 0.1: 
                master_cat['member'][counter] = 1         # member=1 for FIELD
                for ii in range(len(field[1])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of field objects by cluster
                        field[1][ii]+=1
            elif abs(master_cat['z_clusterspec'][counter]) > z_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) < z_cutoff[1]: #Q_cutoff[1]: #
                master_cat['member'][counter] = 2         # member=2 for FALSE POSITIVE
                for ii in range(len(pos[1])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of false pos by cluster
                        pos[1][ii]+=1
            elif abs(master_cat['z_clusterspec'][counter]) < z_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) > z_cutoff[1]: #Q_cutoff[1]: #
                master_cat['member'][counter] = 3         # member=3 for FALSE NEGATIVE
                for ii in range(len(neg[1])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of false neg by cluster
                        neg[1][ii]+=1
            elif abs(master_cat['z_clusterspec'][counter]) < z_cutoff[0] and abs(master_cat['z_clusterphot'][counter]) < z_cutoff[1]: #Q_cutoff[1]: #
                master_cat['member'][counter] = 0         # member=0 for cluster MEMBERS
                for ii in range(len(mem[1])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of cluster members by cluster
                        mem[1][ii]+=1
            else: 
                for ii in range(len(lost_due_to_buffer[1])):
                    if master_cat['cluster'][counter] == (ii+1):  # keep track of cluster members by cluster
                        lost_due_to_buffer[1][ii]+=1
    else: other_member+=1
#
#
#
#
#
###################     PROGRAM END     ###################