# Created on Thu Jul 02 16:20:11 2020
#
####
#
### This used to Section (5) of "master_data*.py", but I needed to use it in the Variational Analysis loop as well, so instead of copying & pasting the code I put it in this script, which "master_data*.py" calls twice. 
#
#
###################     PROGRAM START     ###################
#
#
### Quiescent sub-sample: segregate into secure member, field
# initialize arrays, format is row_1=SF, row_2=Q
mem_phot = np.array([[0]*6]*2)    # for tracking phot members/field galaxies, by cluster
field_phot = np.array([[0]*6]*2)
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
                if master_cat[counter]['z_peak'] >0.7 or master_cat[counter]['z_peak'] <0.2:
                    #
                    #pass            ### EDIT HERE if you want to analyze the field same as well
                    #
                    #master_cat[counter]['member'] = 4       # memfield outlier
                    for ii in range(len(field_outliers[0])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field outlier galaxies by cluster
                            field_outliers[0][ii]+=1
                else:
                    master_cat[counter]['member'] = 1               #phot SF field sample
                    for ii in range(len(field_phot[0])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field galaxies by cluster
                            field_phot[0][ii]+=1
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
                if master_cat[counter]['z_peak'] >0.7 or master_cat[counter]['z_peak'] <0.2:
                   #
                   #pass   ##   EDIT HERE TO ANALYZE FIELD SAMPLE
                   #
                   #master_cat[counter]['member'] = 4       # memfield outlier
                    for ii in range(len(field_outliers[1])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field outlier galaxies by cluster
                            field_outliers[1][ii]+=1
                else:
                    master_cat[counter]['member'] = 1               #phot SF field sample
                    for ii in range(len(field_phot[1])):
                        if master_cat['cluster'][counter] == (ii+1):  # keep track of field galaxies by cluster
                            field_phot[1][ii]+=1
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