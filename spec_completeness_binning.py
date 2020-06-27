## SECTION (3.2): calculate SPECTROSCOPIC COMPLETENESS correction. basically, look at all the false positives/false negatives, and sort them by type (i.e. SF/Q). then bin them (i.e. make histograms of false pos/neg for each of SF/Q). take their ratio of false pos to false neg, and plot that ratio. it is the correction factor to be applied to the photometric subsample
#
#
SF_pos = []
SF_neg = []
Q_pos = []
Q_neg = []
pos_by_cluster = np.array([[0]*6]*2)    #for tracking false pos/neg by cluster; row_1=SF, row_2=Q
neg_by_cluster = np.array([[0]*6]*2)
objects_below_lim_mass = np.array([0]*6)    # for tracking objects below the limiting mass of each cluster
#
for counter in range(len(master_cat)):
    for ii in range(len(limiting_mass)):
        if master_cat['cluster'][counter] == (ii+1):           # only look at objects above the limiting mass for each cluster
            if master_cat['lmass'][counter] > limiting_mass[ii]:      
                if master_cat['type'][counter] == 1:      # type=1 for SF
                    if master_cat['member'][counter] == 2:     # member=2 for false pos
                        SF_pos.append(master_cat['lmass'][counter])
                        for ii in range(len(pos_by_cluster[0])):
                            if master_cat['cluster'][ii] == (ii+1):
                                pos_by_cluster[0]+=1           # track false pos for SF
                    elif master_cat['member'][counter] == 3:   # member=3 for false neg
                        SF_neg.append(master_cat['lmass'][counter])
                        for ii in range(len(pos_by_cluster[0])):
                            if master_cat['cluster'][ii] == (ii+1):
                                neg_by_cluster[0]+=1           # track false neg for SF
                elif master_cat['type'][counter] == 2:     # type=2 for Q
                    if master_cat['member'][counter] == 2:     # member=2 for false pos
                        Q_pos.append(master_cat['lmass'][counter])
                        for ii in range(len(pos_by_cluster[1])):
                            if master_cat['cluster'][ii] == (ii+1):
                                pos_by_cluster[1]+=1           # track false pos for Q
                    elif master_cat['member'][counter] == 3:   # member=3 for false neg
                        Q_neg.append(master_cat['lmass'][counter])
                        for ii in range(len(neg_by_cluster[1])):
                            if master_cat['cluster'][ii] == (ii+1):
                                neg_by_cluster[1]+=1           # track false neg for Q
            else: 
                objects_below_lim_mass[ii]+=1
#
#
#
## to be used in building strings throughout program
space = ' '   
#
## sort false pos/neg lists in ascending order
SF_pos = np.sort(SF_pos)
SF_neg = np.sort(SF_neg)
Q_pos = np.sort(Q_pos)
Q_neg = np.sort(Q_neg)
#
#
#
# SPEC. BINNING: iterate through different number of histogram bins to see which yields a set of corrections closest in general to ~1; this has now been expanded to test: diag_flag4 = 1, symmetric binning; diag_flag4 = 2, asymmetric binning, equal number of objects in each bin;
#
#
### METHOD 1: bin SF & Q in SYMMETRIC BINS, then compute false pos/neg fractions by mass bin for correction factors. 
#
## number of histogram mass bins to try for spec. completeness correction
num_bins_to_try = [2,3,4,5,6,7,8]   
# write a loop that interatively uses a different number of bins in the histrogram, to isolate the largest number for which all bins have at least one entry; NOTE: the lines that stop the loop have been commented out, to investigate the relative fraction of false pos/neg for each different # of bins
## open a file to print to
method = 1
# open a file for method 1
f = open('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/diagnostic_outputs/spec_binning/symmetric_bins_%s'%z_cutoff[0]+'_spec_cutoff_%s'%z_cutoff[1]+'_phot_cutoff.txt','w+')

for number in range(len(num_bins_to_try)):
    #
    # make histograms
    SF_pos_hist, bins_SF = np.histogram(SF_pos, bins=num_bins_to_try[number], range=range2)
    SF_neg_hist, bins_SF = np.histogram(SF_neg, bins=num_bins_to_try[number], range=range2)
    Q_pos_hist, bins_Q = np.histogram(Q_pos, bins=num_bins_to_try[number], range=range2)
    Q_neg_hist, bins_Q = np.histogram(Q_neg, bins=num_bins_to_try[number], range=range2)
    #
    SF_ratio = SF_pos_hist / SF_neg_hist
    Q_ratio = Q_pos_hist / Q_neg_hist
    for jj in range(len(SF_pos_hist)):
        if SF_pos_hist[jj] == 0 and SF_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
            SF_ratio[jj] = 1
    for jj in range(len(Q_pos_hist)):
        if Q_pos_hist[jj] == 0 and Q_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
            Q_ratio[jj] = 1
    #
    ## compute variance of SF/Q ratios from 1
    SF_var = np.sum((1 - SF_ratio)**2)
    Q_var = np.sum((1 - Q_ratio)**2)
    total_var = SF_var + Q_var
    ## prepare what to write to file
    bin_entry1 = str(z_cutoff[0])+space+str(z_cutoff[1])+space+'SF'+space+str(method)+space+str(num_bins_to_try[number])+space+'%.3f'%SF_var+space+'%.3f'%total_var
    bin_entry2 = str(z_cutoff[0])+space+str(z_cutoff[1])+space+' Q'+space+str(method)+space+str(num_bins_to_try[number])+space+'%.3f'%Q_var+space+'%.3f'%total_var
    writer = '%s'%str(bin_entry1)+'\n'+'%s'%str(bin_entry2)+'\n#\n'
    f.write(writer)
##        
f.close()
#
#
#
#
#
### METHOD 2: ASYMMETRIC BINS for bins with ~EQUAL NUMBER OF OBJECTS in each bin (for each of the SF/Q false pos/neg pairs), return the asymmetric bin edges and compute the midpoints b/w the *_pos/*_neg bin edges. Use these midpoints as the bin edges for a new histogram for the *_pos/*_neg lists, and re-compute the false pos/neg ratio for each of SF/Q. Then print relevant data to a file.
#
method = 2



## find index corresponding to the 25th, 50th, 75th & 100th percentiles
num_bins_to_try = [4,5,6,7,8,9]      # set number of asymmetric bins edges, bins will be one less; run through in a loop and PRINT TO OUTPUT FILE
## open a file to print to
f = open('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/diagnostic_outputs/spec_binning/asymmetric_equal_num_bins_%s'%z_cutoff[0]+'_spec_cutoff_%s'%z_cutoff[1]+'_phot_cutoff.txt','w+')
for number in range(len(num_bins_to_try)):
    ## compute the index corresponding to the first evenly-space bin edge
    SF_pos_index = int(np.ceil(len(SF_pos)/num_bins_to_try[number]))
    SF_neg_index = int(np.ceil(len(SF_neg)/num_bins_to_try[number]))
    Q_pos_index = int(np.ceil(len(Q_pos)/num_bins_to_try[number]))
    Q_neg_index = int(np.ceil(len(Q_neg)/num_bins_to_try[number]))
    ## build arrays over range [7.3,12.3] with 4 bins, where each bin has ~equal # of objects in it; then compute the midpoints between the false pos bin edges and the false neg bin edges
    num_bins_SF = [[],[]]     # [pos, neg]
    num_bins_Q = [[],[]]     # [pos, neg]
    #a=0
    #b=0
    for ii in range(num_bins_to_try[number]-1):   # -1 b/c there are 1 fewer bins than bin edges, which is what is stored in 'num_bins_to_try'
        #if ii != (num_bins_to_try[number]-2):      # don't fill the last cell in array, in case the index is too large and raises an error. I'm replacing it below in the 'else' statement anyways
        num_bins_SF[0].append(SF_pos[(ii*SF_pos_index)])
        num_bins_SF[1].append(SF_neg[(ii*SF_neg_index)])
        num_bins_Q[0].append(Q_pos[(ii*Q_pos_index)])
        num_bins_Q[1].append(Q_neg[(ii*Q_neg_index)])
    #        a+=1
    num_bins_SF[0][-1] = (SF_pos[-1])
    num_bins_SF[1][-1] = (SF_neg[-1])
    num_bins_Q[0][-1] = (Q_pos[-1])
    num_bins_Q[1][-1] = (Q_neg[-1])
    #        b+=1
    #print(a)
    #print(b)
    num_bins_SF = np.array(num_bins_SF)
    num_bins_Q = np.array(num_bins_Q)
    print('# bins: %s'%(num_bins_to_try[number]-1))
    print('length of num_bins_* array: %s'%len(num_bins_SF[0]))
    #print(('num_bins_SF: %s'%num_bins_SF)
    #
    ## compute midpoints between false pos/neg bin edges
    bin_edge_means_SF = np.mean(num_bins_SF,axis=0)
    bin_edge_means_Q = np.mean(num_bins_Q,axis=0)
    ## set first/last entry as limits of mass range for smf
    bin_edge_means_SF[0] = range2[0]         # reminder: row [0] = false pos; row [1] = false neg
    bin_edge_means_SF[-1] = range2[-1]
    bin_edge_means_Q[0] = range2[0]
    bin_edge_means_Q[-1] = range2[-1]
    print('length of bin_edge_means_* array: %s'%len(bin_edge_means_SF))
    #print('bin_edge_means_SF: %s'%bin_edge_means_SF)
    ## build new histograms
    SF_pos_hist, bins_SF = np.histogram(SF_pos, bins=bin_edge_means_SF, range=range2)
    SF_neg_hist, bins_SF = np.histogram(SF_neg, bins=bin_edge_means_SF, range=range2)
    Q_pos_hist, bins_Q = np.histogram(Q_pos, bins=bin_edge_means_Q, range=range2)
    Q_neg_hist, bins_Q = np.histogram(Q_neg, bins=bin_edge_means_Q, range=range2)
    ## find ratios, set ratio==1 for bins with no false pos or false neg
    SF_ratio = SF_pos_hist / SF_neg_hist
    Q_ratio = Q_pos_hist / Q_neg_hist
    for jj in range(len(SF_pos_hist)):
        if SF_pos_hist[jj] == 0 and SF_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
            SF_ratio[jj] = 1
    for jj in range(len(Q_pos_hist)):
        if Q_pos_hist[jj] == 0 and Q_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
            Q_ratio[jj] = 1
    #
    ## prepare what to write to file
    bin_entry1 = '\n\n****************\nFOR %s'%(num_bins_to_try[number]-2)+' BINS\n\nSF_pos bin edges: '+str(num_bins_SF[0])+'\nSF_neg bin edges: '+str(num_bins_SF[1])+'\n\nQ_pos bin edges: '+str(num_bins_Q[0])+'\nQ_neg bin edges: '+str(num_bins_Q[1])
    bin_entry2 = '\nHistogram based on midpoints of the above\nSF:\nBINS: '+str(bins_SF)+'\nSF pos: '+str(SF_pos_hist)+'\nSF neg: '+str(SF_neg_hist)+'\nSF false pos/neg ratio: '+str(SF_ratio)+'\n\nQ:\nBINS: '+str(bins_Q)+'\nQ pos: '+str(Q_pos_hist)+'\nQ neg: '+str(Q_neg_hist)+'\nQ false pos/neg ratio: '+str(Q_ratio)
    # write it
    if number == 0:      # setup header of document
        asterisks = '*********************************************************\n*********************************************************'
        header1 = "\n\n### This file shows first a sorted list of false pos/neg \nfor the SF/Q samples. It then lists bin edges for \nASYMMETRIC BINNING, reporting an equal number of objects\nin each bin. It then computes the midpoint between the\nbin edges of the i'th false pos bin and the i'th false neg\nbin. Histograms are then computed for false pos/neg using\nthose midpoints as the histogram bin edges for both, and the\ncorresponding false pos/neg RATIO is printed.\nNOTE: bins which have a count of 0 for both false pos \n& false neg have a ratio set equal to 1.\nYou'll just have to make up your own mind from there###\n"
        header2 = '\nSF_pos sorted list: \n'+str(SF_pos)+'\n\nSF_neg sorted list: \n'+str(SF_neg)+'\n\nQ_pos sorted list: \n'+str(Q_pos)+'\n\nQ_neg sorted list: \n'+str(Q_neg)+'\n'

        writer = asterisks+'%s\n'%header1+asterisks+'\n\n\nSORTED LISTS%s'%header2+str(bin_entry1)+'\n%s'%str(bin_entry2)+'\n'
        f.write(writer)
    else: 
        writer = '\n%s'%str(bin_entry1)+'\n'+'\n%s'%str(bin_entry2)+'\n'
        f.write(writer)
##        
f.close()
#





#
#
### METHOD 1: ASYMMETRIC BINS for bins with ~EQUAL AMOUNTS OF MASS in each bin (for each of the SF/Q false pos/neg pairs), return the asymmetric bin edges and compute the midpoints b/w the *_pos/*_neg bin edges. Use these midpoints as the bin edges for a new histogram for the *_pos/*_neg lists, and re-compute the false pos/neg ratio for each of SF/Q. Then print relevant data to a file.
#
method = 3


## find index corresponding to the 25th, 50th, 75th & 100th percentiles
num_bins_to_try = [3,4,5,6,7,8,9]      # set number of asymmetric bins edges, bins will be one less; run through in a loop and PRINT TO OUTPUT FILE
## open a file to print to
for number in range(len(num_bins_to_try)):
    #
    ## SYMMETRIC BINNING - METHOD 1
    f = open('/Users/gsarrouh/Documents/Programs/Python/nserc17/working_data/diagnostic_outputs/spec_binning/asymmetric_equal_mass_bins_%s'%z_cutoff[0]+'_spec_cutoff_%s'%z_cutoff[1]+'_phot_cutoff.txt','w+')

    ## compute index corresponding to (1/num_bins_to_try[number])% of the mass for all clusters, e.g. if num_bins_to_try[number] = 4, compute how much mass constitutes 25% of the total cluster mass
    SF_pos_sum_index = np.sum(SF_pos)/num_bins_to_try[number] 
    SF_neg_sum_index = np.sum(SF_neg)/num_bins_to_try[number]
    Q_pos_sum_index = np.sum(Q_pos)/num_bins_to_try[number]
    Q_neg_sum_index = np.sum(Q_neg)/num_bins_to_try[number]
    ## build arrays over range [7.3,12.3] with num_bins_to_try[number] bins, where each bin has ~equal amounts of mass in it; then compute the midpoints between the false pos bin edges and the false neg bin edges
    num_bins_SF_index = [[],[]]     # [pos, neg]
    num_bins_Q_index = [[],[]]     # [pos, neg]
    # look through SF/Q false pos/neg lists one galaxy at a time, summing their mass until you reach the amount stored in *_sum_index (e.g. SF_pos_sum_index for SF_pos list), then record that index as a bin edge
    #num_bins_SF[0].append(SF_pos[0])
    ## SF_posfor ii in range(num_bins_to_try[number]):
    for ii in range(num_bins_to_try[number]):
        mass_sum = 0
        for jj in range(len(SF_pos)):
            mass_sum = mass_sum + SF_pos[jj]
            if mass_sum >= ii*SF_pos_sum_index:
                print('\nSF_pos\nmass sum : %s'%mass_sum,'\nmass limit index equals: %s'%(ii*SF_pos_sum_index),'jj: %s'%jj)
                num_bins_SF_index[0].append(jj)
                break
    ## SF_neg
    for ii in range(num_bins_to_try[number]):
        mass_sum = 0
        for jj in range(len(SF_neg)):
            mass_sum = mass_sum + SF_neg[jj]
            if mass_sum >= ii*SF_neg_sum_index:
                print('\nSF_neg\nmass sum : %s'%mass_sum,'\nmass limit index equals: %s'%(ii*SF_neg_sum_index),'jj: %s'%jj)
                num_bins_SF_index[1].append(jj)
                break
    ## Q_pos
    for ii in range(num_bins_to_try[number]):
        mass_sum = 0
        for jj in range(len(Q_pos)):
            mass_sum = mass_sum + Q_pos[jj]
            if mass_sum >= ii*Q_pos_sum_index:
                print('\nQ_pos\nmass sum : %s'%mass_sum,'\nmass limit index equals: %s'%(ii*Q_pos_sum_index),'jj: %s'%jj)
                num_bins_Q_index[0].append(jj)
                break
    ## Q_neg
    for ii in range(num_bins_to_try[number]):
        mass_sum = 0
        for jj in range(len(Q_neg)):
            mass_sum = mass_sum + Q_neg[jj]
            if mass_sum >= ii*Q_neg_sum_index:
                print('\nQ_neg\nmass sum : %s'%mass_sum,'\nmass limit index equals: %s'%(ii*Q_neg_sum_index),'jj: %s'%jj)
                num_bins_Q_index[1].append(jj)
                break
    #print(num_bins_SF)  # diagnostic
    #a=0       # diagnostic
    #b=0       # diagnostic
    ## convert to arrays
    num_bins_SF_index = np.array(num_bins_SF_index)
    num_bins_Q_index = np.array(num_bins_Q_index)
    print('# bins: %s'%(num_bins_to_try[number]-1))
    print('length of num_bins_* array: %s'%len(num_bins_SF_index[0]),'\nlength of SF array: %s'%len(SF_pos),'\nlength of Q array: %s'%len(Q_pos))
    print('indices corresponding to bin edges\nSF: %s'%num_bins_SF_index,'\nQ: %s'%num_bins_Q_index)
    #
    ##
    num_bins_SF = np.empty_like(num_bins_SF_index,dtype='float32')
    num_bins_Q = np.empty_like(num_bins_Q_index,dtype='float32')
    #
    for ii in range(len(num_bins_SF[0])):
        num_bins_SF[0][ii] = SF_pos[num_bins_SF_index[0][ii]]
        num_bins_SF[1][ii] = SF_neg[num_bins_SF_index[1][ii]]
        num_bins_Q[0][ii] = Q_pos[num_bins_Q_index[0][ii]]
        num_bins_Q[1][ii] = Q_neg[num_bins_Q_index[1][ii]]
    ## set first/last entry as limits of mass range for smf
    num_bins_SF[0][0] = range2[0]         # reminder: row [0] = false pos; row [1] = false neg
    num_bins_SF[1][0] = range2[0]
    num_bins_Q[0][0] = range2[0]
    num_bins_Q[1][0] = range2[0]
    num_bins_SF[0][-1] = range2[-1] 
    num_bins_SF[1][-1] = range2[-1]
    num_bins_Q[0][-1] = range2[-1]
    num_bins_Q[1][-1] = range2[-1]
    #
    print('\nSF bin edges are ([pos],[neg]): %s'%num_bins_SF,'\nQ bin edges are: %s'%num_bins_Q)
    ## compute midpoints between false pos/neg bin edges
    bin_edge_means_SF = np.mean(num_bins_SF,axis=0)
    bin_edge_means_Q = np.mean(num_bins_Q,axis=0)
    #
    print('length of bin_edge_means_* array: %s'%len(bin_edge_means_SF))
    #print('bin_edge_means_SF: %s'%bin_edge_means_SF)
    ## build new histograms
    SF_pos_hist, bins_SF = np.histogram(SF_pos, bins=bin_edge_means_SF, range=range2)
    SF_neg_hist, bins_SF = np.histogram(SF_neg, bins=bin_edge_means_SF, range=range2)
    Q_pos_hist, bins_Q = np.histogram(Q_pos, bins=bin_edge_means_Q, range=range2)
    Q_neg_hist, bins_Q = np.histogram(Q_neg, bins=bin_edge_means_Q, range=range2)
    ## find ratios, set ratio==1 for bins with no false pos or false neg
    SF_ratio = SF_pos_hist / SF_neg_hist
    Q_ratio = Q_pos_hist / Q_neg_hist
    for jj in range(len(SF_pos_hist)):
        if SF_pos_hist[jj] == 0 and SF_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
            SF_ratio[jj] = 1
    for jj in range(len(Q_pos_hist)):
        if Q_pos_hist[jj] == 0 and Q_neg_hist[jj] == 0:      # if both lists = 0 for the same bin, that's fine!
            Q_ratio[jj] = 1
    #
    ## prepare what to write to file
    bin_entry1 = '\n\n****************\nFOR %s'%(num_bins_to_try[number]-1)+' BINS\n\nSF_pos bin edges: '+str(num_bins_SF[0])+'\nSF_neg bin edges: '+str(num_bins_SF[1])+'\n\nQ_pos bin edges: '+str(num_bins_Q[0])+'\nQ_neg bin edges: '+str(num_bins_Q[1])
    bin_entry2 = '\nHistogram based on midpoints of the above\nSF:\nBINS: '+str(bins_SF)+'\nSF pos: '+str(SF_pos_hist)+'\nSF neg: '+str(SF_neg_hist)+'\nSF false pos/neg ratio: '+str(SF_ratio)+'\n\nQ:\nBINS: '+str(bins_Q)+'\nQ pos: '+str(Q_pos_hist)+'\nQ neg: '+str(Q_neg_hist)+'\nQ false pos/neg ratio: '+str(Q_ratio)
    # write it
    if number == 0:      # setup header of document
        asterisks = '*********************************************************\n*********************************************************'
        header1 = "\n\n### This file shows first a sorted list of false pos/neg \nfor the SF/Q samples. It then lists bin edges for \nASYMMETRIC BINNING, reporting an equal amount of mass\nin each bin. It then computes the midpoint between the\nbin edges of the i'th false pos bin and the i'th false neg\nbin. Histograms are then computed for false pos/neg using\nthose midpoints as the histogram bin edges for both, and the\ncorresponding false pos/neg RATIO is printed.\nNOTE: bins which have a count of 0 for both false pos \n& false neg have a ratio set equal to 1.\nYou'll just have to make up your own mind from there###\n"
        header2 = '\nSF_pos sorted list: \n'+str(SF_pos)+'\n\nSF_neg sorted list: \n'+str(SF_neg)+'\n\nQ_pos sorted list: \n'+str(Q_pos)+'\n\nQ_neg sorted list: \n'+str(Q_neg)+'\n'

        writer = asterisks+'%s\n'%header1+asterisks+'\n\n\nSORTED LISTS%s'%header2+str(bin_entry1)+'\n%s'%str(bin_entry2)+'\n'
        f.write(writer)
    else: 
        writer = '\n%s'%str(bin_entry1)+'\n'+'\n%s'%str(bin_entry2)+'\n'
        f.write(writer)
##        
f.close()
#