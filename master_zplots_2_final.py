#Created on Fri Jun 26 05:20:20 2020
#
#
### WHAT THIS PROGRAM DOES:
### This script reads in all data for the Hubble Frontier Fields images and prepares data for plotting and analysis. Key information is summarized in the tables: 
#
#  
### This program creates plots for the master data catalogue, containing:
### Fig. 1: photometric v specroscopic redshift;   Fig. 2: delta-z (member/field/false pos/false neg segregation);   Fig. 3: UVJ diagram( "colour-colour" plot)
#
#
#
### Section summary:
#
### PROGRAM START
#
### Section (0) FLAG: preliminary section to set which plots get created
### Section (1) add PLOT_FLAG_1: create plot comparing photometric v specroscopic redshift (Fig. 1)
### Section (2) add PLOT_FLAG_2: create "delta-z" plot, scatterplot of member/field/false pos/false neg by type (SF/Q - Fig. 2)
### Section (3) add PLOT_FLAG_3: create UVJ diagram per Shipley et al 2018 to segregate between type (SF/Q - i.e. a "colour-colour" plot)
#
### PROGRAM END
#
#
#
###################     PROGRAM START
#
## TIME_FLAG: START
## superior time_flag which supercedes all others and times the entire program
time_flag = 0     # track & print time to execute current section
#
if time_flag == 1:
    start_time = time.time()
#
# this next line is specific to jupyter notebook, and allows for Figure editting in a GUI, instead of inline
###%matplotlib qt
#  
# import modules
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
import astropy
from astropy.table import Table
from astropy.table import Column
import time
#
## definitions
#
## define a function to compute the Y-arrays for the bounds defining spectroscopic outliers
def delz_bound(x):
    y_upper_bound = 0.15*(1+x) + x
    y_lower_bound = -0.15*(1+x) + x
    return y_upper_bound, y_lower_bound
#
#
#
## SECTION (0): set PLOT_FLAGS    # 0=off (i.e. don't make plot), 1=on (i.e. make plot)
### MAY NEED TO EDIT ### plot_flag_*/time_flag_*/diag_flag_*
#
plot_flag_1 = 0           # Fig.1 - z_spec v z_phot plot
plot_flag_2 = 0           # Fig.2 - del_z plot
plot_flag_3 = 1           # Fig.3 - UVJ diagram
plot_flag_4 = 0           # diagnostic plots in greyscale;  0=off;  1=on, plot above log_10(Msol) > threshold;  2=on, plot everything; NOTE: plot_flag_3 must be turned on for this plot to work
#
time_flag_1 = 0
time_flag_2 = 0
time_flag_3 = 0     # 0=off;   1=on;   track & print time to execute current section
#
diag_flag_1 = 0           # Section 1: z_spec v z_phot - compare spec. subsample to "master_data*.py"
diag_flag_2 = 0           #
diag_flag_3 = 0           # Section 3: UVJ diagram
#
## SECTION (1): create plot comparing photometric ('z_peak') v specroscopic redshift ('z_spec') (Fig. 1)
#
if plot_flag_1 == 1:
    #
    ## PLOT 1: make z_spec vs z_phot plot
    #
    ## Make figure
    fig, ax1 = plt.subplots(1,1)
    #fig.suptitle('Horizontally stacked subplots')
    ## SF plot
    ## scatter the points
    #ax1.scatter(SF_members[0],SF_members[1],s=(SF_members[2])*10,c='g', marker='o',alpha=0.4, linewidths = 0, label='Member')
    #ax1.scatter(SF_fields[0],SF_fields[1],s=(SF_fields[2])*10,c='b', marker='o',alpha=0.4, linewidths = 0, label='Field')
    #
    #
    ## plot spec subsample using [ 'sub'==1 (phot+spec subsample) & 'sub'!=4 (stars) ] as identifier for loop
    mem_scatter = []               # list del_z values for caluclating stats later
    outlier_scatter = []
    SF_scatter = []
    Q_scatter = []
    stars = 0             # count total # of outliers
    count_mem = 0               # count the rest
    count_outlier = 0 
    for counter in range(len(master_cat)):
        if master_cat['sub'][counter] == 1:    # sub=1 for (spec+phot) subsample
            if master_cat['type'][counter] == 3:    # type=3 for outliers
                outz = ax1.scatter(master_cat['z_spec'][counter],master_cat['z_peak'][counter],c='r',marker='v',alpha=0.8,linewidths = 0)
                outlier_scatter.append(np.abs(master_cat['del_z'][counter]))
                count_outlier+=1
                #print('Outlier: %s'%master_cat['type'][counter])
            elif master_cat['type'][counter] !=0 :
                memz = ax1.scatter(master_cat['z_spec'][counter],master_cat['z_peak'][counter],c='b', marker='^',alpha=0.5,linewidths = 0)
                count_mem+=1
                if master_cat['type'][counter] == 1:         # type=1 for SF
                    SF_scatter.append(np.abs(master_cat['del_z'][counter]))
                elif master_cat['type'][counter] == 2:       # type=1 for Q
                    Q_scatter.append(np.abs(master_cat['del_z'][counter]))
                mem_scatter.append(np.abs(master_cat['del_z'][counter]))
                #print('Non-outlier: %s'%master_cat['type'][counter])
            else:
                if master_cat['type'][counter] == 0:      #type=0 for stars
                    stars+=1
    #
    ## compute some stats
    SF_scatter = np.array(SF_scatter)
    Q_scatter = np.array(Q_scatter)
    mem_scatter = np.array(mem_scatter)
    outlier_scatter = np.array(outlier_scatter)             # convert to array for math operations
    abs_scatter = mem_scatter
    mean_abs_delz = np.round(np.array([np.median(mem_scatter),np.mean(outlier_scatter),np.mean(abs_scatter)]),decimals=2)   # compute mean of |del_z|; [SF, Q, total]
    std_abs_delz = np.round(np.array([np.std(mem_scatter),np.std(outlier_scatter),np.std(abs_scatter)]),decimals=2)   # compute std dev of |del_z|; [SF, Q, total]
    outlier_fraction = np.round((len(outlier_scatter) / (len(outlier_scatter)+len(mem_scatter)))*100,decimals=2)
    #
    #
    ### MAY NEED TO EDIT ### diag_flag_1
    ##  compare outliers & |del_z| stats identified in master_data_7_final.py with those identified above
    diag_flag_1 = 1             # 0=off (don't display diagnostic); 1=on (display diagnostic table)
    #
    if diag_flag_1 == 1:
        print('Data preparation file finds the following\nOUTLIERS total: %s' % np.sum(outliers),'\nOutlier fraction: %s' % (np.sum(outliers)/np.sum(both)),'\n|del_z| median: %s'%delz_median,'\n|del_z| scatter: %s\n'%delz_scatter)
        print('This plotting program found the following\nOUTLIERS total: %s' %len(outlier_scatter),'\nOutlier fraction: %s' %outlier_fraction,'\n|del_z| median: %s'%mean_abs_delz[0],'\n|del_z| scatter (excluding outliers): %s\n'%std_abs_delz[0])
        #print('FULL MEAN: %s'%np.mean(full_scatter),'\nFULL SCATTER: %s'%np.std(full_scatter))
        print('DIFFERENCES\nOutlier total: %s'%(np.sum(outliers)-count_outlier),'\nOutlier fraction: %s'%((np.sum(outliers)/np.sum(both))-outlier_fraction))#,'\n|del_z| mean: %s'%(delz_mean-mean_abs_delz[0]),'\n|del_z| scatter: %s\n'%(delz_scatter-std_abs_delz[0]))
        print('# of stars: %s'%stars,'\nNon-outlier count: %s'%count_outlier,'\nNon-outlier count: %s'%count_mem)
        print('SF scatter: %s'%np.std(SF_scatter),'\nQ scatter: %s'%np.std(Q_scatter))
    #
    ## derive the bound which separate outliers in terms of plotting arrays for X and Y
    x_bound = np.arange(0.0,1.601,0.001)
    #
    y_upper_bound, y_lower_bound = delz_bound(x_bound)
    #
    ## construct a string to plot the outlier fraction & std dev ('scatter')
    string = 'Median |$\Delta$z|: %.3f'%mean_abs_delz[0]+'\n$\sigma_{z}$ = %.3f'%std_abs_delz[0]+'\nOutliers: ~%.2f'%outlier_fraction+'%'
    ## add text to plot
    ax1.text(0.15,1.085,string)
    ax1.plot([0,2],[0,2],'--k', linewidth=1)
    ax1.plot(x_bound,y_upper_bound,':k', linewidth=1)
    ax1.plot(x_bound,y_lower_bound,':k', linewidth=1)
    ax1.set_xlabel("$z_{spec}$")
    ax1.set_xscale('linear')
    ax1.set_xlim(0,1.5)
    ax1.set_ylabel("$z_{phot}$")
    ax1.set_yscale('linear')
    ax1.set_ylim(0,1.5)
    #ax1.set_title("$z_{spec} vs. z_{phot}$")
    ax1.legend((memz,outz),('Spec. subsample: %i'%count_mem,'Outliers: %i'%count_outlier),scatterpoints=1,loc='upper left', frameon=False)
    #plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
    ax1.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on')
    ax1.minorticks_on()
    plt.show()
    #
#
#
#
#
#
## SECTION (2): del_z(spec) vs del_z(phot)
#
if plot_flag_2 == 1:
    #
    SF_counting_array = np.array([0]*4)
    Q_counting_array = np.array([0]*4)
    SF_members = [ [], [], [] ]                      # [ [z_clusterspec (x)], [z_clusterphot (y)], [mass (size)]]
    SF_fields = [ [], [], [] ]
    SF_poss = [ [], [], [] ]
    SF_negs = [ [], [], [] ]
    Q_members = [ [], [], [] ]                      # [ [z_clusterspec (x)], [z_clusterphot (y)], [mass (size)] ]
    Q_fields = [ [], [], [] ]
    Q_poss = [ [], [], [] ]
    Q_negs = [ [], [], [] ]
    #
    for counter in range(len(master_cat)):
        if master_cat['sub'][counter] == 1:                # only look at spec. sub-sample
            if master_cat[counter]['member'] == 0:                                     # member=0: cluster members
                if master_cat['type'][counter] == 1:       # type=1: SF
                    SF_members[0].append(master_cat[counter]['z_clusterspec'])
                    SF_members[1].append(master_cat[counter]['z_clusterphot'])
                    SF_members[2].append(master_cat[counter]['lmass'])
                    SF_counting_array[0]+=1
                elif master_cat['type'][counter] == 2:       # type=2: Q
                    Q_members[0].append(master_cat[counter]['z_clusterspec'])
                    Q_members[1].append(master_cat[counter]['z_clusterphot'])
                    Q_members[2].append(master_cat[counter]['lmass'])
                    Q_counting_array[0]+=1
            elif master_cat[counter]['member'] == 1:                                   # member=1: field
                if master_cat['type'][counter] == 1:       # type=1: SF
                    SF_fields[0].append(master_cat[counter]['z_clusterspec'])
                    SF_fields[1].append(master_cat[counter]['z_clusterphot'])
                    SF_fields[2].append(master_cat[counter]['lmass'])
                    SF_counting_array[1]+=1
                elif master_cat['type'][counter] == 2:       # type=2: Q
                    Q_fields[0].append(master_cat[counter]['z_clusterspec'])
                    Q_fields[1].append(master_cat[counter]['z_clusterphot'])
                    Q_fields[2].append(master_cat[counter]['lmass'])                    
                    Q_counting_array[1]+=1
            elif master_cat[counter]['member'] == 2:                                   # member=2: false pos
                if master_cat['type'][counter] == 1:       # type=1: SF
                    SF_poss[0].append(master_cat[counter]['z_clusterspec'])
                    SF_poss[1].append(master_cat[counter]['z_clusterphot'])
                    SF_poss[2].append(master_cat[counter]['lmass'])
                    SF_counting_array[2]+=1
                elif master_cat['type'][counter] == 2:       # type=2: Q
                    Q_poss[0].append(master_cat[counter]['z_clusterspec'])
                    Q_poss[1].append(master_cat[counter]['z_clusterphot'])
                    Q_poss[2].append(master_cat[counter]['lmass'])
                    Q_counting_array[2]+=1
            elif master_cat[counter]['member'] == 3:                                   # member=3: false neg
                if master_cat['type'][counter] == 1:       # type=1: SF
                    SF_negs[0].append(master_cat[counter]['z_clusterspec'])
                    SF_negs[1].append(master_cat[counter]['z_clusterphot'])
                    SF_negs[2].append(master_cat[counter]['lmass'])
                    SF_counting_array[3]+=1
                elif master_cat['type'][counter] == 2:       # type=2: Q
                    Q_negs[0].append(master_cat[counter]['z_clusterspec'])
                    Q_negs[1].append(master_cat[counter]['z_clusterphot'])
                    Q_negs[2].append(master_cat[counter]['lmass'])
                    Q_counting_array[3]+=1
    #
    #
    ## convert to arrays
    SF_members = np.array(SF_members)
    SF_fields = np.array(SF_fields)
    SF_poss = np.array(SF_poss)
    SF_negs = np.array(SF_negs)
    Q_members = np.array(Q_members)
    Q_fields = np.array(Q_fields)
    Q_poss = np.array(Q_poss)
    Q_negs = np.array(Q_negs)
    #
    ## print information for user - diagnostic check
    if diag_flag_2 == 1:
        print('\n"master_data*.py" REPORTS:\nSF - members: %s'%np.sum(mem[0]),'; fields: %s'%np.sum(field[0]),'; false pos: %s'%np.sum(pos[0]),'; false neg: %s'%np.sum(neg[0]))
        print('Q - members: %s'%np.sum(mem[1]),'; fields: %s'%np.sum(field[1]),'; false pos: %s'%np.sum(pos[1]),'; false neg: %s'%np.sum(neg[1]))
        print('\n"master_zplots*.py" REPORTS:\nSF - members: %s'%SF_counting_array[0],'; fields: %s'%SF_counting_array[1],'; false pos: %s'%SF_counting_array[2],'; false neg: %s'%SF_counting_array[3])
        print('Q - members: %s'%Q_counting_array[0],'; fields: %s'%Q_counting_array[1],'; false pos: %s'%Q_counting_array[2],'; false neg: %s'%Q_counting_array[3])
    #    
    ## Make figure
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.subplots_adjust(wspace=0,hspace=0)
    ## SF plot
    ## scatter the points
    ax1.scatter(SF_members[0],SF_members[1],s=55,c='g', marker='P',alpha=0.4, linewidths = 0, label='Member')
    ax1.scatter(SF_fields[0],SF_fields[1],s=55,c='b', marker='P',alpha=0.4, linewidths = 0, label='Field')
    ax1.scatter(SF_poss[0],SF_poss[1],s=55,c='firebrick', marker='X',alpha=0.4, linewidths = 0, label='False pos.')
    ax1.scatter(SF_negs[0],SF_negs[1],s=55,c='darkviolet', marker='X',alpha=0.4, linewidths = 0, label='False neg.')
    # add z_cutoff limits
    ax1.plot([-0.5,1],[z_cutoff[1],z_cutoff[1]],':k', linewidth=1)  # horizontal cuts (phot)
    ax1.plot([-0.5,1],[-z_cutoff[1],-z_cutoff[1]],':k', linewidth=1)
    ax1.plot([-z_cutoff[0],-z_cutoff[0]],[-0.5,1],':k', linewidth=1)  #vertical cuts
    ax1.plot([z_cutoff[0],z_cutoff[0]],[-0.5,1],':k', linewidth=1)  
    # add details
    ax1.set_title('Star-forming')
    ax1.set_xlabel('$(z_{spec} - z_{cluster})/(1+z_{spec})$')
    ax1.set_xlim(-0.25,0.25)
    ax1.set_ylabel('$(z_{phot} - z_{cluster})/(1+z_{phot})$')
    ax1.set_ylim(-0.25,0.25)
    ax1.grid(b=False)
    ax1.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=False, labelleft=True, labeltop=False,labelbottom=True)
    ax1.minorticks_on()
    ax1.legend(loc='upper left', frameon=False)
    #
    ## now do the Q plot
    ## scatter the points
    ax2.scatter(Q_members[0],Q_members[1],s=55,c='g', marker='P',alpha=0.4, linewidths = 0)
    ax2.scatter(Q_fields[0],Q_fields[1],s=55,c='b', marker='P',alpha=0.4, linewidths = 0)
    ax2.scatter(Q_poss[0],Q_poss[1],s=55,c='firebrick', marker='X',alpha=0.4, linewidths = 0)
    ax2.scatter(Q_negs[0],Q_negs[1],s=55,c='darkviolet', marker='X',alpha=0.4, linewidths = 0)
    # add z_cutoff limits
    ax2.plot([-0.5,1],[z_cutoff[1],z_cutoff[1]],':k', linewidth=1)  # horizontal cuts (phot)
    ax2.plot([-0.5,1],[-z_cutoff[1],-z_cutoff[1]],':k', linewidth=1)
    ax2.plot([-z_cutoff[0],-z_cutoff[0]],[-0.5,1],':k', linewidth=1)  #vertical cuts
    ax2.plot([z_cutoff[0],z_cutoff[0]],[-0.5,1],':k', linewidth=1)  
    # add details
    ax2.set_title('Quiescent')
    ax2.set_xlabel('$(z_{spec} - z_{cluster})/(1+z_{spec})$')
    ax2.set_xlim(-0.25,0.25)
    ax2.set_ylabel('$(z_{phot} - z_{cluster})/(1+z_{phot})$')
    ax2.yaxis.set_label_position('right')
    ax2.set_ylim(-0.25,0.25)
    ax2.grid(b=False)
    ax2.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=True, labelleft=False, labeltop=False,labelbottom=True)
    ax2.minorticks_on()                        
    #
    #
#####    
#    
#   
#
#
#
# SECTION 3:
#
## UVJ Diagram: colour-colour plots; modify to add in photometric sub-sample
#
if plot_flag_3 ==1:
    #
    #
    ## TIME_FLAG_3 START
    #
    #
    if time_flag_3 == 1 or time_flag == 2:
        start_time = time.time()
    #
    #
    counting_array_uvj = np.array([0]*5)         # counts number of galaxies at each condition; diagnostic check
    SF_array = [ [], [], [] ]                # stores: [ [(V-J)], [(U-V)], [mass] ]
    Q_array = [ [], [], [] ]    
    total_array = [ [], [], [] ]
    field_array = [ [], [], [] ]
    ## cluster
    for counter in range(len(master_cat)):
        if master_cat['member'][counter] == 1:                   # field galaxies
            counting_array_uvj[0]+=1
            field_array[0].append(master_cat['vj'][counter])
            field_array[1].append(master_cat['uv'][counter])
            field_array[2].append(master_cat['lmass'][counter])
            total_array[0].append(master_cat['vj'][counter])
            total_array[1].append(master_cat['uv'][counter])
            total_array[2].append(master_cat['lmass'][counter])
        elif master_cat['member'][counter] == 0:                 # members
            if master_cat['type'][counter] == 1:
                counting_array_uvj[1]+=1                         # count for SF
                SF_array[0].append(master_cat['vj'][counter])
                SF_array[1].append(master_cat['uv'][counter])
                SF_array[2].append(master_cat['lmass'][counter])
                total_array[0].append(master_cat['vj'][counter])
                total_array[1].append(master_cat['uv'][counter])
                total_array[2].append(master_cat['lmass'][counter])
            elif master_cat['type'][counter]== 2:
                counting_array_uvj[2]+=1                         # count for Q
                Q_array[0].append(master_cat['vj'][counter])
                Q_array[1].append(master_cat['uv'][counter])
                Q_array[2].append(master_cat['lmass'][counter])
                total_array[0].append(master_cat['vj'][counter])
                total_array[1].append(master_cat['uv'][counter])
                total_array[2].append(master_cat['lmass'][counter])
            else:
                counting_array_uvj[3] +=1
        else:
            counting_array_uvj[4] +=1
    #
    ## field
    field_array_par = [ [], [], [] ]
    counting_array_par = np.array([0]*2)
    #
    for counter in range(len(master_cat_par)):
        if master_cat_par['member'][counter] == 1:                   # field galaxies
            counting_array_par[0]+=1
            field_array_par[0].append(master_cat_par['vj'][counter])
            field_array_par[1].append(master_cat_par['uv'][counter])
            field_array_par[2].append(master_cat_par['lmass'][counter])
        else:
            counting_array_par[1] +=1
    #
    ## convert lists to arrays
    SF_array = np.array(SF_array)
    Q_array = np.array(Q_array)
    total_array = np.array(total_array)
    field_array = np.array(field_array)
    field_array_par = np.array(field_array_par)
    #
    #
    ## The following section is for plotting points above a certain mass limit
    #    
    if plot_flag_4 ==1:
        total_array_hi_mass = [ [], [], [] ]
        field_array_hi_mass = [ [], [], [] ]
        field_array_hi_mass_par = [ [], [], [] ]
        #
        #
        mass_threshold = [8,range2[1]]                        # mass threshold in units of log_10(Msol)
        counting_array_hi_mass = np.array([0]*2)
        #
        for ii in range(len(total_array[0])):
            if total_array[2][ii] >= mass_threshold[0]:
                total_array_hi_mass[0].append(total_array[0][ii])
                total_array_hi_mass[1].append(total_array[1][ii])
                total_array_hi_mass[2].append(total_array[2][ii])
                counting_array_hi_mass[0]+=1
        for ii in range(len(field_array[0])):
            if field_array[2][ii] >= mass_threshold[0]:
                field_array_hi_mass[0].append(field_array[0][ii])
                field_array_hi_mass[1].append(field_array[1][ii])
                field_array_hi_mass[2].append(field_array[2][ii])
                counting_array_hi_mass[1]+=1
        for ii in range(len(field_array_par[0])):
            if field_array_par[2][ii] >= mass_threshold[0]:
                field_array_hi_mass_par[0].append(field_array_par[0][ii])
                field_array_hi_mass_par[1].append(field_array_par[1][ii])
                field_array_hi_mass_par[2].append(field_array_par[2][ii])
                counting_array_hi_mass[1]+=1
        #
        #
        ## convert lists to arrays
        total_array_hi_mass = np.array(total_array_hi_mass)
        field_array_hi_mass = np.array(field_array_hi_mass)
        field_array_hi_mass_par = np.array(field_array_hi_mass_par)
        #
        #
        # definitions for the axes
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        spacing = 0.005
        #
        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom + height + spacing, width, 0.2]
        rect_histy = [left + width + spacing, bottom, 0.2, height]
        #
        # start with a rectangular Figure
        plt.figure(figsize=(8, 8))
        #
        plt.suptitle('~.025 < z < ~0.6, $log_{10}$(M/M$_{\odot})$ > %s'%mass_threshold[0])
        #
        ax_scatter = plt.axes(rect_scatter)
        ax_scatter.tick_params(direction='in', top=True, right=True)
        ax_histx = plt.axes(rect_histx)
        ax_histx.tick_params(direction='in', labelbottom=False)
        ax_histy = plt.axes(rect_histy)
        ax_histy.tick_params(direction='in', labelleft=False)
        #
        # the scatter plot:\
        ax_scatter.scatter(field_array_hi_mass[0],field_array_hi_mass[1],s=100,c='darkgrey', marker='.', alpha=0.6, linewidths=0)
        ax_scatter.scatter(field_array_hi_mass_par[0],field_array_hi_mass_par[1],s=100,c='darkgrey', marker='.', alpha=0.6, linewidths=0)
        ax_scatter.scatter(total_array_hi_mass[0],total_array_hi_mass[1],s=100,c='slategrey', marker='.', alpha=0.6, linewidths=0)        
        #
        # now determine nice limits by hand:
        binwidth = 0.05
        lim = np.ceil(np.abs([total_array[0],total_array[1]]).max() / binwidth) * binwidth
        ax_scatter.set_xlim(-1.5,2)
        ax_scatter.set_ylim(0,2.5)
        #
        bins_x = np.arange(-1.5, 2, binwidth)
        bins_y = np.arange(0, 2.5, binwidth)
        ax_histx.hist(np.concatenate([total_array_hi_mass[0],field_array_hi_mass[0],field_array_hi_mass_par[0]]), bins=bins_x, color='slategrey')
        ax_histy.hist(np.concatenate([total_array_hi_mass[1],field_array_hi_mass[1],field_array_hi_mass_par[1]]), bins=bins_y, color='slategrey', orientation='horizontal')
        #
        ax_histx.set_xlim(ax_scatter.get_xlim())
        ax_histy.set_ylim(ax_scatter.get_ylim())
        #
        #ax_scatter.plot([-5,0.75],[1.3,1.3],'-k',[0.75,3],[1.3,3],'-k', linewidth=1) #overlay boundary cutoff for SF/Passive
        ax_scatter.set_xlabel('$(V-J)_{rest}$',fontsize=20)
        #ax_scatter.set_xlim(-1.5,2)
        ax_scatter.set_ylabel('$(U-V)_{rest}$', fontsize=20)
        #ax_scatter.set_ylim(0,2.5)
        ax_scatter.grid(b=False)
        ax_scatter.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=False, labelleft=True, labeltop=False)
        ax_scatter.minorticks_on()
        #
        ax_scatter.set_title('~0.25 < z < ~0.6')
        ax_scatter.text(-1.3,2.3,'Members: %s'%counting_array_hi_mass[0],c='k',fontsize=15)
        ax_scatter.text(-1.3,2.15,'Field (clu+par): %s'%(counting_array_hi_mass[1]),c='k',fontsize=15)
        #
        plt.show()
    #####
    elif plot_flag_4 ==2:
        #
        # definitions for the axes
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        spacing = 0.005
        #
        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom + height + spacing, width, 0.2]
        rect_histy = [left + width + spacing, bottom, 0.2, height]
        #
        # start with a rectangular Figure
        plt.figure(figsize=(8, 8))
        #
        plt.suptitle('~.025 < z < ~0.6')
        #
        ax_scatter = plt.axes(rect_scatter)
        ax_scatter.tick_params(direction='in', top=True, right=True)
        ax_histx = plt.axes(rect_histx)
        ax_histx.tick_params(direction='in', labelbottom=False)
        ax_histy = plt.axes(rect_histy)
        ax_histy.tick_params(direction='in', labelleft=False)
        #
        # the scatter plot:
        ax_scatter.scatter(field_array[0],field_array[1],s=100,c='darkgrey', marker='.', alpha=0.6, linewidths=0)
        ax_scatter.scatter(field_array_par[0],field_array_par[1],s=100,c='darkgrey', marker='.', alpha=0.6, linewidths=0)
        ax_scatter.scatter(total_array[0],total_array[1],s=100,c='slategrey', marker='.', alpha=0.6, linewidths=0)
        #
        # now determine nice limits by hand:
        binwidth = 0.05
        lim = np.ceil(np.abs([total_array[0],total_array[1]]).max() / binwidth) * binwidth
        ax_scatter.set_xlim(-1.5,2)
        ax_scatter.set_ylim(0,2.5)
        #
        bins_x = np.arange(-1.5, 2, binwidth)
        bins_y = np.arange(0, 2.5, binwidth)
        ax_histx.hist(np.concatenate([total_array[0],field_array[0],field_array_par[0]]), bins=bins_x, color='slategrey')
        ax_histy.hist(np.concatenate([total_array[1],field_array[1],field_array_par[1]]), bins=bins_y, color='slategrey', orientation='horizontal')
        #
        ax_histx.set_xlim(ax_scatter.get_xlim())
        ax_histy.set_ylim(ax_scatter.get_ylim())
        #
        #ax_scatter.plot([-5,0.75],[1.3,1.3],'-k',[0.75,3],[1.3,3],'-k', linewidth=1) #overlay boundary cutoff for SF/Passive
        ax_scatter.set_xlabel('$(V-J)_{rest}$',fontsize=20)
        #ax_scatter.set_xlim(-1.5,2)
        ax_scatter.set_ylabel('$(U-V)_{rest}$', fontsize=20)
        #ax_scatter.set_ylim(0,2.5)
        ax_scatter.grid(b=False)
        ax_scatter.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=False, labelleft=True, labeltop=False)
        ax_scatter.minorticks_on()
        #
        ax_scatter.set_title('~0.25 < z < ~0.6')
        ax_scatter.text(-1.3,2.3,'SF: %s'%counting_array_uvj[1],c='k',fontsize=15)
        ax_scatter.text(-1.3,2.15,'Q: %s'%counting_array_uvj[2],c='k',fontsize=15)
        ax_scatter.text(-1.3,2.0,'Field (clu+par): %s'%(counting_array_uvj[0]+counting_array_par[0]),c='k',fontsize=15)
        #
        plt.show()
        #
        #
    else:
        #
        fig, ax = plt.subplots()
        # print diagnostic information
        print('Field: %s'%counting_array[0],'\nSF members: %s'%counting_array[1],'\nQ members: %s'%counting_array[2],'\nOther than SF/Q (members): %s'%counting_array[3],'\nOther than member/field: %s'%counting_array[4])
        #
        c_map = plt.get_cmap("CMRmap")        # set colorbar map
        sm =  ScalarMappable(cmap=c_map)
        sm.set_array([])
        sm.set_clim([min(min(total_array[2]),min(field_array[2]),min(field_array_par[2])),max(max(total_array[2]),max(field_array[2]),max(field_array_par[2]))])
        cbar = fig.colorbar(sm, ax=ax)
        cbar.ax.set_title('$log_{10}$(M/M$_{\odot})$')
        #
        #
        ax.scatter(field_array[0],field_array[1],s=(field_array[2]*1.5),c='grey', marker='s', alpha=0.35, linewidths=0)
        ax.scatter(field_array_par[0],field_array_par[1],s=(field_array_par[2]*1.5),c='grey', marker='s', alpha=0.35, linewidths=0)
        ax.scatter(total_array[0],total_array[1],c=total_array[2],cmap= c_map, marker='.', alpha=0.8, linewidths=0)
        #ax.scatter(SF_array[0],SF_array[1],s=(SF_array[2]*1.5),cmap= c_map, marker='*', alpha=0.8, linewidths=0)
        #ax.scatter(Q_array[0],Q_array[1],s=(Q_array[2]*1.5),cmap= c_map, marker='.', alpha=0.4, linewidths=0)
        ax.plot([-5,0.75],[1.3,1.3],'-k',[0.75,3],[1.3,3],'-k', linewidth=1) #overlay boundary cutoff for SF/Passive
        ax.set_xlabel('$(V-J)_{rest}$',fontsize=20)
        ax.set_xlim(-1.5,2)
        ax.set_ylabel('$(U-V)_{rest}$', fontsize=20)
        ax.set_ylim(0,2.5)
        ax.grid(b=False)
        ax.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=False, labelleft=True)
        ax.minorticks_on()
        #
        ax.set_title('~0.25 < z < ~0.6')
        ax.text(-1.3,2.3,'SF: %s'%counting_array_uvj[1],c='k',fontsize=15)
        ax.text(-1.3,2.15,'Q: %s'%counting_array_uvj[2],c='k',fontsize=15)
        ax.text(-1.3,2.0,'Field (clu+par): %s'%(counting_array_uvj[0]+counting_array_par[0]),c='k',fontsize=15)
        #
        plt.show()
        #
    #
    #
## TIME_FLAG END
#
if time_flag_3 == 1 or time_flag == 2:
    print('\n"master_zplots*.py" Section 3 (UVJ diagram) complete.\n\nIt took: %s seconds produce this figure.\n\n' % (time.time() - start_time))
    #
#
#
#
#
## TIME_FLAG END
#
if time_flag == 1:
    print('Program "master_zplots_2_final.py" took: %s seconds to run.\n\n' % (time.time() - start_time))
#
print('\n\n"master_zplots_2_final.py"  terminated successfully.\n')
#
#
#
#
#                        
###### PROGRAM END ######