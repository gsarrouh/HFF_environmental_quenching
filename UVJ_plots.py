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
### Section (1) aggregate cluster/field galaxies into arrays for plotting
### Section (2) add PLOT_FLAG_1: create UVJ diagram per Shipley et al 2018 to segregate between type (SF/Q - i.e. a "colour-colour" plot)
### Section (3) add PLOT_FLAG_2: diagnostic plots: 0=off;  1=on, plot everything above/below "mass_threshold";  2=on, plot everything (all masses)
#
### PROGRAM END
#
#
#
###################     PROGRAM START
#
#  
# import modules
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
import time
#
## definitions
#
## 
#
## MAY NEED TO EDIT
## set mass threshold for plot_flag_2
mass_threshold = [9,range2[1]]                        # mass threshold in units of log_10(Msol)
#
## SECTION (0): set PLOT_FLAGS    # 0=off (i.e. don't make plot), 1=on (i.e. make plot)
### MAY NEED TO EDIT ### plot_flag_*/time_flag_*/diag_flag_*
#
plot_flag_1 = 1           # Fig.1 - UVJ diagram
plot_flag_2 = 1           # Fig.2-4 - diagnostic plots in greyscale;  0=off;  1=on, plot above log_10(Msol) >/< threshold;  2=on, plot everything; 
#
lines_flag = 0            # include boundary line b/w SF/Q populations in UVJ plot
#
diag_flag_1 = 0           # Section 1: UVJ diagram
#
#
#
# SECTION 1: aggregate lists
#
#
#
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
#
## UVJ Diagram: colour-colour plots; modify to add in photometric sub-sample
#
if (plot_flag_1 == 1 and project_plot_flag ==2) or project_plot_flag == 1:
    if project_plot_flag == 0:
        pass
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
        #
        if lines_flag == 1:
            ax.plot([-5,0.75],[1.3,1.3],'-k',[0.75,3],[1.3,3],'-k', linewidth=1) #overlay boundary cutoff for SF/Passive
        #
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
    #
#
## The following section is for plotting points above a certain mass limit
#    
if (plot_flag_2 == 1 and project_plot_flag ==2) or project_plot_flag == 1:
    if project_plot_flag == 0:
        pass
    else:
        #
        ## create plot for all galaxies ABOVE "mass_threshold"
        total_array_hi_mass = [ [], [], [] ]
        field_array_hi_mass = [ [], [], [] ]
        field_array_hi_mass_par = [ [], [], [] ]
        #
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
        if lines_flag == 1:
            ax_scatter.plot([-5,0.75],[1.3,1.3],'-k',[0.75,3],[1.3,3],'-k', linewidth=1) #overlay boundary cutoff for SF/Passive
        #
        ax_scatter.set_xlabel('$(V-J)_{rest}$',fontsize=20)
        #ax_scatter.set_xlim(-1.5,2)
        ax_scatter.set_ylabel('$(U-V)_{rest}$', fontsize=20)
        #ax_scatter.set_ylim(0,2.5)
        ax_scatter.grid(b=False)
        ax_scatter.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=False, labelleft=True, labeltop=False)
        ax_scatter.minorticks_on()
        #
        ax_scatter.set_title('~0.25 < z < ~0.6')
        ax_scatter.text(-1.3,2.3,'Members: %s'%counting_array_hi_mass[0],c='slategrey',fontsize=15)
        ax_scatter.text(-1.3,2.15,'Field (clu+par): %s'%(counting_array_hi_mass[1]),c='darkgrey',fontsize=15)
        #
        plt.show()
        #
        #
        ## now do the same for all objects BELOW "mass_threshold"
        #
        total_array_lo_mass = [ [], [], [] ]
        field_array_lo_mass = [ [], [], [] ]
        field_array_lo_mass_par = [ [], [], [] ]
        #
        counting_array_lo_mass = np.array([0]*2)
        #
        for ii in range(len(total_array[0])):
            if total_array[2][ii] <= mass_threshold[0]:
                total_array_lo_mass[0].append(total_array[0][ii])
                total_array_lo_mass[1].append(total_array[1][ii])
                total_array_lo_mass[2].append(total_array[2][ii])
                counting_array_lo_mass[0]+=1
        for ii in range(len(field_array[0])):
            if field_array[2][ii] <= mass_threshold[0]:
                field_array_lo_mass[0].append(field_array[0][ii])
                field_array_lo_mass[1].append(field_array[1][ii])
                field_array_lo_mass[2].append(field_array[2][ii])
                counting_array_lo_mass[1]+=1
        for ii in range(len(field_array_par[0])):
            if field_array_par[2][ii] <= mass_threshold[0]:
                field_array_lo_mass_par[0].append(field_array_par[0][ii])
                field_array_lo_mass_par[1].append(field_array_par[1][ii])
                field_array_lo_mass_par[2].append(field_array_par[2][ii])
                counting_array_lo_mass[1]+=1
        #
        #
        ## convert lists to arrays
        total_array_lo_mass = np.array(total_array_lo_mass)
        field_array_lo_mass = np.array(field_array_lo_mass)
        field_array_lo_mass_par = np.array(field_array_lo_mass_par)
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
        plt.suptitle('~.025 < z < ~0.6, $log_{10}$(M/M$_{\odot})$ < %s'%mass_threshold[0])
        #
        ax_scatter = plt.axes(rect_scatter)
        ax_scatter.tick_params(direction='in', top=True, right=True)
        ax_histx = plt.axes(rect_histx)
        ax_histx.tick_params(direction='in', labelbottom=False)
        ax_histy = plt.axes(rect_histy)
        ax_histy.tick_params(direction='in', labelleft=False)
        #
        # the scatter plot:\
        ax_scatter.scatter(field_array_lo_mass[0],field_array_lo_mass[1],s=100,c='darkgrey', marker='.', alpha=0.6, linewidths=0)
        ax_scatter.scatter(field_array_lo_mass_par[0],field_array_lo_mass_par[1],s=100,c='darkgrey', marker='.', alpha=0.6, linewidths=0)
        ax_scatter.scatter(total_array_lo_mass[0],total_array_lo_mass[1],s=100,c='slategrey', marker='.', alpha=0.6, linewidths=0)        
        #
        # now determine nice limits by hand:
        binwidth = 0.05
        lim = np.ceil(np.abs([total_array[0],total_array[1]]).max() / binwidth) * binwidth
        ax_scatter.set_xlim(-1.5,2)
        ax_scatter.set_ylim(0,2.5)
        #
        bins_x = np.arange(-1.5, 2, binwidth)
        bins_y = np.arange(0, 2.5, binwidth)
        ax_histx.hist(np.concatenate([total_array_lo_mass[0],field_array_lo_mass[0],field_array_lo_mass_par[0]]), bins=bins_x, color='slategrey')
        ax_histy.hist(np.concatenate([total_array_lo_mass[1],field_array_lo_mass[1],field_array_lo_mass_par[1]]), bins=bins_y, color='slategrey', orientation='horizontal')
        #
        ax_histx.set_xlim(ax_scatter.get_xlim())
        ax_histy.set_ylim(ax_scatter.get_ylim())
        #
        if lines_flag == 1:
            ax_scatter.plot([-5,0.75],[1.3,1.3],'-k',[0.75,3],[1.3,3],'-k', linewidth=1) #overlay boundary cutoff for SF/Passive
        #
        ax_scatter.set_xlabel('$(V-J)_{rest}$',fontsize=20)
        #ax_scatter.set_xlim(-1.5,2)
        ax_scatter.set_ylabel('$(U-V)_{rest}$', fontsize=20)
        #ax_scatter.set_ylim(0,2.5)
        ax_scatter.grid(b=False)
        ax_scatter.tick_params(axis='both', which='both',direction='in',color='k',top=True,right=True,labelright=False, labelleft=True, labeltop=False)
        ax_scatter.minorticks_on()
        #
        ax_scatter.set_title('~0.25 < z < ~0.6')
        ax_scatter.text(-1.3,2.3,'Members: %s'%counting_array_lo_mass[0],c='slategrey',fontsize=15)
        ax_scatter.text(-1.3,2.15,'Field (clu+par): %s'%(counting_array_lo_mass[1]),c='darkgrey',fontsize=15)
        #
        plt.show()
        #
#####
elif (plot_flag_2 == 2 and project_plot_flag ==2) or project_plot_flag == 1:
    if project_plot_flag == 0:
        pass
    else:
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
        ax_scatter.text(-1.3,2.3,'SF: %s'%counting_array_uvj[1],c='slategrey',fontsize=15)
        ax_scatter.text(-1.3,2.15,'Q: %s'%counting_array_uvj[2],c='darkgrey',fontsize=15)
        ax_scatter.text(-1.3,2.0,'Field (clu+par): %s'%(counting_array_uvj[0]+counting_array_par[0]),c='k',fontsize=15)
        #
        plt.show()
        #
        #
        
    
#
#
#
print('\n\n"UVJ_plots.py"  terminated successfully.\n')
#
#
#
#
#                        
###### PROGRAM END ######