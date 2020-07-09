# Created Sat Jul 04 16:42:27 2020
#
## This is a diagnostic file to investigate the distribution of del_z spec(phot) = z_spec(z_phot) - z_cluster / 1+z_spec(z_phot). It will produce histograms of the entire population of del_z's calculated separately for phot. & spec. 
#
###################     PROGRAM START     ###################
#
## modules
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#
## PLOT FLAG: 0=off, supresses output of figure;  1=on, enables plotting of figure
plot_flag = 0
#
## function definitions
#
## define a function to compute the mid-points of the bins from a histogram
def midbins(bins):
    size = len(bins)-1
    x_midbins = np.empty([size,1],dtype='float64')
    for x in range(size):
        x_midbins[x] = (bins[x] + bins[(x+1)])/2
    return x_midbins
#
#
## initialize storage arrays
delz_phot = [[],[],[],[],[],[]]
delz_spec = [[],[],[],[],[],[]]
count_outliers = 0
#
for counter in range(len(master_cat)):
    if master_cat['sub'][counter] == 2:         # identifies PHOT-ONLY subsample
        for ii in range(len(delz_phot)):
            if master_cat['cluster'][counter] == (ii+1):
                delz_phot[ii].append(master_cat['z_clusterphot'][counter])
    elif master_cat['sub'][counter] == 1:       # identifies SPEC-ONLY subsample
        if master_cat['type'][counter] == 3:
            count_outliers+=1
        else:
            for ii in range(len(delz_phot)):
                if master_cat['cluster'][counter] == (ii+1):
                    delz_spec[ii].append(master_cat['z_clusterspec'][counter])
#########
bins_phot1 = np.arange(0.0,0.95,0.05)
bins_phot2 = np.arange(0.0,0.11,0.01)
#bins_phot = [-0.5,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.10,0.15,0.20,0.25,0.3,0.5,1.0]

#
## collape for plotting purposes
delz_phot_plot = []
delz_spec_plot = []
for ii in range(len(delz_phot)):
    for jj in range(len(delz_phot[ii])):
        delz_phot_plot.append(delz_phot[ii][jj])
    for kk in range(len(delz_spec[ii])):    
        delz_spec_plot.append(delz_spec[ii][kk])
    #
####
#
## Visualize del_z distribution
#
## SPEC
plt.figure()
n, bins, patches = plt.hist(x=np.abs(delz_spec_plot),bins=bins_phot1,color='deepskyblue',edgecolor='steelblue',alpha=0.7,rwidth=0.95,log=True)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('$z_{spec}$ - $z_{cluster}$ / 1 + $z_{spec}$',fontsize=12)
plt.ylabel('# count',fontsize=12)
plt.title("|${\Delta}$z|$_{spec}$",fontsize=15)
plt.show()
#
## PHOT
plt.figure()
n, bins, patches = plt.hist(x=np.abs(delz_phot_plot),bins=bins_phot1,color='deepskyblue',edgecolor='steelblue',alpha=0.7,rwidth=0.95,log=True)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('$z_{phot}$ - $z_{cluster}$ / 1 + $z_{phot}$',fontsize=12)
plt.ylabel('# count',fontsize=12)
plt.title("|${\Delta}$z|$_{phot}$",fontsize=15)
plt.show()
#
#
if plot_flag == 1:
    ## PHOT zoomed in on - < del_z < 0.1
    plt.figure()
    n, bins, patches = plt.hist(x=np.abs(delz_phot_plot),bins=bins_phot2,color='deepskyblue',edgecolor='steelblue',alpha=0.7,rwidth=0.95,log=False)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('$z_{phot}$ - $z_{cluster}$ / 1 + $z_{phot}$',fontsize=12)
    plt.ylabel('# count',fontsize=12)
    plt.title("|${\Delta}$z|$_{phot}$",fontsize=15)
    plt.show()


#
print('# phot: %s'%len(delz_phot_plot),'\n# spec: %s'%len(delz_spec_plot),'\nTotal: %s'%(len(delz_phot_plot)+len(delz_spec_plot)))
#
if plot_flag == 1:
    ## Visualize by cluster
    #
    nrows=2
    ncols=3
    #
    fig, axs = plt.subplots(nrows,ncols,sharex=True,sharey=True,tight_layout=True)
    fig.subplots_adjust(wspace=0,hspace=0)
    #fig.suptitle('abs(DEL_Z PHOT)',fontsize=12)
    #
    counter = 0
    #
    for jj in range(nrows):
        for kk in range(ncols):
            if counter < len(delz_phot):
            #
                ax = axs[jj][kk]
                #
                ax.hist(np.abs(delz_phot[counter]),bins=bins_phot1,color='deepskyblue',edgecolor='steelblue',alpha=0.7)
                ax.grid(axis='y', alpha=0.75)
                ax.set_xlim([min(bins_phot1),max(bins_phot1)])
                ax.text(0.05,700,cluster_names[counter],fontsize=9)
                counter+=1
                if jj == 1:      # put xlabels on 1st row
                    ax.set_xlabel('$z_{phot}$ - $z_{cluster}$ / 1 + $z_{phot}$',fontsize=10)
                if kk == 0:
                    ax.set_ylabel('# count',fontsize=10)
    #####            
    #            
    #            
    plt.show()
    #
#
#
#
#
#
###################     PROGRAM END     ###################