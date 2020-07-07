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
## function definitions
#
#
BCG_threshold = 12
##
num_small_BCG = np.array([0]*6)
num_BCG = np.array([0]*6)
BCG_spec = np.array([[0]*6]*2)       # row1=SF;  row2=Q; cols=clusters
BCG_phot = np.array([[0]*6]*2)
BCG_SF = np.array([0]*6)
BCG_Q =  np.array([0]*6)
BCG_delz = [ [], [], [], [], [], [] ]
small_BCG_delz = [ [], [], [], [], [], [] ]
#
#
for counter in range(len(master_cat)):
    if master_cat['lmass'][counter] >= BCG_threshold:
    #if master_cat['flag_F160W'][counter] == 4:
        #master_cat['member'][counter] = 5              # member=5 for BCGs so they don't contaminate the "member=0" cluster member sample
        for BCG in range(len(num_BCG)):
            if master_cat['cluster'][counter] == (BCG+1):    # track massive galaxies by cluster
                num_BCG[BCG]+=1
                BCG_delz[BCG].append(master_cat['z_clusterphot'][counter])
                if master_cat['type'][counter] == 1:         # SF sample
                    BCG_SF[BCG]+=1    
                    if master_cat['sub'][counter] == 1:         # spec subsample
                        BCG_spec[0][BCG]+=1
                    elif master_cat['sub'][counter] == 2:         # spec subsample
                        BCG_phot[0][BCG]+=1
                elif master_cat['type'][counter] == 2:         # Q sample
                    BCG_Q[BCG]+=1
                    if master_cat['sub'][counter] == 1:         # spec subsample
                        BCG_spec[1][BCG]+=1
                    elif master_cat['sub'][counter] == 2:         # spec subsample
                        BCG_phot[1][BCG]+=1
    elif master_cat['lmass'][counter] >= 11 and master_cat['lmass'][counter] < BCG_threshold:      # track "small bCGs"
        for BCG in range(len(num_BCG)):
            if master_cat['cluster'][counter] == (BCG+1):    # track massive galaxies by cluster
                small_BCG_delz[BCG].append(master_cat['z_clusterphot'][counter])
                num_small_BCG[BCG]+=1
    #############
#####
#
#bins_phot = [-0.5,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.10,0.15,0.20,0.25,0.3,0.5,1.0]
#bins_phot = [0.0,0.05,0.10,0.15,0.20,0.25,0.3,0.5,1.0]
bins_phot = np.arange(0.0,0.505,0.01)
#
## collape for plotting purposes
BCG_delz_plot = []
for ii in range(len(BCG_delz)):
    for jj in range(len(BCG_delz[ii])):
        BCG_delz_plot.append(np.abs(BCG_delz[ii][jj]))
    #
## Now look through list BCG_delz and pick out all the galaxies w/ BAD redshift estimates (quantify them), and then get some stats on the remaining galaxies
num_bad_BCG = 0
good_BCG = []
#
for ii in range(len(BCG_delz_plot)):
    if np.abs(BCG_delz_plot[ii]) == 99.0:
        num_bad_BCG+=1
    else:
        good_BCG.append(BCG_delz_plot[ii])
        #
#
total_BCG = len(BCG_delz_plot)
#
string = 'z=-99.0: %s'%num_bad_BCG+'/%s'%total_BCG
print('z=-99 (not in Parent sample): %s'%num_bad_BCG,'\nOther (in Parent sample)%s'%len(good_BCG))
#
## Visualize del_z distribution
# PHOT
plt.figure()
n, bins, patches = plt.hist(x=np.abs(good_BCG),bins=bins_phot,color='deepskyblue',edgecolor='steelblue',alpha=0.7, rwidth=1)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('$z_{phot}$ - $z_{cluster}$ / 1 + $z_{phot}$',fontsize=12)
plt.ylabel('# count',fontsize=12)
plt.title('"|${\Delta}$z|$_{phot}$" bCG',fontsize=15)
ax.text(0.2,3,string,fontsize=9)
plt.show()
# SPEC
#plt.figure()
#n, bins, patches = plt.hist(x=BCG_spec_plot,bins=bins_phot,color='darkcyan',alpha=0.7, rwidth=1)
#plt.grid(axis='y', alpha=0.75)
#plt.xlabel('$z_{spec}$ - $z_{cluster}$ / 1 + $z_{spec}$',fontsize=12)
#plt.ylabel('# count',fontsize=12)
#plt.title('abs(DEL_Z SPEC)',fontsize=15)
#plt.show()
#
#
#print('# phot: %s'%len(delz_phot_plot),'\n# spec: %s'%len(delz_spec_plot),'\nTotal: %s'%(len(delz_phot_plot)+len(delz_spec_plot)))
#
## Visualize by cluster
#
#
nrows=2
ncols=3
#plot_titles = ['SF','Q']
#
fig, axs = plt.subplots(nrows,ncols,sharex=True,sharey=True,tight_layout=True)
#fig.suptitle('abs(DEL_Z PHOT)',fontsize=12)
#
counter = 0
#
for jj in range(nrows):
    for kk in range(ncols):
        if counter < len(BCG_delz_plot):
        #
            ax = axs[jj][kk]
            #
            ax.hist(np.abs(BCG_delz[counter]),bins=bins_phot,color='deepskyblue',edgecolor='steelblue',alpha=0.7)
            ax.grid(axis='y', alpha=0.75)
            ax.set_xlim([min(bins_phot),max(bins_phot)])
            ax.title.set_text(cluster_names[kk])
            #ax.text(0.05,2500,cluster_names[counter],fontsize=9)
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