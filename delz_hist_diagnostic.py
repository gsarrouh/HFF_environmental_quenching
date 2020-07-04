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
import astropy
from astropy.table import Table
from astropy.table import Column
from scipy.optimize import curve_fit
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
delz_phot = []
delz_spec = []
#
for counter in range(len(master_cat)):
    if master_cat['sub'][counter] == 2:
        delz_phot.append(np.abs(master_cat['z_clusterphot'][counter]))
    elif master_cat['sub'][counter] == 1:
        delz_spec.append(np.abs(master_cat['z_clusterspec'][counter]))
#########
bins_phot = np.arange(0,1.1,0.1)
#
## Visualize
# PHOT
plt.figure()
n, bins, patches = plt.hist(x=delz_phot, bins=bins_phot, color='#0504aa',alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('$z_{phot}$ - $z_{cluster}$ / 1 + $z_{phot}$',fontsize=12)
plt.ylabel('# count',fontsize=12)
plt.title('abs(DEL_Z PHOT)',fontsize=15)
plt.show()
# PHOT
plt.figure()
n, bins, patches = plt.hist(x=delz_spec,bins=bins_phot,color='#0504aa',alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('$z_{spec}$ - $z_{cluster}$ / 1 + $z_{spec}$',fontsize=12)
plt.ylabel('# count',fontsize=12)
plt.title('abs(DEL_Z SPEC)',fontsize=15)
plt.show()
#
#
print('# phot: %s'%len(delz_phot),'\n# spec: %s'%len(delz_spec),'\nTotal: %s'%(len(delz_phot)+len(delz_spec)))
#
#
#
###################     PROGRAM END     ###################