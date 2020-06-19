#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 09:17:20 2017

@author: gsarrouh
"""
###This program creates plots for cluster macs0416, containing:
###photometric v specroscopic redshift, delta-z, & colour-colour plots

#import numpy as np
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import astropy
#from astropy.table import Table
#from astropy.table import Column

plt.clf()
#make z_spec vs z_phot plot
zees = plt.figure(num=1)
total_zspec = Column(np.concatenate((SF_spec['z_spec'],SF_field['z_spec'],SF_fpos['z_spec'],SF_fneg['z_spec'],Q_spec['z_spec'],Q_field['z_spec'],Q_fpos['z_spec'],Q_fneg['z_spec'])))
total_zphot = Column(np.concatenate((SF_spec['z_peak'],SF_field['z_peak'],SF_fpos['z_peak'],SF_fneg['z_peak'],Q_spec['z_peak'],Q_field['z_peak'],Q_fpos['z_peak'],Q_fneg['z_peak'])))
out_spec = Column(np.concatenate((SF_spec_out['z_spec'],Q_spec_out['z_spec'])))
out_phot = Column(np.concatenate((SF_spec_out['z_peak'],Q_spec_out['z_peak'])))

zvz = zees.add_subplot(111)
memz = zvz.scatter(total_zspec,total_zphot, c='b', marker='*', linewidths = 0)
outz = zvz.scatter(out_spec,out_spec, c='r', marker='x', linewidths = 0)
lin1 = zvz.plot([0,1],[0,1],':k', linewidth=1)
plt.xlabel('z_spec')
plt.xscale('linear')
plt.xlim(0,1)
plt.ylabel('z_photo')
plt.yscale('linear')
plt.ylim(0,1)
plt.title('z_spec vs. z_photo')
plt.legend((memz,outz),('spec population','outliers'),scatterpoints=1,loc='upper left')
#plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
#plt.tick_params(axis='y', direction='in',color='k')
#plt.show()

#cluster_spec vs cluster_phot plot, side-by-side
delzees = plt.figure(num=2)

delzSF = delzees.add_subplot(121) #SF plot
SFz = delzSF.scatter(SF_spec['z_clusterspec'],SF_spec['z_clusterphot'], c='g', marker='+', linewidths = 0)
SFz_field = delzSF.scatter(SF_field['z_clusterspec'],SF_field['z_clusterphot'], c='b', marker='+', linewidths = 0)
SFz_fpos = delzSF.scatter(SF_fpos['z_clusterspec'],SF_fpos['z_clusterphot'], c='r', marker='x', linewidths = 0)
SFz_fneg = delzSF.scatter(SF_fneg['z_clusterspec'],SF_fneg['z_clusterphot'], c='m', marker='x', linewidths = 0)
plt.xlabel('(z_spec - z_cluster)/(1+z_spec)')
plt.xscale('linear')
plt.xlim(-0.125,0.125)
plt.tick_params(axis='both', which='both',direction='in',color='k')
plt.ylabel('(z_phot - z_cluster)/(1+z_phot)')
delzSF.yaxis.tick_left()
plt.yscale('linear')
plt.ylim(-0.15,0.15)
lin = delzSF.plot([-0.5,1],[0.05,0.05],':k', linewidth=1)  # horizontal cuts
lin = delzSF.plot([-0.5,1],[-0.05,-0.05],':k', linewidth=1)
lin = delzSF.plot([-0.01,-0.01],[-0.5,1],':k', linewidth=1)  #vertical cuts
lin = delzSF.plot([0.01,0.01],[-0.5,1],':k', linewidth=1)  
plt.legend((SFz,SFz_field,SFz_fpos,SFz_fneg),('secure member','secure field','false positive','false negative'),scatterpoints=1,loc='upper left', fontsize=10)
plt.title('Star-Forming')
#plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
plt.tick_params(axis='y', which='both', direction='in',color='k')

delzPass = delzees.add_subplot(122) #Passive plot
Passz = delzPass.scatter(Q_spec['z_clusterspec'],Q_spec['z_clusterphot'], c='g', marker='+', linewidths = 0)
Passz_field = delzPass.scatter(Q_field['z_clusterspec'],Q_field['z_clusterphot'], c='b', marker='+', linewidths = 0)
Passz_fpos = delzPass.scatter(Q_fpos['z_clusterspec'],Q_fpos['z_clusterphot'], c='r', marker='x', linewidths = 0)
Passz_fneg = delzPass.scatter(Q_fneg['z_clusterspec'],Q_fneg['z_clusterphot'], c='m', marker='x', linewidths = 0)
plt.xlabel('(z_spec - z_cluster)/(1+z_spec)')
plt.xscale('linear')
plt.xlim(-0.125,0.125)
plt.tick_params(axis='both', which='both',direction='in',color='k')
plt.ylabel('(z_phot - z_cluster)/(1+z_phot)')
delzPass.yaxis.tick_right()
delzPass.yaxis.set_label_position("right")
plt.yscale('linear')
plt.ylim(-0.15,0.15)
lin = delzPass.plot([-0.5,1],[0.05,0.05],':k', linewidth=1)  # horizontal cuts
lin = delzPass.plot([-0.5,1],[-0.05,-0.05],':k', linewidth=1)
lin = delzPass.plot([-0.01,-0.01],[-0.5,1],':k', linewidth=1)  #vertical cuts
lin = delzPass.plot([0.01,0.01],[-0.5,1],':k', linewidth=1)  
#plt.legend((SFz,Passz),('del_z < 0.1','del_z > 0.1'),scatterpoints=1,loc='lower left', fontsize=10)
plt.title('Passive')
#plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
plt.tick_params(axis='y', which='both', direction='in',color='k')

plt.show()

#colour-colour plots; modify to add in photometric sub-sample
colcol = plt.figure(num=3)
cvc = colcol.add_subplot(111)
SF_vj = Column(np.concatenate((SF_spec['vj'],SF_field['vj'],SF_fpos['vj'],SF_fneg['vj'])))
SF_uv = Column(np.concatenate((SF_spec['uv'],SF_field['uv'],SF_fpos['uv'],SF_fneg['uv'])))
Q_vj = Column(np.concatenate((Q_spec['vj'],Q_field['vj'],Q_fpos['vj'],Q_fneg['vj'])))
Q_uv = Column(np.concatenate((Q_spec['uv'],Q_field['uv'],Q_fpos['uv'],Q_fneg['uv'])))
ccSF = cvc.scatter(SF_vj,SF_uv, c = 'b', marker = '*', linewidths = 0)
ccPass = cvc.scatter(Q_vj,Q_uv, c = 'r', marker = '.', linewidths = 0)
#ccout = cvc.plot(SF['vj'],SF['uv'], c = )
bounds = cvc.plot([0,0.7],[1.3,1.3],'-k',[0.7,1.6],[1.3,2],'-k',[1.6,1.6],[2,2.5],'-k', linewidth=1) #overlay boundary cutoff for SF/Passive
plt.xscale('linear')
plt.xlabel('(V-J)rest')
plt.xlim(0,2)
plt.yscale('linear')
plt.ylabel('(U-V)rest')
plt.ylim(0,2.5)
plt.title('V-J vs U-V')
plt.text(0.25,2,'Passive',fontsize=10)
plt.text(1.5,0.75,'Star-Forming',fontsize=10)
#plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')

plt.show()