#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 09:17:20 2017

@author: gsarrouh
"""
###This program creates plots for the master data catalogue, containing:
###photometric v specroscopic redshift, delta-z, & colour-colour plots
plt.close()
plt.clf()
#
## PLOT 1: make z_spec vs z_phot plot
#
zees = plt.figure(num=1)
zvz = zees.add_subplot(111)
##plot spec subsample using master_cat['sub']==1 as identifier for loop
#plotzs = np.array
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['sub'] != 1 or master_cat[counter]['type'] ==0:      #skip no data and phot subsample & stars
        counter +=1
    elif abs(master_cat[counter]['del_z']) > 0.15:
        outz = zvz.scatter(master_cat[counter]['z_spec'],master_cat[counter]['z_peak'],c='r', marker='v', linewidths = 0)
        counter +=1
    else:
        memz = zvz.scatter(master_cat[counter]['z_spec'],master_cat[counter]['z_peak'],c='b', marker='^', linewidths = 0)
        counter +=1
plt.plot([0,2],[0,2],':k', linewidth=1)
plt.xlabel("$z^{spec}$")
plt.xscale('linear')
plt.xlim(0,1.25)
plt.ylabel('z$phot$')
plt.yscale('linear')
plt.ylim(0,1.25)
plt.title('z_spec vs. z_photo')
plt.legend((memz,outz),('spec population','outliers'),scatterpoints=1,loc='upper left', frameon=False)
#plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = ':')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on')
plt.minorticks_on()
plt.show()
#
#
## PLOT 2: del_z(spec) vs del_z(phot)
#
delzees = plt.figure(num=2)
delzSF = delzees.add_subplot(121) #SF plot
#
SF_membership = ([0]*4) # [(memebrs),(field),(false pos),(false neg)]
Q_membership = ([0]*4) # [(memebrs),(field),(false pos),(false neg)]
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['sub'] == 1 and master_cat[counter]['type']==1:      #only look at SF spec sub-sample, type 1 = SF
        if master_cat[counter]['member'] == 0:        #secure cluster member
            SF_membership[0] +=1
            Smem = delzSF.scatter(master_cat[counter]['z_clusterspec'],master_cat[counter]['z_clusterphot'],c='g', marker='+', linewidths = 0)
        elif master_cat[counter]['member'] == 1:        #secure field
            SF_membership[1] +=1
            Sfield = delzSF.scatter(master_cat[counter]['z_clusterspec'],master_cat[counter]['z_clusterphot'],c='b', marker='+', linewidths = 0)
        elif master_cat[counter]['member'] == 2:        #false pos
            SF_membership[2] +=1
            Spos = delzSF.scatter(master_cat[counter]['z_clusterspec'],master_cat[counter]['z_clusterphot'],c='r', marker='x', linewidths = 0)
        elif master_cat[counter]['member'] == 3:        #false neg
            SF_membership[3] +=1
            Sneg = delzSF.scatter(master_cat[counter]['z_clusterspec'],master_cat[counter]['z_clusterphot'],c='m', marker='x', linewidths = 0)
    counter +=1
#
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
plt.legend((Smem,Sfield,Spos,Sneg),('secure member','secure field','false positive','false negative'),scatterpoints=1,loc='upper left', fontsize=8, frameon=False)
plt.title('Star-Forming')
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off',labelleft='on')

#
delzPass = delzees.add_subplot(122) #Passive plot
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['sub'] == 1 and master_cat[counter]['type']==2:      #only look at Q spec sub-sample, type 2 = Q
        if master_cat[counter]['member'] == 0:        #secure cluster member
            Q_membership[0] +=1
            Qmem = delzPass.scatter(master_cat[counter]['z_clusterspec'],master_cat[counter]['z_clusterphot'],c='g', marker='+', linewidths = 0)
        elif master_cat[counter]['member'] == 1:        #secure field
            Q_membership[1] +=1
            Qfield = delzPass.scatter(master_cat[counter]['z_clusterspec'],master_cat[counter]['z_clusterphot'],c='b', marker='+', linewidths = 0)
        elif master_cat[counter]['member'] == 2:        #false pos
            Q_membership[2] +=1
            Qpos = delzPass.scatter(master_cat[counter]['z_clusterspec'],master_cat[counter]['z_clusterphot'],c='r', marker='x', linewidths = 0)
        elif master_cat[counter]['member'] == 3:        #false neg
            Q_membership[3] +=1
            Qneg = delzPass.scatter(master_cat[counter]['z_clusterspec'],master_cat[counter]['z_clusterphot'],c='m', marker='x', linewidths = 0)
    counter +=1
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
plt.minorticks_on()
plt.tick_params(axis='y', which='both', direction='in',color='k',top='on',right='on',labelright='on',labelleft='off')

plt.show()

#colour-colour plots; modify to add in photometric sub-sample
check = 0
aa = 0
bb = 0          #counting variables to check plotting code
cc = 0
colcol = plt.figure(num=3)
cvc = colcol.add_subplot(111)
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['member'] ==0:
        if master_cat[counter]['type'] ==1:
            aa +=1      #count for SF
            SF = cvc.scatter(master_cat[counter]['vj'],master_cat[counter]['uv'], c='b', marker='*', linewidths=0)
        elif master_cat[counter]['type']==2:
            bb +=1    #count for Q
            Q = cvc.scatter(master_cat[counter]['vj'],master_cat[counter]['uv'], c='r', marker='.', linewidths=0)
        else:
            cc +=1
    else:
        check +=1
    counter +=1
#    if master_cat[counter]['sub'] ==1 or master_cat[counter]['sub'] == 2:   #look at spec & phot only, not 'nodata'
#        if master_cat[counter]['type']==1:
#            aa +=1      #count for SF
#            SF = cvc.scatter(master_cat[counter]['vj'],master_cat[counter]['uv'], c='b', marker='*', linewidths=0)
#        elif master_cat[counter]['type']==2:
#            bb +=1    #count for Q
#            Q = cvc.scatter(master_cat[counter]['vj'],master_cat[counter]['uv'], c='r', marker='.', linewidths=0)
#        else:
#            dd +=1 #count for stars & outliers
#    else:
#        check +=1  #count for nodata
#    counter +=1   
bounds = cvc.plot([0,0.8],[1.3,1.3],'-k',[0.8,1.6],[1.3,2],'-k',[1.6,1.6],[2,2.5],'-k', linewidth=1) #overlay boundary cutoff for SF/Passive
plt.xscale('linear')
plt.xlabel('(V-J)rest')
plt.xlim(0,2)
plt.yscale('linear')
plt.ylabel('(U-V)rest')
plt.ylim(0,2.5)
plt.title('V-J vs U-V')
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on')
plt.minorticks_on()
plt.text(0.25,2,'Passive',fontsize=10)
plt.text(1.5,0.75,'Star-Forming',fontsize=10)
#plt.grid(b=True, which='major', axis='both', color = 'k', linestyle = '--')

plt.show()