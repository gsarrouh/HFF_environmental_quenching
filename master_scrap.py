#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 17:24:27 2017

@author: gsarrouh
"""
aa=[]
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['type'] == 1 or master_cat[counter]['type'] == 2:  #SF & Q     
        if master_cat[counter]['sub'] == 1:                 #objects w/ both spec & phot only
            aa.append(master_cat[counter]['del_z'])
    counter +=1
##
a = 0
b = 0
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['member'] == 3 and master_cat[counter]['type'] ==1:
        if master_cat[counter]['lmass'] > 7.4 :
            a+=1
        else: b+=1
    counter +=1
#
###########################################        
err_s = f_s/serr    #calculate flux / delta_fluc = y = 10^(-Ax + B)

A,B = np.polyfit(smag, np.log(err_s), 1)   # calculate polynomial
f = np.poly1d([A,B])
smag_new = np.linspace(min(smag), max(smag), 200)  # calculate new x's and y's
err_s_new = f_spec(smag_new)

def func(x,a,b,c):
    return a * np.power(10,b*x)+c

sfit, s_cov = curve_fit(func, smag, err_s,  p0=(100, 0.,500))
pfit, p_cov = curve_fit(func, pmag, err_p,  p0=(50, 0.,1000))

#sfit, s_cov = curve_fit(lambda x,a,b,c: np.log10(a)+c,smag,err_s, p0=(1, 0.,10000))
#pfit, p_cov = curve_fit(lambda x,a,b,c: a*np.log10(-b*x)+c,pmag,err_p,  p0=(1, 0.,10000))
s_yfit = [func(x,sfit[0],sfit[1],sfit[2]) for x in smag]
p_yfit = [func(x,pfit[0],pfit[1],pfit[2]) for x in pmag]#
########################################

X = [1,2,3,4,5,6]
C = [[0]*6]
counter = 0
size = len(master_cat)
while counter < size:
    for x in X:
        if master_cat[counter]['cluster'] == x:# for x in X:
            C[(x-1)] +=1
        counter+=1


plt.close()
e_F160 = plt.figure(num=2)
plt.suptitle('Flux/Error vs. Apparent Magnitude')
spec_error = e_F160.add_subplot(121)
e_spec = spec_error.scatter(smag,err_s,c='g',marker='.',linewidth=0.5)
fit_spec = plt.scatter(smag, s_yfit,c='k',marker='.',linewidth=0.1)#,label='Fitted Curve')
plt.plot([0,40],[5,5],':r', linewidth=1)
#plt.plot(smag_new,err_s_new,'-k',linewidth=1)
plt.xlabel("$m_{F160W}$")
plt.xscale('linear')
#plt.gca().invert_xaxis()
plt.xlim(15,40)
plt.ylabel("$flux_{F160W}/error$")
plt.yscale('log')
#plt.ylim(5,15)
plt.title('Spec')
#plt.legend((e_spec,e_phot),('spec','phot'),scatterpoints=1,loc='best', frameon=False)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='off')
#
phot_error = e_F160.add_subplot(122)
e_phot = phot_error.scatter(pmag,err_p,c='m',marker='.',linewidth=0.5)
fit_phot = phot_error.scatter(pmag, p_yfit,c='k',marker='.',linewidth=0.1)#,label='Fitted Curve')
plt.plot([0,40],[5,5],':r', linewidth=0.8)
#plt.plot(pmag_new,err_p_new,'-k',linewidth=1)
plt.xlabel("$m_{F160W}$")
plt.xscale('linear')
#plt.gca().invert_xaxis()
plt.xlim(15,40)
plt.ylabel("$flux_{F160W}/error$")
phot_error.yaxis.set_label_position("right")
plt.yscale('log')
#plt.ylim(5,15)
plt.title('Phot')
#plt.legend((e_spec,e_phot),('spec','phot'),scatterpoints=1,loc='best', frameon=False)
plt.grid(b=True, which='both', axis='both', color='k', linestyle = ':',linewidth=0.15)
plt.minorticks_on()
plt.tick_params(axis='both', which='both',direction='in',color='k',top='on',right='on',labelright='on', labelleft='off')
#
plt.show()  

x = np.linspace(0,1,100)
y = np.exp(x)

plt.close()
plt.plot(x,y,'--k',linewidth=0.5)


###########################################
aa=0
bb=0
cc=0
dd=0
ee=0
ff=0
gg=0
hh=0
ss=0
other = 0
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['type'] == 0 or master_cat[counter]['type'] == 3 or master_cat[counter]['sub'] == 0:  
        ss+=1
    elif master_cat[counter]['f_F160W'] < 0 and master_cat[counter]['sub'] == 3:
        cc+=1
    elif master_cat[counter]['f_F160W'] < 0:# 
        aa+=1
    elif master_cat[counter]['sub'] == 3:
        bb+=1
 #       elif master_cat[counter]['z_peak'] <0.3:
 #           dd+=1
 #       else: ee+=1
    else: ff+=1
    counter +=1
#
SF_midbins1 = np.empty_like(SF_pos_hist1, dtype='float64')
#
size = len(binsSF1)-1
counter = 0
while counter < size:
    SF_midbins1[counter] = (binsSF1[counter] + binsSF1[counter+1])/2
    counter +=1
#
SF1 = Table([SF_midbins1,SF_pos_hist1,SF_neg_hist1],names=('SF Bins','SF pos','SF neg'))
#

A = [0,9,8,4,2]
i = 0
for x in A:
    if i < 2:
        A[i] = 100
        i +=1
    elif i >= 2 and i <4:
        A[i] = 50
        i +=1
    elif i >= 4:
        A[i] = 25
        i +=1




#        
SF_total = 0
SF_small = 0
Q_total = 0
Q_small = 0
counter = 0
size = len(master_cat)
while counter < size:
    if master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==1:  
        SF_total +=1
        if master_cat[counter]['lmass'] < 7.3:
            SF_small +=1
    elif master_cat[counter]['member'] == 0 and master_cat[counter]['type'] ==2:  
        Q_total +=1
        if master_cat[counter]['lmass'] < 7.3:
            Q_small +=1
    counter +=1