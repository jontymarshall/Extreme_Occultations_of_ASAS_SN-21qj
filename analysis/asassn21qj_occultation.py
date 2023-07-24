#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 13:21:33 2022

@author: jonty
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii

direc = '/Users/jonty/mydata/asassn21qj/'
target = 'asassn21qj'
filename = 'asassn21qj_magnitudes.csv'
table = ascii.read(filename,comment='#')

tob = table['Time'].data
mag = table['Mag'].data
err = table['MagErr'].data
filt = table['Filter'].data

mtgt = 13.7
gband = np.where(filt =='gp')

tg = tob[gband]
mg = mag[gband]
ug = err[gband]

toc = []
depth = []
for i in range(1,len(mg)-1):
    
    dm_prev = (mg[i-1] - mg[i]) / np.sqrt(ug[i-1]**2 + ug[i]**2)
    dm_next = (mg[i+1] - mg[i]) / np.sqrt(ug[i+1]**2 + ug[i]**2)
    print(tg[i],dm_prev,dm_next)
    if dm_prev < -3 and dm_next < -3 :
        toc.append(tg[i])
        depth.append(1. - 10**(-0.4*(mg[i] - 13.7)))
#        print(tg[i],mg[i-1],mg[i],mg[i+1])

#Plot up results
fig, ax1 = plt.subplots()
ax1.set_xlabel(r"time (MJD)")
ax1.set_ylabel(r"$g^{\prime}$ (mag)")
ax1.invert_yaxis()
ax1.errorbar(tg,mg,xerr=None,yerr=ug,linestyle='',marker='o',mfc='gray',mec='darkgray',ecolor='gray',color='gray',alpha=0.5,label="sdss gp")
for i in range(0,len(toc)):
    ax1.plot([toc[i],toc[i]],[19.0,20.0],marker='',linestyle='-',color='black')
ax1.plot([np.min(tg),np.max(tg)],[13.7,13.7],linestyle=':',marker='',color='black',alpha=1.0)
ax1.set_xlim(np.min(tg) - 10,np.max(tg) + 10)
ax1.set_ylim(21.,13.)
ax1.legend(loc='lower right')
fig.savefig(direc + target + '_lcogt_occultations_plot.png',dpi=200)
plt.show()
plt.close()
plt.clf()

n, bins, patches = plt.hist(np.log10(depth), 10, density=True, facecolor='b', alpha=0.75)
plt.xlabel('Depth')
plt.ylabel('Number')
plt.xlim(0., -1.)
plt.ylim(0, 10.)
plt.show()