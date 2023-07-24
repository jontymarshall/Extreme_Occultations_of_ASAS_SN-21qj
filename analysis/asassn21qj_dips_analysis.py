#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 16:04:00 2022

@author: jonty
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy.interpolate import interp1d

direc = '/Users/jonty/mydata/asassn21qj/'

#ASASSN 21qj ASASSN observations
asassn_data = ascii.read(direc+'ASASSN_lightcurve_ful_aperturephot.csv',guess=False,delimiter=',',data_start=2)


asassn_hjd = np.asarray(asassn_data['HJD'].data)
asassn_mag = np.asarray(asassn_data['mag'].data)
asassn_err = np.asarray(asassn_data['mag_err'].data)
asassn_fil = np.asarray(asassn_data['Filter'].data)


#for i in range(0,len(asassn_mag)):
#    print(i,asassn_mag[i],asassn_err[i])

good = np.where((abs(asassn_mag / asassn_err) > 3.) & (asassn_mag > 12) )


asassn_mjd = asassn_hjd - 2400000.0 #HJD -> MJD

asassn_mjd = asassn_mjd[good]
asassn_mag = asassn_mag[good]
asassn_err = asassn_err[good]
asassn_fil = asassn_fil[good]

asassn_v = np.where(asassn_fil == 'V')

asassn_mag[asassn_v] = asassn_mag[asassn_v] + 0.65*0.65 - 0.12

#ASASSN 21qj LCO observations
lco_data = ascii.read(direc+'asassn21qj_lco_magnitudes_gp_combined.csv',guess=False,delimiter=',',data_start=2) #DDT

lco_mjd = lco_data['Time'].data
lco_mag = lco_data['Mag'].data
lco_err = lco_data['MagErr'].data
lco_fil = lco_data['Filter'].data

t_all = np.append(asassn_mjd,lco_mjd)
g_all = np.append(asassn_mag,lco_mag)
u_all = np.append(asassn_err,lco_err)


#Identify significant dips in the lightcurve
toc = []
depth = []
for i in range(1,len(g_all)-1):
    
    dm_prev = (g_all[i-1] - g_all[i]) / np.sqrt(u_all[i-1]**2 + u_all[i]**2)
    dm_next = (g_all[i+1] - g_all[i]) / np.sqrt(u_all[i+1]**2 + u_all[i]**2)
#    print(t_all[i],dm_prev,dm_next)
    if dm_prev < -3 and dm_next < -3 :
        toc.append(t_all[i])
        depth.append(1. - 10**(-0.4*(g_all[i] - 13.8)))
#        print(tg[i],mg[i-1],mg[i],mg[i+1])

#Plot up results
#fig, ax1 = plt.subplots()
#ax1.set_xlabel(r"time (MJD)")
#ax1.set_ylabel(r"$g^{\prime}$ (mag)")
#ax1.invert_yaxis()
#ax1.errorbar(t_all,g_all,xerr=None,yerr=u_all,linestyle='',marker='o',mec='white',mfc='green',ecolor='green',alpha=0.5,label="sdss gp")
#for i in range(0,len(toc)):
#    ax1.plot([toc[i],toc[i]],[19.0,19.0],marker='^',linestyle='',color='black')
#ax1.plot([np.min(g_all),np.max(g_all)],[13.8,13.8],linestyle=':',marker='',color='black',alpha=1.0)
#ax1.set_xlim(59300,59900)
#ax1.set_ylim(21.,13.)
#ax1.legend(loc='lower right')
#plt.tight_layout()
#fig.savefig(direc + 'asassn21qj_occultations_plot.pdf',dpi=200)
#plt.show()
#plt.close()
#plt.clf()

#n, bins, patches = plt.hist(np.log10(depth), 10, density=False, facecolor='g', alpha=0.75)
#plt.xlabel('Depth')
#plt.ylabel('Number')
#plt.xlim(0., -1.)
#plt.ylim(0, 10.)
#plt.show()

#Lomb Scargle periodogram of photometry
from astropy.timeseries import LombScargle
freq_asn, power_asn = LombScargle(asassn_mjd, asassn_mag, asassn_err).autopower(minimum_frequency=(1/130),maximum_frequency=(1/10))
freq_lco, power_lco = LombScargle(lco_mjd, lco_mag, lco_err).autopower(minimum_frequency=(1/130),maximum_frequency=(1/10))
freq_all, power_all = LombScargle(t_all, g_all, u_all).autopower(minimum_frequency=(1/130),maximum_frequency=(1/10))

print(power_asn.max(),1/freq_asn[np.where(power_asn == power_asn.max())])
print(power_lco.max(),1/freq_lco[np.where(power_lco == power_lco.max())])
print(power_all.max(),1/freq_all[np.where(power_all == power_all.max())])

probabilities = [0.1, 0.05, 0.01]

print(LombScargle(asassn_mjd, asassn_mag, asassn_err).false_alarm_level(probabilities))
print(LombScargle(asassn_mjd, asassn_mag, asassn_err).false_alarm_probability(np.max(power_asn)))

print(LombScargle(lco_mjd, lco_mag, lco_err).false_alarm_level(probabilities))
print(LombScargle(lco_mjd, lco_mag, lco_err).false_alarm_probability(np.max(power_lco)))

print(LombScargle(t_all, g_all, u_all).false_alarm_level(probabilities))
print(LombScargle(t_all, g_all, u_all).false_alarm_probability(np.max(power_all)))

plt.plot(1/freq_lco, power_lco,color='blue',label='LCO only')
plt.plot(1/freq_all,power_all,color='green',label='Both')
plt.plot(1/freq_asn, power_asn,color='orange',label='ASASSN only')
plt.legend()
plt.xlabel(r"Period (days)",fontsize=18)
plt.ylabel(r"Magnitude",fontsize=18)

#Model overall extinction with asymmetric Gaussian
stdv_long  = 100
stdv_short = 50
depth = 1.0
minimum = 59450
t_cloud = np.linspace(59000,60000,4000)
m_cloud = 13.8 + depth*np.e**(-0.5*((minimum-t_cloud)/stdv_short)**2)
tail = np.where(t_cloud >= minimum)
m_cloud[tail] = 13.8 + depth*np.e**(-0.5*((minimum-t_cloud[tail])/stdv_long)**2)

#Model narrow feature at 59476.5 days
stdv = 0.6
depth = 6.0
minimum = 59476.1
t_dip = np.linspace(59000,60000,4000)
m_dip = m_cloud + depth*np.e**(-0.5*((minimum-t_dip)/stdv)**2)

times = np.sort(toc)
tdips = t_cloud
mdips = m_cloud
depths = [0,6.0,2.0,9.5,1,
		  0.3,1.1,0.9,2.0,1.8,
		  1.2,1.3,0.2,0.4,0.3,
		  1,0.8,1.4,1.2,1,
		  1.3,0.8,2.2,0.5,0.9,
		  0.7,0.6,0.9,1.75,0.75,
		  0.4,0.2,0,0,0]
widths = [1,0.6,0.6,0.1,2.0,
		  1,1,1,1.0,0.8,
		  1,1,1,1,1,
		  1,1,1,1,1,
		  1,1,0.6,1,0.6,
		  1,1,1,0.5,1.2,
		  1.8,1,1,1,1]

for i in range(0,len(times)):
	#print(toc[i],depths[i],widths[i])
	mdips += depths[i]*np.e**(-0.5*((times[i]-t_cloud)/widths[i])**2)

print(np.min(mdips),np.max(mdips))

f = interp1d(tdips,mdips)

interp = np.where((t_all >=59000)&(t_all < 60000))
residuals = g_all[interp] - f(t_all[interp])

print(np.min(residuals),np.max(residuals))

fig,ax = plt.subplots(1,1)
fig.set_size_inches((12,6))
ax.set_xlabel(r"Modified Julian Date",fontsize=18)
ax.set_ylabel(r"Magnitude",fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
ax.set_xlim((59300,59900))
ax.set_ylim((21,13.5))
ax.errorbar(t_all,g_all,xerr=None,yerr=u_all,marker='.',linestyle='',mec='white',mfc='green',ecolor='green',alpha=0.5)
ax.plot([59000,60000],[13.8,13.8],marker=None,linestyle='--',color='green')
ax.plot(t_cloud,m_cloud,marker=None,linestyle=':',color='black')
ax.plot(tdips,mdips,marker=None,linestyle='-',color='black')
for i in range(0,len(toc)):
    ax.plot([toc[i],toc[i]],[19.0,19.0],marker='^',linestyle='',color='forestgreen',alpha=0.25)
plt.tight_layout()
fig.savefig(direc+'asassn21qj_dimming.pdf',dpi=200)
plt.show()
plt.clf()