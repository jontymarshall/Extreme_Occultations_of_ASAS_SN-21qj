#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 16:04:00 2022

@author: jonty
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from astropy.io import ascii

direc = '/Users/jonty/mydata/asassn21qj/'

#ASASSN 21qj ASASSN observations
#asassn_data = ascii.read(direc+'ASASSN_21qj_ASASSN_light_curve.csv',guess=False,delimiter=',',data_start=2)
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
lco_data1 = ascii.read(direc+'asassn21qj_magnitudes.csv',guess=False,delimiter=',',data_start=2) #DDT
lco_data2 = ascii.read(direc+'../little_dippers/ASSASN-21qj_magnitudes.csv',guess=False,delimiter=',',data_start=2) #Survey

lco_mjd = np.append(lco_data1['Time'].data, lco_data2['Time'].data)
lco_mag = np.append(lco_data1['Mag'].data, lco_data2['Mag'].data)
lco_err = np.append(lco_data1['MagErr'].data, lco_data2['MagErr'].data)
lco_fil = np.append(lco_data1['Filter'].data, lco_data2['Filter'].data)

green  = np.where(lco_fil == 'gp')
red    = np.where(lco_fil == 'rp')
orange = np.where(lco_fil == 'ip')
brown  = np.where(lco_fil == 'zs')

#Make single panel LCO plot

fig,ax = plt.subplots(1,1)
ax.set_xlabel(r"Modified Julian Date")
ax.set_ylabel(r"Magnitude")
ax.invert_yaxis()


mjdmin = 59450
mjdmax = 60000

ax.set_xlim((mjdmin,mjdmax))
ax.set_ylim((21,13))

ax.errorbar(asassn_mjd,asassn_mag,yerr=asassn_err,marker='s',markersize=4,linestyle='',mec='green',mfc='lightgreen',ecolor='green',elinewidth=0.5,alpha=1)
ax.errorbar(lco_mjd[green],lco_mag[green],yerr=lco_err[green],marker='o',markersize=4,linestyle='',mec='green',mfc='lightgreen',ecolor='green',elinewidth=0.5,alpha=1,label=r'$g^{\prime}$')
ax.errorbar(lco_mjd[red],lco_mag[red],yerr=lco_err[red],marker='o',linestyle='',markersize=4,mec='darkred',mfc='red',ecolor='darkred',elinewidth=0.5,alpha=1,label=r'$r^{\prime}$')
ax.errorbar(lco_mjd[orange],lco_mag[orange],yerr=lco_err[orange],marker='o',markersize=4,linestyle='',mec='blue',mfc='dodgerblue',ecolor='blue',elinewidth=0.5,alpha=1,label=r'$i^{\prime}$')
ax.errorbar(lco_mjd[brown],lco_mag[brown],yerr=lco_err[brown],marker='o',markersize=4,linestyle='',mec='darkorange',mfc='orange',ecolor='orange',elinewidth=0.5,alpha=1,label=r'$z^{\prime}$')
ax.plot([59400,59800],[13.8,13.8],linestyle='--',marker='',color='green',alpha=0.25)
ax.plot([59400,59800],[13.37,13.37],linestyle='--',marker='',color='red',alpha=0.25)
ax.plot([59400,59800],[13.25,13.25],linestyle='--',marker='',color='blue',alpha=0.25)
ax.plot([59400,59800],[13.22,13.22],linestyle='--',marker='',color='orange',alpha=0.25)
ax.plot([59476.1,59476.1],[21,13],linestyle='--',marker='',color='grey',alpha=0.25)
ax.plot([59507,59507],[21,13],linestyle='--',marker='',color='grey',alpha=0.25)
ax.plot([59537.9,59537.9],[21,13],linestyle='--',marker='',color='grey',alpha=0.25)
ax.plot([59568.8,59568.8],[21,13],linestyle='--',marker='',color='grey',alpha=0.25)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(16)

ax.legend()


plt.tight_layout()
fig.savefig(direc+'asassn21qj_timeseries_closeup_R1.pdf',dpi=200)
plt.close()

#ASASSN 21qj NEOWISE observations

neo_data2 = ascii.read(direc+'asassn21qj_neowise_data_complete.csv',guess=False,delimiter=',',data_start=2) #NEOWISE all magnitudes

neo_all_mjd = neo_data2["mjd"].data

neo_all_w1_mag = neo_data2["w1mpro"].data
neo_all_w1_err = neo_data2["w1sigmpro"].data

neo_all_w2_mag = neo_data2["w2mpro"].data
neo_all_w2_err = neo_data2["w2sigmpro"].data

latest = np.where((neo_all_mjd > 59600)&(neo_all_mjd < 59750))

neo_w1_mag_latest = np.mean(neo_all_w1_mag[latest]) 
neo_w1_err_latest = np.mean(neo_all_w1_err[latest])

neo_w2_mag_latest = np.mean(neo_all_w2_mag[latest])
neo_w2_err_latest = np.mean(neo_all_w2_err[latest])


neo_data1 = ascii.read(direc+'neowise_avg_asassn21qj_magnitudes.csv',guess=False,delimiter=',',data_start=2) #NEOWISE avg magnitudes

neo_avg_mjd = neo_data1["MJD"].data

neo_avg_w1_mag = neo_data1["W1Mag"].data
neo_avg_w1_err = neo_data1["W1Err"].data

neo_avg_w2_mag = neo_data1["W2Mag"].data
neo_avg_w2_err = neo_data1["W2Err"].data

latest = np.where((neo_all_mjd > 59600)&(neo_all_mjd < 59750))

neo_w1_mag_latest = np.mean(neo_all_w1_mag[latest]) 
neo_w1_err_latest = np.mean(neo_all_w1_err[latest])

neo_w2_mag_latest = np.mean(neo_all_w2_mag[latest])
neo_w2_err_latest = np.mean(neo_all_w2_err[latest])

neo_avg_mjd = np.append(neo_avg_mjd,np.mean(neo_all_mjd[latest]))

neo_avg_w1_mag = np.append(neo_avg_w1_mag,neo_w1_mag_latest)
neo_avg_w1_err = np.append(neo_avg_w1_err,neo_w1_err_latest)

neo_avg_w2_mag = np.append(neo_avg_w2_mag,neo_w2_mag_latest)
neo_avg_w2_err = np.append(neo_avg_w2_err,neo_w2_err_latest)


latest = np.where((neo_all_mjd > 59750))

neo_w1_mag_latest = np.mean(neo_all_w1_mag[latest]) 
neo_w1_err_latest = np.mean(neo_all_w1_err[latest])

neo_w2_mag_latest = np.mean(neo_all_w2_mag[latest])
neo_w2_err_latest = np.mean(neo_all_w2_err[latest])

neo_avg_mjd = np.append(neo_avg_mjd,np.mean(neo_all_mjd[latest]))

neo_avg_w1_mag = np.append(neo_avg_w1_mag,neo_w1_mag_latest)
neo_avg_w1_err = np.append(neo_avg_w1_err,neo_w1_err_latest)

neo_avg_w2_mag = np.append(neo_avg_w2_mag,neo_w2_mag_latest)
neo_avg_w2_err = np.append(neo_avg_w2_err,neo_w2_err_latest)

#Make 3 panel plot

fig,(ax1,ax2,ax3,ax4) = plt.subplots(4,1,sharex=True, facecolor='w',gridspec_kw={'height_ratios': [2.5, 2, 1, 1]})
fig.set_size_inches((8,8))
fig.set_dpi(200)
ax4.set_xlabel(r"Modified Julian Date")

#ax1.set_ylabel(r"Magnitude")
ax2.set_ylabel(r"Magnitude")
#ax3.set_ylabel(r"Magnitude")
ax1.invert_yaxis()
ax2.invert_yaxis()
ax3.invert_yaxis()
ax4.invert_yaxis()

mjdmin = 57250
mjdmax = 60000

ax1.set_xlim((mjdmin,mjdmax))
ax1.set_ylim((18,13))

ax1.errorbar(asassn_mjd,asassn_mag,yerr=asassn_err,marker='s',linestyle='',markersize=2,mec='green',mfc='lightgreen',ecolor='green',alpha=0.5,label=r'$V$, $g^{\prime}$')
ax1.errorbar(lco_mjd[green],lco_mag[green],yerr=lco_err[green],marker='o',linestyle='',markersize=2,mec='lightgreen',mfc='lightgreen',ecolor='lightgreen',alpha=0.5)
#ax1.errorbar(lco_mjd[red],lco_mag[red],yerr=lco_err[red],marker='o',linestyle='',mec='red',mfc='red',ecolor='darkred',alpha=0.1)
#ax1.errorbar(lco_mjd[orange],lco_mag[orange],yerr=lco_err[orange],marker='o',linestyle='',mec='orange',mfc='orange',ecolor='orange',alpha=0.1)
#ax1.errorbar(lco_mjd[brown],lco_mag[brown],yerr=lco_err[brown],marker='o',linestyle='',mec='brown',mfc='brown',ecolor='brown',alpha=0.1)
ax1.plot([np.min(neo_all_mjd),np.max(lco_mjd)],[13.8,13.8],linestyle='--',color='green')
lb = 3*np.mean(asassn_err)
box_x = [57250,57250,59750,59750,57250]
box_y = [13.8-lb,13.8+lb,13.8+lb,13.8-lb,13.8-lb]
ax1.fill(box_x,box_y, fill=True, color="green",alpha=0.1)

ax1.add_patch(plt.Rectangle((.3, .8), .12, .1, ls="-", ec="k", fc="none",
                           transform=ax1.transAxes))
ax1.legend(loc='upper left')

axins = inset_axes(ax1, bbox_to_anchor=(.1, .1, .3, .4), bbox_transform=ax1.transAxes,width=4, height=1, loc=3)
axins.set_xlim((58000,58300))
axins.set_ylim((14,13.5))
axins.plot([np.min(neo_all_mjd),np.max(lco_mjd)],[13.8,13.8],linestyle='--',color='green')
axins.errorbar(asassn_mjd,asassn_mag,yerr=asassn_err,marker='s',linestyle='',markersize=2,mec='green',mfc='lightgreen',ecolor='green',alpha=0.5)
box_x = [mjdmin,mjdmin,mjdmax,mjdmax,mjdmin]
box_y = [13.8-lb,13.8+lb,13.8+lb,13.8-lb,13.8-lb]
axins.fill(box_x,box_y, fill=True, color="green",alpha=0.1)

atlas_data = ascii.read(direc+'atlas_forced_photometry_asassn21qj.txt')

atlas_tim = atlas_data["MJD"]
atlas_mag = atlas_data["m"]
atlas_err = atlas_data["dm"]
atlas_fil = atlas_data["F"]

good = np.where((atlas_fil == 'o')&(atlas_mag > 13.)&(atlas_mag < 20.)&(atlas_err < 0.1))

atlas_mjd = atlas_tim[good]
atlas_c_mag = atlas_mag[good]
atlas_c_err = atlas_err[good]

ax2.set_ylim(15,13)
lb = 3*np.mean(atlas_c_err)
noslope = np.where((atlas_mjd <= 58500))
cmag = np.mean(atlas_c_mag[noslope])
ax2.errorbar(atlas_mjd,atlas_c_mag,yerr=atlas_c_err,marker='o',linestyle='',markersize=2,mec='blue',mfc='dodgerblue',ecolor='blue',alpha=0.5,label=r'ATLAS cyan')
ax2.plot([57000,60000],[cmag,cmag],linestyle='--',color='blue')
box_x = [mjdmin,mjdmin,mjdmax,mjdmax,mjdmin]
box_y = [cmag-lb,cmag+lb,cmag+lb,cmag-lb,cmag-lb]
ax2.fill(box_x,box_y, fill=True, color="dodgerblue",alpha=0.1)
ax2.legend(loc='upper left')

ax3.set_ylim((12,10.5))

ax3.errorbar(neo_all_mjd,neo_all_w1_mag,yerr=neo_all_w1_err,linestyle='',markersize=4,marker='o',mec='orange',mfc='white',ecolor='orange',alpha=0.25)
ax3.errorbar(neo_avg_mjd,neo_avg_w1_mag,yerr=neo_avg_w1_err,linestyle='',markersize=4,marker='o',mec='orange',mfc='yellow',ecolor='orange',label=r'NEOWISE W1')
ax3.plot([np.min(neo_all_mjd),np.max(lco_mjd)],[11.681,11.681],linestyle='--',color='orange')

box_x = [mjdmin,mjdmin,mjdmax,mjdmax,mjdmin]
box_y = [11.531,11.831,11.831,11.531,11.531]
ax3.fill(box_x,box_y, fill=True, color="orange",alpha=0.1)
ax3.legend(loc='upper left')

ax4.set_ylim((12,10.5))

ax4.errorbar(neo_all_mjd,neo_all_w2_mag,yerr=neo_all_w2_err,linestyle='',markersize=4,marker='o',mec='red',mfc='white',ecolor='red',alpha=0.25)
ax4.errorbar(neo_avg_mjd,neo_avg_w2_mag,yerr=neo_avg_w2_err,linestyle='',markersize=4,marker='o',mec='darkred',mfc='red',ecolor='red',label=r'NEOWISE W2')
ax4.plot([np.min(neo_all_mjd),np.max(lco_mjd)],[11.749,11.749],linestyle='--',color='red')

box_x = [mjdmin,mjdmin,mjdmax,mjdmax,mjdmin]
box_y = [11.6,11.9,11.9,11.6,11.6]
ax4.fill(box_x,box_y, fill=True, color="red",alpha=0.1)
ax4.legend(loc='upper left')

plt.tight_layout()
fig.savefig(direc+'asassn21qj_multipanel_plot_R1b.pdf',dpi=300)
plt.close()

# Colour-colour plot
lco_colours = ascii.read(direc+'asassn21qj_colours.csv',guess=False,delimiter=',',data_start=2)

lco_gr = lco_colours["g-r"].data
lco_ri = lco_colours["r-i"].data
lco_ugr = lco_colours["Err g-r"].data
lco_uri = lco_colours["Err r-i"].data

fig,ax = plt.subplots(1,1)
ax.set_xlabel(r"$g^{\prime}$ - $r^{\prime}$ (mag)")
ax.set_ylabel(r"$r^{\prime}$ - $i^{\prime}$ (mag)")

ax.set_xlim((0,2))
ax.set_ylim((0,1))

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(16)

ax.errorbar(lco_gr,lco_ri,xerr=np.sqrt(lco_ugr**2 + 0.02**2),yerr=np.sqrt(lco_uri**2 + 0.02**2),marker='o',markersize=6,linestyle='',mec='darkgreen',mfc='green',ecolor='darkgreen',elinewidth=0.5,alpha=0.5)
ax.plot([0.43],[0.12],marker='*',markersize=36,linestyle='',mec='orange',mfc='yellow')
plt.tight_layout()
fig.savefig(direc+'asassn21qj_colour_plot_R1b.pdf',dpi=200)
plt.close()

