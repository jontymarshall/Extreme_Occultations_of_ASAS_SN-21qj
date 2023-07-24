#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 16:37:41 2022

@author: jonty
"""

import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
from scipy import interpolate
from astropy.io import ascii
import time
import os
from scipy.interpolate import interp1d
import scipy.odr as odr
import json

#constants
h = 6.626e-34
c = 299792458.0 # m/s
k = 1.38e-23
sb = 5.67e-8 #
au     = 1.495978707e11 # m 
pc     = 3.0857e16 # m
lsol   = 3.828e26 # W
rsol   = 6.96342e8 # m
MEarth = 5.97237e24 # kg
G = 6.67e-11
Rsol = 6.96342e8 # m
Msol = 1.99e30 #kg
Lsol = 3.828e26 #W
pc = 3.086e16 #m
au = 1.496e11 #m
um = 1.e-6 #m


direc = '/Users/jonty/mydata/asassn21qj/single_size/'

#Read in colours from LCO observations
filename = '../asassn21qj_colours.csv'
table = ascii.read(filename,comment='#')

mg = table['mg'].data
mr = table['mr'].data
mi = table['mi'].data
mz = table['mz'].data

gmr = table['g-r'].data
ugmr = table['Err g-r'].data

rmi = table['r-i'].data
urmi = table['Err r-i'].data

rmz = table['r-z'].data
urmz = table['Err r-z'].data

good = np.where(rmz != -99)

#Determine reddening slopes in g-r, r-i  and g-r, r-z space for material fitting.
def linfit_ri(B, x):
    return B[0]*(x + 0.44) + (B[1] + 0.11)
def linfit_rz(B, x):
    return B[0]*(x + 0.44) + (B[1] + 0.14)

linear = odr.Model(linfit_ri)
gr_ri_data = odr.RealData(gmr, rmi, sx=ugmr, sy=urmi)
gr_ri_odr = odr.ODR(gr_ri_data, linear, beta0=[1., 2.])
gr_ri_output = gr_ri_odr.run()
gr_ri_output.pprint()

x_fit = np.linspace(gmr[0], gmr[-1], 1000)
y_fit = linfit_ri(gr_ri_output.beta, x_fit)

plt.errorbar(gmr, rmi, xerr=ugmr, yerr=urmi, linestyle='None', marker='x')
plt.plot(x_fit, y_fit)

#plt.show()

linear = odr.Model(linfit_rz)
gr_rz_data = odr.RealData(gmr[good], rmz[good], sx=ugmr[good], sy=urmz[good])
gr_rz_odr = odr.ODR(gr_rz_data, linear,beta0=[1., 2.])
gr_rz_output = gr_rz_odr.run()
gr_rz_output.pprint()

x_fit = np.linspace(gmr[0], gmr[-1], 1000)
y_fit = linfit_rz(gr_rz_output.beta, x_fit)

plt.errorbar(gmr[good], rmz[good], xerr=ugmr[good], yerr=urmz[good], linestyle='None', marker='x')
plt.plot(x_fit, y_fit)

#plt.show()

#exit()

#wavelengths of observation for the LCO filterbands griz
wave = [0.4770,0.6231,0.7625,0.9134] #um

#brightness of a G2V star in mJy at 550pc
mtgt = 13.7
mags = np.asarray([mtgt,mtgt-0.44,mtgt-0.55,mtgt-0.58])
flux = 3631*10**(-0.4*np.asarray(mags))

distance = 550.*pc

#range of grain sizes for log-normal distributions
compositions = ['astrosil','pyr-mg100','pyr-mg80','pyr-mg60','pyr-mg40','pyr-c-mg96','ol-c-mg100','ol-mg50','c-z','c-gra','sio2','sic','fe-c']

compositions_neat = ['AstroSil','Pyroxene 100% Mg','Pyroxene 80% Mg','Pyroxene 60% Mg','Pyroxene 40% Mg','Pyroxene Cryst. 96% Mg','Olivine Cryst. 100% Mg','Olivine 50% Mg','Carbon','Graphite',r'SiO$_{2}$','SiC','Fe']

lines = ["--","--","--","--"]
linecycler = cycle(lines)
#colors = ["red","orange","yellow","green","blue","purple"]
#colorcycler = cycle(colors)

sizes = np.linspace(0.10,1.20,111)
width = 0.5

mclump = 2e-18*MEarth

def getRho():

    with open(direc+'dustkappa.dat') as f:
        lines = f.readlines()

    rho = None
    for line in lines:
        #print(line)
        linesplit = line.split()
        #print(linesplit)
        if len(linesplit) > 4 and linesplit[1] == 'grain':
            rho = float(linesplit[3])
            return rho

csq_limit= 20000.0
index = 0
for composition in compositions:
    composition_neat = compositions_neat[index]

    output_sizes_all = []
    output_slope_all = []
    output_sizes = []
    output_slope = []
    output_csq = []
    output_csq_all = []

    for size in sizes:

        print(size,composition)

        #Generate Qabs for range of grain sizes using optool
        os.system("optool "+composition+" 1.0 -a 0.01 10.0 "+str(size)+":"+str(width)+" -l 0.2 2.0 -fmax 0")

        #Read in dustkappa.dat
        optconst = ascii.read(direc+"dustkappa.dat",comment="#",data_start=2) 

        wavelength = optconst["col1"].data
        qabs = optconst["col2"].data
        qsca = optconst["col3"].data

        qext = (qabs + qsca)*1e-1 #cm^2/g : convert from cgs -> mks 

        rho = getRho()

        mgrain = (4./3.)*np.pi*(rho*1e-3)*(size*um)**3 #kg
        agrain = np.pi*(size*um)**2

        tau = qext*agrain*(mclump/mgrain)/(rsol**2) #number density (N/V) * distance cancels out factor of distance

        ext = np.e**(-1.*tau) #ratio I/Io = e**-tau, tau = sigma*N*l
        
        f = interp1d(wavelength,ext)
        extinction = f(wave)

        delta = -2.5*np.log(extinction) #-2.5*np.log10((flux[0])/3631) + 2.5*np.log10((flux[0]*extinction[0])/3631)

        col_gr = (mags[0] + delta[0]) - (mags[1] + delta[1])
        col_ri = (mags[1] + delta[1]) - (mags[2] + delta[2])
        col_rz = (mags[1] + delta[1]) - (mags[3] + delta[2])
        
        g2v_gr = -2.5*np.log(flux[0]/flux[1])
        g2v_ri = -2.5*np.log(flux[1]/flux[2])

        x_fit = gmr
        ri_param = [col_ri/col_gr,gr_ri_output.beta[1]]
        y_ri_fit = linfit_ri(ri_param, x_fit)
        rz_param = [col_rz/col_gr,gr_rz_output.beta[1]]
        y_rz_fit = linfit_rz(rz_param, x_fit)

        csq_ri = np.sum((rmi - y_ri_fit)**2/urmi**2)/ (len(gmr) - 3)
        
        csq_rz = np.sum(np.sum((rmz[good] - y_rz_fit[good])**2/urmz[good]**2) / (len(gmr) - 3))

        csq = 0.0
        csq = csq_rz + csq_ri

        #all results
        output_sizes_all.append(size)
        output_slope_all.append([col_gr,col_ri,col_rz])
        output_csq_all.append(csq)

        #if extinction slope is correct then ratios should be same number for a given material+ grain size
        if csq < csq_limit :
            print(size,csq,col_ri,col_rz)
            output_sizes.append(size)
            output_slope.append([col_gr,col_ri,col_rz])
            output_csq.append(csq)

        #if abs(slopegr_large - sloperi_large) < 0.05 and abs(slopegr_large - sloperz_large) < 0.05 :
        #    print("large", size,slopegr_large,sloperi_large,sloperz_large)
        
        #print("Ratio of observed gradient dx, dy to slope of grains: ",1.14 / col_gr, 0.51 / col_ri) #get numbers to match to fit general extinction slope
        #print("Size of grain, graidents",sgrain,3.896 / col_gr, 0.79 / col_ri) #get numbers to match to fit deepest extinction slope
    
    filename="asassn21qj_single_size_"+composition_neat+".json"

    output = {'species': composition,
              'name': composition_neat,
              'sizes':output_sizes_all,
              'slopes':output_slope_all,
              'csq':output_csq_all
              }
    with open('./'+filename, 'w') as payload:
        json.dump(output, payload)

    #Print results
    print("Grain sizes closely matching extinction:")
    print(output_sizes,output_csq)



    #Plot up results
    #all results
    fig, ax1 = plt.subplots()
    ax1.set_xlabel(r"Colour g$^{\prime}$ - r$^{\prime}$ (mag)",fontsize=16)
    ax1.set_ylabel(r"Colour r$^{\prime}$ - i$^{\prime}$ (mag)",fontsize=16)
    ax1.invert_yaxis()

    ax1.errorbar(gmr,rmi,xerr=ugmr,yerr=urmi,linestyle='',marker='o',mfc='blue',mec='darkblue',ecolor='darkblue',color='blue',alpha=0.5,label="gr-ri")
    for i in range(0,len(output_sizes_all)):
        scaling = 1.14 / output_slope_all[i][0]
        alpha = 0.25
        #if output_csq_all[i] < csq_limit: 
        #    alpha = 1.0
        ax1.plot([0.44,0.44+1.14],[0.11,0.11+scaling*output_slope_all[i][1]],next(linecycler),color='dodgerblue',alpha=alpha)
    best = np.argmin(output_csq_all)
    scaling = 1.14 / output_slope_all[best][0]
    ax1.plot([0.44,0.44+1.14],[0.11,0.11+scaling*output_slope_all[best][1]],linestyle='-',color='black',label=r"{:5.2f}$~\mu$m".format(output_sizes_all[best]),alpha=1.0)
    ax1.plot([0.44],[0.11],marker='*',markersize=24,mfc='yellow',mec='orange')
    ax1.annotate(composition_neat,(0.05,0.90),fontsize=20)
    ax1.set_xlim(0.0,2.0)
    ax1.set_ylim(0.0,1.0)
    ax1.tick_params(axis='both',labelsize=16)
    ax1.legend(loc='lower right')
    fig.savefig(direc + 'asassn21qj_lcogt_gr-ri_plot_'+composition+'_single_size_all.pdf',dpi=200)
    #plt.show()
    plt.tight_layout()
    plt.close()
    plt.clf()

    fig, ax1 = plt.subplots()
    good = np.where(rmz != -99)
    ax1.set_xlabel(r"Colour g$^{\prime}$ - r$^{\prime}$ (mag)",fontsize=16)
    ax1.set_ylabel(r"Colour r$^{\prime}$ - z$^{\prime}$ (mag)",fontsize=16)
    ax1.invert_yaxis()

    ax1.errorbar(gmr[good],rmz[good],xerr=ugmr[good],yerr=urmz[good],linestyle='',marker='o',mfc='red',mec='darkred',ecolor='red',color='red',alpha=0.5,label="gr-rz")
    for i in range(0,len(output_sizes_all)):
        alpha = 0.25
        #if output_csq_all[i] < csq_limit: 
        #    alpha = 1.0
        scaling = 1.14 / output_slope_all[i][0]
        ax1.plot([0.44,0.44+1.14],[0.14,0.14+scaling*output_slope_all[i][2]],next(linecycler),color='orange',alpha=alpha)
    best = np.argmin(output_csq_all)
    scaling = 1.14 / output_slope_all[best][0]
    ax1.plot([0.44,0.44+1.14],[0.14,0.14+scaling*output_slope_all[best][2]],linestyle='-',color='black',label=r"{:5.2f}$~\mu$m".format(output_sizes_all[best]),alpha=1.0)
    ax1.plot([0.44],[0.14],marker='*',markersize=24,mfc='yellow',mec='orange')
    ax1.annotate(composition_neat,(0.05,0.90),fontsize=20)
    ax1.set_xlim(0.0,2.0)
    ax1.set_ylim(0.0,1.0)
    ax1.tick_params(axis='both',labelsize=16)
    ax1.legend(loc='lower right')
    plt.tight_layout()
    fig.savefig(direc + 'asassn21qj_lcogt_gr-rz_plot_'+composition+'_single_size_all.pdf',dpi=200)
    #plt.show()
    plt.close()
    plt.clf()


    #results that vaugely match the reddening vectors
    if len(output_sizes) > 0:
        fig, ax1 = plt.subplots()
        ax1.set_xlabel(r"Colour g$^{\prime}$ - r$^{\prime}$ (mag)",fontsize=16)
        ax1.set_ylabel(r"Colour r$^{\prime}$ - i$^{\prime}$ (mag)",fontsize=16)
        ax1.invert_yaxis()

        ax1.errorbar(gmr,rmi,xerr=ugmr,yerr=urmi,linestyle='',marker='o',mfc='blue',mec='darkblue',ecolor='darkblue',color='blue',alpha=0.5,label="gr-ri")
        for i in range(0,len(output_sizes)):
            alpha = 0.25
            #if output_csq[i] < csq_limit: 
            #    alpha = 1.0
            scaling = 1.14 / output_slope[i][0]
            ax1.plot([0.44,0.44+1.14],[0.11,0.11+scaling*output_slope[i][1]],next(linecycler),color='dodgerblue',alpha=alpha)
        best = np.argmin(output_csq)
        scaling = 1.14 / output_slope[best][0]
        ax1.plot([0.44,0.44+1.14],[0.11,0.11+scaling*output_slope[best][1]],linestyle='-',color='black',label=r"{:5.2f}$~\mu$m".format(output_sizes[best]),alpha=1.0)
        ax1.plot([0.44],[0.11],marker='*',markersize=24,mfc='yellow',mec='orange')
        ax1.annotate(composition_neat,(0.05,0.90),fontsize=20)
        ax1.set_xlim(0.0,2.0)
        ax1.set_ylim(0.0,1.0)
        ax1.tick_params(axis='both',labelsize=16)
        ax1.legend(loc='lower right')
        plt.tight_layout()
        fig.savefig(direc + 'asassn21qj_lcogt_gr-ri_plot_'+composition+'_single_size.pdf',dpi=200)
        #plt.show()
        plt.close()
        plt.clf()

        fig, ax1 = plt.subplots()
        good = np.where(rmz != -99)
        ax1.set_xlabel(r"Colour g$^{\prime}$ - r$^{\prime}$ (mag)",fontsize=16)
        ax1.set_ylabel(r"Colour r$^{\prime}$ - z$^{\prime}$ (mag)",fontsize=16)
        ax1.invert_yaxis()

        ax1.errorbar(gmr[good],rmz[good],xerr=ugmr[good],yerr=urmz[good],linestyle='',marker='o',mfc='red',mec='darkred',ecolor='red',color='red',alpha=0.5,label="gr-rz")
        for i in range(0,len(output_sizes)):
            alpha = 0.25
            #if output_csq[i] < csq_limit: 
            #    alpha = 1.0
            scaling = 1.14 / output_slope[i][0]
            ax1.plot([0.44,0.44+1.14],[0.14,0.14+scaling*output_slope[i][2]],next(linecycler),color='orange',alpha=alpha)
        best = np.argmin(output_csq)
        scaling = 1.14 / output_slope[best][0]
        ax1.plot([0.44,0.44+1.14],[0.14,0.14+scaling*output_slope[best][2]],linestyle='-',color='black',label=r"{:5.2f}$~\mu$m".format(output_sizes[best]),alpha=1.0)
        ax1.plot([0.44],[0.14],marker='*',markersize=24,mfc='yellow',mec='orange')
        ax1.annotate(composition_neat,(0.05,0.90),fontsize=20)
        ax1.set_xlim(0.0,2.0)
        ax1.set_ylim(0.0,1.0)
        ax1.tick_params(axis='both',labelsize=16)
        ax1.legend(loc='lower right')
        plt.tight_layout()
        fig.savefig(direc + 'asassn21qj_lcogt_gr-rz_plot_'+composition+'_single_size.pdf',dpi=200)
        #plt.show()
        plt.close()
        plt.clf()

    index += 1
