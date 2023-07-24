#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 16:37:41 2022

@author: jonty
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from astropy.io import ascii
import time
import os
from scipy.interpolate import interp1d

direc = '/Users/jonty/mydata/asassn21qj/'

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

#read in LCO observations
lco_obs = ascii.read(direc+"asassn21qj_magnitudes.csv",comment="#",data_start=2) 

mjd = lco_obs["Time"]
mag = lco_obs["Mag"]
unc = lco_obs["MagErr"]
fil = lco_obs["Filter"]

sdssg = np.where(fil == 'gp')

mjd = mjd[sdssg]
mag = mag[sdssg]
unc = unc[sdssg]

#wavelengths of observation for the LCO filterbands griz
wave = [0.4770,0.6231,0.7625,0.9134] #um

#brightness of a G2V star in mJy
mtgt = 13.8 #unocculted magnitude of ASAS-SN 21qj
mags = np.asarray([mtgt,mtgt-0.44,mtgt-0.55,mtgt-0.58])
flux = 3631*10**(-0.4*np.asarray(mags))

distance = 567.*pc

#range of grain sizes for log-normal distributions
composition = "pyr-mg80"
size    = 0.7
width   = 0.5
rho     = 2.90 # gcm^-3

os.system("optool "+composition+" 0.7 -a 0.01 10.0 "+str(size)+":"+str(width)+" -l 0.2 2.0 -fmax 0")
optconst = ascii.read(direc+"dustkappa.dat",comment="#",data_start=2) 

wavelength = optconst["col1"].data
qabs = optconst["col2"].data
qsca = optconst["col3"].data

qext = (qabs + qsca)*1e1 # 10^4 m^2/ 10-3 kg : convert from cgs -> mks 



mgrain = (4./3.)*np.pi*(rho*1e3)*(size*um)**3 #kg
agrain = np.pi*(size*um)**2

f = interp1d(wavelength,qext)
qext = f(wave[0])
msum = 0.
mlo = 0.
mhi = 0.
msum_lo = 0.
msum_hi = 0.
for i in range(0,len(mag)):

    fobs = 3631*10**(-0.4*mag[i])
    fhi = 3631*10**(-0.4*(mag[i] - np.sqrt(unc[i]**2 + 0.03**2))) #NB because magnitudes are abhorrent to God.
    flo = 3631*10**(-0.4*(mag[i] + np.sqrt(unc[i]**2 + 0.03**2)))

    ext = ( fobs / flux[0] ) #np.e**(-1.*tau) #ratio I = Io * e**-tau , tau = sigma*N*l = cm2 cm-3 cm

    ext_lo = ( fhi / flux[0] ) #NB: high flux == low extinction, hence the switch here
    ext_hi = ( flo / flux[0] )

    tau = -1*np.log(ext)

    tau_lo = -1*np.log(ext_lo)
    tau_hi = -1*np.log(ext_hi)

    #tau = kappa * rho * x = Qext * (mass / volume) * distance = Qext * mass / rsol**2

    mclump = tau * np.pi * rsol**2 / qext
    
    print(mjd[i],mag[i],mclump)

    mlo = tau_lo * np.pi * rsol**2 / qext
    mhi = tau_hi * np.pi * rsol**2 / qext

    #print(mjd[i],mag[i],fobs/flux[0],tau,np.log10(mclump/MEarth))

    if(mjd[i] > 59476.) and (mjd[i] < 59508.): #range between two deepest eclipses
        time = 0.5 * (mjd[i+1] - mjd[i]) + 0.5*(mjd[i] - mjd[i-1])  
        msum += ((time*0.03918*au)/(2*rsol))*(mclump/MEarth) # number of independent stellar pencil beams between measurements
        msum_lo += ((time*0.03918*au)/(2*rsol))*(mlo/MEarth)
        msum_hi += ((time*0.03918*au)/(2*rsol))*(mhi/MEarth)

print("Mtotal between MJD 59476 and MJD 59508: ",msum," M_Earth or ",msum*MEarth," kg")
print("Lower mass: ",msum - msum_lo)
print("Upper mass: ",msum_hi - msum)