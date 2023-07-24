#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 16:20:29 2022

@author: jonty
"""

import numpy as np
import glob
from astropy.io import ascii, fits
from astropy import units as u 
from astropy.wcs import WCS
import subprocess
import matplotlib.pyplot as plt
import datetime
from photutils.aperture import CircularAperture, CircularAnnulus
from photutils import aperture_photometry

zpnt = [25.897,25.675,25.112,24.33] #griz zero points (mags)
kcor = [0.20,0.09,0.02,0.02] #griz airmass corrections 

#Definitions
direc = '/Users/jonty/mydata/asassn21qj/'

def convert_coords(coords):
    """
    
    Parameters
    ----------
    coords : Str list
        String list of coordinates in hms dms 
    
    Returns
    -------
    ra : flt list
        Right Ascension in decimal degrees
    dec : flt list
        Declination in decimal degrees
    """
    
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    
    c = SkyCoord(coords,unit=(u.hourangle, u.deg))
    ra = c.ra.degree
    dec = c.dec.degree
    
    return ra, dec

#First check for new observations
check_data = True
if check_data == True: 
    subprocess.call("python download_data.py -sdate 2022-03-01 -edate 2023-04-01", shell=True)

#Then read in the sample properties.
target = "Asassn-21qj"
radec  = "08 15 23.3 -38 59 23.3"
mtgt   = 13.80
ra,dec = convert_coords(radec)
print(ra,dec)
#List of subdirectories in raw/ folder containing observations
data_direc = direc + 'raw/'
dates = glob.glob(data_direc + "*/")

#Next loop over all targets in the sample, and all days for which there are 
#observations to create the light curves.

tobs = []
mobs = []
uobs = []
filt = []

for date in dates:
    #Identify all the files in each date folder, open them iteratively.
    files = glob.glob(date + '*.fits.fz')
    #print(date)
    #Open each fits image in turn, check to see if target names match
    #then append date, magnitude, and uncertainty to list of observations
    #for the target.
    for file in files:
        #print(file)
        hdul = fits.open(file)
        
        w    = WCS(hdul["SCI"].header)
        data = hdul['SCI'].data
        erro = hdul['ERR'].data
        #print(hdul.info())
        #print(hdul[1].header)
        obsblock = hdul["SCI"].header['OBJECT']
        #print(target,obsblock)            
        try:
            #Read in observations catalogue for the image (BANZAI)
            cat = hdul["CAT"].data
            
            #Find closest match to target position within the field
            cat_ra  = np.asarray(cat['ra'])
            cat_dec = np.asarray(cat['dec'])

            #Do magnitude values exist? 
            cat_mag = cat['mag']
            cat_magerr = cat['magerr']

            separation = 3600.*np.sqrt((cat_ra*np.cos((np.pi/180.)*cat_dec) - ra*np.cos((np.pi/180.)*dec))**2 + (cat_dec - dec)**2)
            nearest = np.argmin(separation)
            
            #print(separation[nearest])
            #Append values to list of observations
            #print(hdul["SCI"].header['MJD-OBS'],separation[nearest],cat['mag'][nearest],cat['magerr'][nearest])
            #print(file)
            if separation[nearest] < 2.0 : 
                tobs.append(float(hdul["SCI"].header['MJD-OBS']))
                mobs.append(float(cat['mag'][nearest]))
                uobs.append(float(cat['magerr'][nearest]))
                filt.append(hdul["SCI"].header['FILTER'])

            if separation[nearest] >= 2.0 and separation[nearest] < 60.0:
                print("No target found in catalogue at location, doing aperture photometry")
                src_coords = np.asarray([[ra,dec]])
                pix = w.all_world2pix(src_coords,1)
                xy = (pix[0][0],pix[0][1])
                
                aperture = CircularAperture(xy, r=5.)
                annulus  = CircularAnnulus(xy, r_in=15., r_out=20.)
                
                apers = [aperture, annulus]
                
                phot = aperture_photometry(data,apers,error=erro)
                
                bkg_mean = phot['aperture_sum_1'][0] / annulus.area
                bkg_sum = bkg_mean * aperture.area
                
                counts_bkgsub = phot['aperture_sum_0'][0] - bkg_sum
                
                annulus_masks = annulus.to_mask(method='center')
                annulus_data = annulus_masks.multiply(data)
                mask = annulus_masks.data
                annulus_data_1d = annulus_data[mask > 0]
                
                counts_stdev  = np.sqrt(phot['aperture_sum_err_0'][0]**2 + phot['aperture_sum_err_1'][0]**2)
                
                texp = hdul["SCI"].header['EXPTIME']
                amss = hdul["SCI"].header['AIRMASS']
                thjd = hdul["SCI"].header['MJD-OBS']
                band = hdul["SCI"].header['FILTER']

                if band == 'gp':
                    mag = -2.5*np.log10(counts_bkgsub/texp) + (amss*kcor[0]) + zpnt[0]
                    magerr = mag - (-2.5*np.log10((counts_bkgsub+counts_stdev)/texp) + (amss*kcor[0]) + zpnt[0]) 

                if band == 'rp':
                    mag = -2.5*np.log10(counts_bkgsub/texp) + (amss*kcor[1]) + zpnt[1]
                    magerr = mag - (-2.5*np.log10((counts_bkgsub+counts_stdev)/texp) + (amss*kcor[1]) + zpnt[1])

                if band == 'ip':
                    mag = -2.5*np.log10(counts_bkgsub/texp) + (amss*kcor[2]) + zpnt[2]
                    magerr = mag - (-2.5*np.log10((counts_bkgsub+counts_stdev)/texp) + (amss*kcor[2]) + zpnt[2]) 

                if band == 'zs':
                    mag = -2.5*np.log10(counts_bkgsub/texp) + (amss*kcor[3]) + zpnt[3]
                    magerr = mag - (-2.5*np.log10((counts_bkgsub+counts_stdev)/texp) + (amss*kcor[3]) + zpnt[3])

                tobs.append(float(hdul["SCI"].header['MJD-OBS']))
                mobs.append(mag)
                uobs.append(magerr)
                filt.append(hdul["SCI"].header['FILTER'])
            
            if separation[nearest] > 60.:
                print("Nearest target > 60arcsec away, misaligned image")  

        except:
            print("No catalogue found for target ",target," in ",file," doing aperture photometry")
            src_coords = np.asarray([[ra,dec]])
            pix = w.all_world2pix(src_coords,1)
            xy = (pix[0][0],pix[0][1])
            
            aperture = CircularAperture(xy, r=5.)
            annulus  = CircularAnnulus(xy, r_in=15., r_out=20.)
            
            apers = [aperture, annulus]
            
            phot = aperture_photometry(data,apers,error=erro)
            
            bkg_mean = phot['aperture_sum_1'][0] / annulus.area
            bkg_sum = bkg_mean * aperture.area
            
            counts_bkgsub = phot['aperture_sum_0'][0] - bkg_sum
            
            annulus_masks = annulus.to_mask(method='center')
            annulus_data = annulus_masks.multiply(data)
            mask = annulus_masks.data
            annulus_data_1d = annulus_data[mask > 0]
            
            counts_stdev  = np.sqrt(phot['aperture_sum_err_0'][0]**2 + phot['aperture_sum_err_1'][0]**2)
            
            texp = hdul["SCI"].header['EXPTIME']
            amss = hdul["SCI"].header['AIRMASS']
            thjd = hdul["SCI"].header['MJD-OBS']
            band = hdul["SCI"].header['FILTER']

            if band == 'gp':
                mag = -2.5*np.log10(counts_bkgsub/texp) + (amss*kcor[0]) + zpnt[0]
                magerr = mag - (-2.5*np.log10((counts_bkgsub+counts_stdev)/texp) + (amss*kcor[0]) + zpnt[0]) 

            if band == 'rp':
                mag = -2.5*np.log10(counts_bkgsub/texp) + (amss*kcor[1]) + zpnt[1]
                magerr = mag - (-2.5*np.log10((counts_bkgsub+counts_stdev)/texp) + (amss*kcor[1]) + zpnt[1])

            if band == 'ip':
                mag = -2.5*np.log10(counts_bkgsub/texp) + (amss*kcor[2]) + zpnt[2]
                magerr = mag - (-2.5*np.log10((counts_bkgsub+counts_stdev)/texp) + (amss*kcor[2]) + zpnt[2]) 

            if band == 'zs':
                mag = -2.5*np.log10(counts_bkgsub/texp) + (amss*kcor[3]) + zpnt[3]
                magerr = mag - (-2.5*np.log10((counts_bkgsub+counts_stdev)/texp) + (amss*kcor[3]) + zpnt[3])

            tobs.append(float(hdul["SCI"].header['MJD-OBS']))
            mobs.append(mag)
            uobs.append(magerr)
            filt.append(hdul["SCI"].header['FILTER']) 
            #print(file,thjd,mag,magerr)

#sort values into chronological order
tobs = np.asarray(tobs)
mobs = np.asarray(np.nan_to_num(mobs))
uobs = np.asarray(np.nan_to_num(uobs))
filt = np.asarray(filt)

index = tobs.argsort()

tobs = tobs[index]
mobs = mobs[index]
uobs = uobs[index]
filt = filt[index]

#for i in range(len(tobs)):
#    print(tobs[i],mobs[i],uobs[i],filt[i])

#Write out table
from astropy.table import Table
t = Table([tobs, mobs, uobs, filt],names=('Time', 'Mag', 'MagErr','Filter'))
t.meta['comments'] = ['Data from LCO Program DDT2021B-003','PI S. Ertel','Time in MJD','Source brightness and uncertaintiy in Magnitudes','Filters SDSS gp,rp,ip and PANSTARRS Z']
t.write('asassn21qj_magnitudes.csv', format='csv', overwrite=True,comment='#')

#divide into filters
gp = np.where(filt == 'gp')
rp = np.where(filt == 'rp')
ip = np.where(filt == 'ip')
zs = np.where(filt == 'zs')

#Plot up results
fig, ax1 = plt.subplots()
ax1.set_xlabel(r"time (MJD)")
ax1.set_ylabel(r"$g^{\prime}$ (mag)")
ax1.invert_yaxis()
ax1.errorbar(tobs[gp],mobs[gp],xerr=None,yerr=uobs[gp],linestyle='',marker='o',mfc='gray',mec='darkgray',ecolor='gray',color='gray',alpha=0.5,label="sdss gp")
ax1.errorbar(tobs[rp],mobs[rp],xerr=None,yerr=uobs[rp],linestyle='',marker='o',mfc='red',mec='darkred',ecolor='red',color='red',alpha=0.5,label="sdss rp")
ax1.errorbar(tobs[ip],mobs[ip],xerr=None,yerr=uobs[ip],linestyle='',marker='o',mfc='dodgerblue',mec='blue',ecolor='blue',color='blue',alpha=0.5,label="sdss ip")
ax1.errorbar(tobs[zs],mobs[zs],xerr=None,yerr=uobs[zs],linestyle='',marker='o',mfc='pink',mec='purple',ecolor='purple',color='purple',alpha=0.5,label="pstar zs")
ax1.plot([np.min(tobs),np.max(tobs)],[mtgt,mtgt],linestyle=':',marker='',color='black',alpha=1.0)
ax1.set_xlim(np.min(tobs) - 10,np.max(tobs) + 10)
ax1.set_ylim(21.,13.)
ax1.legend(loc='lower right')
fig.savefig(direc + target + '_lcogt_time_series_plot.png',dpi=200)
plt.show()
plt.close()
plt.clf()

#calculate gr and ri colour pairs
tcol = []
mg = []
ug = []
mr = []
ur = []
mi = []
ui = []
mz = []
uz = []
gmr = []
ugmr = []
rmi = []
urmi = []
gmz = []
ugmz = []
rmz = []
urmz = []

for i in range(1,len(tobs)-2):
    if filt[i] == 'rp' and filt[i-1] == 'gp' and filt[i+1] == 'ip':
        tcol.append(tobs[i])
        mg.append(mobs[i-1])
        mr.append(mobs[i])
        mi.append(mobs[i+1])
        ug.append(uobs[i-1])
        ur.append(uobs[i])
        ui.append(uobs[i+1])
        
        gmr.append(mobs[i-1] - mobs[i])
        ugmr.append(np.sqrt(uobs[i-1]**2 + uobs[i]**2))
        rmi.append(mobs[i] - mobs[i+1])
        urmi.append(np.sqrt(uobs[i]**2 + uobs[i+1]**2))
        if i+2 < len(filt) and filt[i+2] == 'zs':
            mz.append(mobs[i+2])
            uz.append(uobs[i+2])
            gmz.append(mobs[i-1] - mobs[i+2])
            ugmz.append(np.sqrt(uobs[i-1]**2 + uobs[i+2]**2))
            rmz.append(mobs[i] - mobs[i+2])
            urmz.append(np.sqrt(uobs[i]**2 + uobs[i+2]**2))
        else:
            mz.append(-99)
            uz.append(-99)
            gmz.append(-99)
            ugmz.append(-99)
            rmz.append(-99)
            urmz.append(-99)
            

tt = Table([tcol,mg,ug,mr,ur,mi,ui,mz,uz, gmr, ugmr, rmi, urmi, gmz, ugmz, rmz, urmz],names=('Time', 'mg','ug','mr','ur','mi','ui','mz','uz','g-r', 'Err g-r','r-i','Err r-i', 'g-z','Err g-z', 'r-z', 'Err r-z'))
tt.meta['comments'] = ['Data from LCO Program DDT2021B-003','PI S. Ertel','Time in MJD','Colours and uncertainties in Magnitudes','Filters SDSS gp,rp,ip']
tt.write('asassn21qj_colours.csv', format='csv', overwrite=True,comment='#')

#Plot up results
fig, ax1 = plt.subplots()
ax1.set_xlabel(r"Colour g$^{\prime}$ - r$^{\prime}$ (mag)")
ax1.set_ylabel(r"Colour r$^{\prime}$ - i$^{\prime}$ (mag)")
ax1.invert_yaxis()
ax1.errorbar(gmr,rmi,xerr=ugmr,yerr=urmi,linestyle='',marker='o',mfc='red',mec='darkred',ecolor='red',color='red',alpha=0.5,label="colour")
ax1.set_xlim(0.0,5.0)
ax1.set_ylim(0.0,1.0)
ax1.legend(loc='lower right')
fig.savefig(direc + target + '_lcogt_colour_plot.png',dpi=200)
plt.show()
plt.close()
plt.clf()

#Plot up results
fig, ax1 = plt.subplots()
ax1.set_xlabel(r"Colour g$^{\prime}$ - r$^{\prime}$ (mag)")
ax1.set_ylabel(r"Colour r$^{\prime}$ - z$^{\prime}$ (mag)")
ax1.invert_yaxis()
ax1.errorbar(gmr,rmz,xerr=ugmr,yerr=urmz,linestyle='',marker='o',mfc='pink',mec='purple',ecolor='purple',color='purple',alpha=0.5,label="colour")
ax1.set_xlim(0.0,5.0)
ax1.set_ylim(0.0,1.5)
ax1.legend(loc='lower right')
fig.savefig(direc + target + '_lcogt_colour_plot2.png',dpi=200)
plt.show()
plt.close()
plt.clf()