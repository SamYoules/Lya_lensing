#!/usr/bin/env python

##-- Plot mean kappa as a function of distance from the centre of the blob.
##- Average kappa in annuli around the blobs (10 blobs in different locations, so we need separate annuli for each one.

import healpy as hp
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import fitsio
import sys


def make_annulus(nside, ra, dec, r1, r2):
    '''Make two discs of different sizes around the blob, and then subtract the
       contents of the inner disc, to return an array of all the healpixel
       indices in the resulting ring.'''

    ra   = np.degrees(ra)
    dec  = np.degrees(dec)
    vec  = hp.ang2vec(ra, dec, lonlat=True)
    rad1 = np.deg2rad(r1)
    rad2 = np.deg2rad(r2)
    hpix_inner = hp.query_disc(nside, vec, rad1, inclusive=True)
    hpix_outer = hp.query_disc(nside, vec, rad2, inclusive=True)
    hpix  = np.asarray(list(set(hpix_outer).difference(hpix_inner)))
    return(hpix)

def make_input():
    '''Make a file containing the mean values of annuli around the input map'''

    kinput = fits.open('/global/cfs/projectdirs/desi/users/syoules/Lensing/kappa_blob.fits')[1].data['T']

    nside = 256
    mean_input = []
    x = np.arange(1,21)*2
    ra = np.pi
    dec = np.pi/2 - 1.1
    for r in range(1,21):
        print(r)
        radius1 = int(r * 2 - 2)
        radius2 = int(r * 2)
        annulus = make_annulus(nside, ra, dec, radius1, radius2)
        mask = []
        for k,l in enumerate(kinput):
            mask.append(k in annulus)
        mask = np.asarray(mask)
        mean_input.append(np.mean(kinput[mask]))

    file = open("/global/cfs/projectdirs/desi/users/syoules/Lensing/mean_blob_input", "wb")
    np.save(file, mean_input)
    file.close


nside = 256
x = np.arange(1,21)*2

#-- Read mean_input
file = open("/global/cfs/projectdirs/desi/users/syoules/Lensing/mean_blob_input", "rb")
mean_input = np.load(file)
file.close

plt.figure()
plt.plot(x, mean_input, label="Input", color='k')

#-- Get auto/cross blob files and their locations in ra and dec
#ra = [np.pi-0.4, np.pi-0.3, np.pi-0.2, np.pi-0.1, np.pi, np.pi+0.1, np.pi+0.2, np.pi+0.3, np.pi+0.4, np.pi+0.5]
#dec = [np.pi/2 - 1.0, np.pi/2 - 1.1, np.pi/2 - 1.2, np.pi/2 - 1.0, np.pi/2 - 1.1, np.pi/2 - 1.2, np.pi/2 - 1.0, np.pi/2 - 1.1, np.pi/2 - 1.2, np.pi/2 - 1.0]

ra = [np.pi-0.2, np.pi, np.pi+0.2, np.pi-0.3, np.pi-0.1, np.pi+0.1, np.pi+0.3, np.pi-0.2, np.pi, np.pi+0.2]
dec = [1.1 - np.pi/2, 1.1 - np.pi/2, 1.1 - np.pi/2, 1.3 - np.pi/2, 1.3 - np.pi/2, 1.3 - np.pi/2, 1.3 - np.pi/2, 1.5 - np.pi/2, 1.5 - np.pi/2, 1.5 - np.pi/2]

sum_auto_all  = np.zeros(20)
sum_cross_all = np.zeros(20)
sum_w_auto_all = np.zeros(20)
sum_w_cross_all = np.zeros(20)

#-- For each realisation of the blob
for i in range(10):
    fa = f'/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Output/Large_Blob/blob_auto{i}.fits.gz'
    fc = f'/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Output/Large_Blob/blob_cross{i}.fits.gz'
    kauto      = fits.open(fa)[1].data['kappa']
    kcross     = fits.open(fc)[1].data['kappa']
    kwa        = fits.open(fa)[1].data['wkappa']
    kwc        = fits.open(fc)[1].data['wkappa']
    sum_auto   = []
    sum_cross  = []
    weighta    = []
    weightc    = []

    #-- For each annulus
    for r in range(1,21):
        print(r)
        radius1 = int(r * 2 - 2)
        radius2 = int(r * 2)
        annulus = make_annulus(nside, ra[i], dec[i], radius1, radius2)
        mask = []
        for k,l in enumerate(kcross):
            mask.append(k in annulus)
        mask = np.asarray(mask)
        sum_cross.append(np.sum( kcross[mask] * kwc[mask] ))
        sum_auto.append(np.sum( kauto[mask] * kwa[mask] ))
        weightc.append(np.sum(kwc[mask]))
        weighta.append(np.sum(kwa[mask]))

    sum_auto  = np.asarray(sum_auto)
    sum_cross = np.asarray(sum_cross)
    weighta   = np.asarray(weighta)
    weightc   = np.asarray(weightc)

    ka = sum_auto/weighta
    kc = sum_cross/weightc
    #plt.plot(x, kc, alpha=0.4, label=i)
    plt.plot(x, ka, color="blue", alpha=0.3)
    plt.plot(x, kc, color="orange", alpha=0.3)
    sum_auto_all += (sum_auto * weighta)
    sum_cross_all += (sum_cross * weightc)
    sum_w_auto_all += (weighta)
    sum_w_cross_all += (weightc)

kaa = sum_auto_all/sum_w_auto_all
kca = sum_cross_all/sum_w_cross_all
plt.plot(x, kaa, label="Auto", color="blue")
plt.plot(x, kca, label="Cross", color="orange")
plt.xlabel('Radius [deg]')
plt.ylabel("kappa")
plt.legend()
plt.savefig("/global/cfs/projectdirs/desi/users/syoules/Lensing/Plots/mean_blob.pdf")
