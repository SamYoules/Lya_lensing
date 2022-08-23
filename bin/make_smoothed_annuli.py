#!/usr/bin/env python

##-- Average kappa in annuli around the blob

import healpy as hp
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import fitsio
import scipy as sp
import sys

def make_annulus(nside, ra, dec, r1, r2):
    "ra, dec and radii are in degrees. Makes two discs of different sizes around the blob, and then subtracts the contents of the inner disc, to return an array of all the healpixels in the resulting ring."

    vec = hp.ang2vec(ra, dec, lonlat=True)
    rad1 = np.deg2rad(r1)
    rad2 = np.deg2rad(r2)
    hpix_inner = hp.query_disc(nside, vec, rad1, inclusive=True)
    hpix_outer = hp.query_disc(nside, vec, rad2, inclusive=True)
    hpix = np.asarray(list(set(hpix_outer).difference(hpix_inner)))
    return(hpix)

#--blob centre is at phi=np.pi, theta=1.1
ra      = 180
dec     = np.degrees(np.pi/2. - 1.1)
nside   = 256

#-- cross, original
#kblob=fits.open('/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Output/blob_cross.fits.gz')[1].data['kappa']
#kblob=fits.open('/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Output/blob_auto.fits.gz')[1].data['kappa']
kblob=fits.open('/global/cfs/projectdirs/desi/users/syoules/Lensing/kappa_blob.fits')[1].data['T']


kmap = np.zeros(kblob.shape)
ring = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
for r in ring:
    print(r)
    radius1 = int(r * 2 - 2)
    radius2 = int(r * 2)
    annulus = make_annulus(nside, ra, dec, radius1, radius2)
    mask = []
    for k,l in enumerate(kblob):
        mask.append(k in annulus)
    mask = np.asarray(mask)
    kmap[mask] = np.mean(kblob[mask])


##- Write map
#mapout = r'/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Output/Cross/kappa_mean_annuli.fits.gz'
#mapout = r'/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Output/Auto/kappa_mean_annuli.fits.gz'
mapout = r'/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Input/kappa_mean_annuli.fits.gz'
out = fitsio.FITS(mapout, 'rw', clobber=True)
head = {}
head['RAD1']=radius1
head['RAD2']=radius2
head['NSIDE']=nside
out.write([kmap], names=['kappa'], header=head)
out.close()

hp.gnomview(kmap, rot=[-180, 30], reso=20)
#plt.savefig(f'/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Output/Cross/gnomview_kappa_mean_annuli.png')
#plt.savefig(f'/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Output/Auto/gnomview_kappa_mean_annuli.png')
plt.savefig(f'/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Input/gnomview_kappa_mean_annuli.png')

hp.mollview(kmap, rot=[-180, 0])
#plt.savefig(f'/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Output/Cross/mollview_kappa_mean_annuli.png')
#plt.savefig(f'/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Output/Auto/mollview_kappa_mean_annuli.png')
plt.savefig('/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Input/mollview_kappa_mean_annuli.png')