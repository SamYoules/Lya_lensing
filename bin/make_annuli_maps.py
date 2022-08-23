#!/usr/bin/env python

##-- Make a ring around the kappa blob and write it as a new file.

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

outdir = sys.argv[1]
band = sys.argv[2]
ring = int(sys.argv[3])

#--blob centre is at phi=np.pi, theta=1.1
ra      = 180
dec     = np.degrees(np.pi/2. - 1.1)
nside   = 256

#-- cross, abs R
#kblob=fits.open('maps/midpoint/blob/kappa_xnoiseless_negrp.fits.gz')[1].data['kappa']
#kw=fits.open('maps/midpoint/blob/kappa_xnoiseless_negrp.fits.gz')[1].data['wkappa']

#-- cross, bands of rp
#kblob=fits.open('maps/midpoint/blob/rp_bands/kappa_xnoiseless_band{}.fits.gz'.format(band))[1].data['kappa']
#kw=fits.open('maps/midpoint/blob/rp_bands/kappa_xnoiseless_band{}.fits.gz'.format(band))[1].data['wkappa']

#-- cross, original
kblob=fits.open('/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Output/blob_cross.fits.gz')[1].data['kappa']
kw=fits.open('/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Output/blob_cross.fits.gz')[1].data['wkappa']

#-- Auto
#kblob=fits.open('/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Output/blob_auto.fits.gz')[1].data['kappa']
#kw=fits.open('/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Output/blob_auto.fits.gz')[1].data['wkappa']

#--Input map
#kblob=fits.open('/global/cfs/projectdirs/desi/users/syoules/Lensing/kappa_blob.fits')[1].data['T']
#alm_in = hp.map2alm(kblob, lmax=3*nside)
#kblob = hp.alm2map(alm_in, nside=nside)

radius1 = int(ring * 2 - 2)
radius2 = int(ring * 2)
annulus = make_annulus(nside, ra, dec, radius1, radius2)
mask = []
kmap = kblob
for k,l in enumerate(kmap):
    mask.append(k in annulus)
mask=np.asarray(mask)
#print(np.sum(mask))
kmap[~mask]=0
kw[~mask]=0
    #hp.mollview(kmap, rot=(180,0,0), title='ring{}'.format(ring)) ##sy
    #plt.show() ##sy

##- Write map
mapout = '/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Output/{}/kappa{}_ring{}.fits.gz'.format(outdir, band, ring)
out = fitsio.FITS(mapout, 'rw', clobber=True)
head = {}
head['RAD1']=radius1
head['RAD2']=radius2
head['NSIDE']=nside
out.write([kmap, kw], names=['kappa', 'wkappa'], header=head)
out.close()

hp.gnomview(kmap, rot=[-180, 30], min=-5, max=5, reso=20)
plt.savefig(f'/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Maps/Output/Cross/gnomview_kappa{band}_ring{ring}.png')
