#!/usr/bin/env python
"""Lenses forests with a blob

Reads delta files, finds new RA and DEC caused by blob, 
and writes new delta files"""

from astropy.io import fits
import numpy as np
import healpy as hp
import sys
import glob
import os
import fitsio
import kappa_lya

#-- Directory names for delta files
indir  = sys.argv[1]
outdir = sys.argv[2]
blobno = sys.argv[3]

#nside=512
nside=256
npix=nside**2*12

#-- Create alpha map with single "blob" lens and write it to a fits file
kappa = kappa_lya.create_blob_kappa(nside=nside, phi0=np.pi+0.5, theta0=1.0)
#kappa.A /=10.
#hp.fitsfunc.write_map('kappa_blob.fits', kappa.A, fits_IDL=False, overwrite=True)
#hp.fitsfunc.write_map(f'2.0-1000/Maps/Input/kappa_blob{blobno}.fits', kappa.A, fits_IDL=False, overwrite=True)

#-- Amend DEC and RA in each of the delta files by the bend angle from map
alldeltas = glob.glob(indir+'/*.fits.gz')
ndel = len(alldeltas)
for i, filename in enumerate(alldeltas):
    hdus = fitsio.FITS(filename)
    print(i, ndel)

    out = fitsio.FITS(outdir+"/"+os.path.basename(filename),'rw',clobber=True)

    for hdu in hdus[1:]:
        header = hdu.read_header()
        ra = header['RA']
        dec = header['DEC']

        #-- Add bend angles to ra and dec
        theta_lens, phi_lens = kappa.displace_objects(np.pi/2-dec, ra) 
        
        #-- Rewrite new delta file with new values
        header['RA'] = phi_lens
        header['DEC'] = np.pi/2-theta_lens
        header['RA0'] = ra
        header['DEC0'] = dec
      
        #-- Re-create columns (maybe there's a better way to do this?) 
        ll = hdu['LOGLAM'][:]
        de = hdu['DELTA'][:]
        we = hdu['WEIGHT'][:]
        co = hdu['CONT'][:] 
        cols=[ll, de, we, co]
        names=['LOGLAM','DELTA','WEIGHT','CONT']
        out.write(cols, names=names, header=header, \
                  extname=str(header['THING_ID']))

    out.close()
