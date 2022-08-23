#!/usr/bin/env python
"""Lenses forests with several different blobs

Reads delta files, finds new RA and DEC caused by blob, 
and writes new delta files"""

from astropy.io import fits
import numpy as np
import healpy as hp
import glob
import os
import fitsio
import kappa_lya

nside=256
npix=nside**2*12

phi = np.array([np.pi-0.4, np.pi-0.3, np.pi-0.2, np.pi-0.1, np.pi, np.pi+0.1, np.pi+0.2, np.pi+0.3, np.pi+0.4, np.pi+0.5])
theta = np.array([1.0, 1.1, 1.2, 1.0, 1.1, 1.2, 1.0, 1.1, 1.2, 1.0])

#ra = phi*1
#dec = np.pi/2 - theta
mock_index = np.arange(phi.size)

for j in mock_index:
    din  = f'/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Deltas/Large_Blob/Unlensed_{j}'
    dout = f'/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Deltas/Large_Blob/Blob{j}'
    kappa = kappa_lya.create_blob_kappa(nside=nside, phi0=phi[j], theta0=theta[j])
    #hp.fitsfunc.write_map('kappa_blob.fits', kappa.A, fits_IDL=False)

    #-- Amend DEC and RA in each of the delta files by the bend angle from map
    alldeltas = glob.glob(din+'/*.fits.gz')
    ndel = len(alldeltas)
    for i, filename in enumerate(alldeltas):
        hdus = fitsio.FITS(filename)
        print(i, ndel)

        out = fitsio.FITS(dout+"/"+os.path.basename(filename),'rw',clobber=True)

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
