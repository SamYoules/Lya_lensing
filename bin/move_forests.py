#!/usr/bin/env python
"""Lenses forests with a Gaussian lensing map

Reads delta files, finds new RA and DEC caused by lensing, and writes new delta
files. Uses seed (between 1 and 100) that corresponds to the input maps to 
produce bend angles with which to lens the deltas."""

import numpy as np
import healpy as hp
import glob
import os
import fitsio
from kappa_lya import *
import argparse

#-- Input arguments

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Move forests with bend angle (alpha) map.')

parser.add_argument('--indir', required=True, type=str, \
           help='folder containing deltas')
parser.add_argument('--outdir', required=True, type=str, \
           help='folder containing lensed deltas') 
parser.add_argument('--mapnumber', required=False, type=int, default=1, \
           help='index number of input map')
args, unknown = parser.parse_known_args()

indir = args.indir
outdir = args.outdir
mapnumber = args.mapnumber


#-- Create angular power spectrum of kappa
theory = Theory()
ell, cell = theory.get_cl_kappa(2.1, kmax=100., nz=100, lmax=10000)
nside=1024
npix=nside**2*12
seed=int(mapnumber)
np.random.seed(seed)
kappa = create_gaussian_kappa(ell, cell, nside=nside, seed=seed)


# Amend DEC and RA in each of the delta files by the bend angle from alpha map
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

        # Add bend angles to ra and dec
        theta_lens, phi_lens = kappa.displace_objects(np.pi/2-dec, ra) 
        
        # Rewrite new delta file with new values
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

