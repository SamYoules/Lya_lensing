# Author: Sam Youles
# 31/5/19

import numpy as np
import healpy as hp
import sys
import glob
import os
import fitsio
from kappa_lya import *
from astropy.table import Table, join
import argparse


#-- Input arguments

if __name__=='__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Compute the convergence (kappa) between qsos and deltas.')

    parser.add_argument('--mapnumber', required=False, type=int, default=1, \
           help='index number of input map')
    parser.add_argument('--zcatin', required=True, type=str, \
           help='original quasar catalogue')
    parser.add_argument('--zcatdir', required=True, type=str, \
           help='folder containing lensed quasar catalogues')
    args, unknown = parser.parse_known_args()

    mapnumber = args.mapnumber
    zcatin = args.zcatin
    zcatdir = args.zcatdir
    zcatname = '{}/zcat_desi_drq_lensed_{}.fits'.format(zcatdir, mapnumber)

    #-- Create angular power spectrum of kappa
    theory = Theory()
    ell, cell = theory.get_cl_kappa(2.1, kmax=100., nz=100, lmax=10000)

    #-- Create kappa map
    nside=1024
    npix=nside**2*12
    seed=int(mapnumber)
    np.random.seed(seed)
    kappa = create_gaussian_kappa(ell, cell, nside=nside, seed=seed)

    #-- Write a new quasar catalogue with new RA & DEC for this input map
    t = Table.read(zcatin)
    theta_lens, phi_lens = kappa.displace_objects \
          (np.pi/2-np.radians(t['DEC']), np.radians(t['RA']))

    #-- Rename RA & DEC columns which hold unlensed position of qso
    t.rename_column('DEC', 'DEC0')
    t.rename_column('RA', 'RA0')

    #-- Add in new columns for lensed position of qso
    t['DEC'] = np.degrees(np.pi/2-theta_lens)
    t['RA'] = np.degrees(phi_lens)

    t.write(zcatname, format='fits', overwrite=True)

