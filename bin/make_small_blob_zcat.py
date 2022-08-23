#!/usr/bin/env python

# SY 3/11/21
# Use small blob to move qsos and write lensed qso catalogues.
# RA and DEC are in degrees in the zcat (but Radians in the deltas)

import numpy as np
import healpy as hp
from astropy.table import Table, join
import kappa_lya

#-- Set resolution for lensing maps
nside=512
npix=nside**2*12

#-- Create alpha map with single "blob" lens
kappa = kappa_lya.create_small_blob_kappa(nside=nside, phi0=np.pi, theta0=1.1)

#-- Amend DEC and RA for each qso by the bend angles (alphas)
t = Table.read('2.0-1000/Catalog/zcat_drq.fits')
theta_lens, phi_lens = kappa.displace_objects \
           (np.pi/2-np.radians(t['DEC']), np.radians(t['RA']))

#-- Rename RA & DEC columns which hold unlensed position of qso
t.rename_column('DEC', 'DEC0')
t.rename_column('RA', 'RA0')

#-- Add in new columns for lensed position of qso
t['DEC'] = np.degrees(np.pi/2-theta_lens)
t['RA'] = np.degrees(phi_lens)

t.write('2.0-1000/Catalog/zcat_drq_small_blob.fits', format='fits', overwrite=True)
