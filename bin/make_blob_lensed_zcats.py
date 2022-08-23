#!/usr/bin/env python

# SY 12/10/21
# Use Blob to move qsos and write lensed qso catalogues.
# RA and DEC are in degrees in the zcat (but Radians in the deltas)

import numpy as np
import healpy as hp
from astropy.table import Table, join
import kappa_lya

#-- Set resolution for lensing maps
nside=256
npix=nside**2*12

phi = np.array([np.pi-0.4, np.pi-0.3, np.pi-0.2, np.pi-0.1, np.pi, np.pi+0.1, np.pi+0.2, np.pi+0.3, np.pi+0.4, np.pi+0.5])
theta = np.array([1.0, 1.1, 1.2, 1.0, 1.1, 1.2, 1.0, 1.1, 1.2, 1.0])

for i in range(10):
    kappa = kappa_lya.create_blob_kappa(nside=nside, phi0=phi[i], theta0=theta[i])
    hp.fitsfunc.write_map(f'2.0-1000/Maps/Input/Large_Blob/kappa_blob{i}.fits', kappa.A, fits_IDL=False, overwrite=True)
    t = Table.read('2.0-1000/Catalog/zcat_drq.fits')
    theta_lens, phi_lens = kappa.displace_objects \
           (np.pi/2-np.radians(t['DEC']), np.radians(t['RA']))
    t.rename_column('DEC', 'DEC0')
    t.rename_column('RA', 'RA0')
    t['DEC'] = np.degrees(np.pi/2-theta_lens)
    t['RA'] = np.degrees(phi_lens)
    t.write(f'2.0-1000/Catalog/Large_Blob/zcat_drq_blob{i}.fits', format='fits', overwrite=True)
