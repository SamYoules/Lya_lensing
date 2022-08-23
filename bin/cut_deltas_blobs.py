#!/usr/bin/env python

## Use a series of lensing maps, each containing a single large Gaussian lens
## Each lens is in a slightly different position so we can use these maps
## to simulate different realisations with the same mocks.
## Get a list of pixels in a ring around each blob and select just these deltas
## Write them to file

import glob
import shutil
import healpy as hp
import numpy as np
nside = 16
alldeltas = glob.glob('*.fits.gz')
phi = np.array([np.pi-0.4, np.pi-0.3, np.pi-0.2, np.pi-0.1, np.pi, np.pi+0.1, np.pi+0.2, np.pi+0.3, np.pi+0.4, np.pi+0.5])
theta = np.array([1.0, 1.1, 1.2, 1.0, 1.1, 1.2, 1.0, 1.1, 1.2, 1.0])

ra = phi*1
dec = np.pi/2 - theta
mock_index = np.arange(phi.size)

for i in mock_index:
    print(i)
    rai   = np.degrees(ra[i])
    deci  = np.degrees(dec[i])
    vec   = hp.ang2vec(rai, deci, lonlat=True)
    rad   = np.deg2rad(40)
    hpix  = hp.query_disc(nside, vec, rad, inclusive=True)
    for d in alldeltas:
        h = d[6:-8]
        j = int(h)
        if j in hpix:
            newfile = f'../Large_Blob/Unlensed_{i}/{d}'
            shutil.copyfile(d, newfile)

