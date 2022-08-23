#!/usr/bin/env python

##-- Get list of pixels in ring around the blob and select just these deltas

import healpy as hp
import numpy as np
import glob
import shutil

#--blob centre is at phi=np.pi, theta=1.1
ra    = 180
dec   = np.degrees(np.pi/2. - 1.1)
nside = 16
vec   = hp.ang2vec(ra, dec, lonlat=True)
rad   = np.deg2rad(40)
hpix  = hp.query_disc(nside, vec, rad, inclusive=True)

alldeltas = glob.glob('*.fits.gz')
for d in alldeltas:
    h = d[6:-8]
    i = int(h)
    if i in hpix:
        newfile = f'/global/project/projectdirs/desi/users/syoules/Lensing/2.0-1000/Deltas/Unlensed_2/{d}'
        shutil.copyfile(d, newfile)


