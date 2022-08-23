#!/usr/bin/env python

"""Make kappa map from weighted average of auto and cross maps"""

import numpy as np
import fitsio
from astropy.io import fits
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Compute kappa for LYAxLYA using midpoint.')

parser.add_argument('--map1',
                    required=True,
                    type=str,
                    help='pathname to kappalist fits file')

parser.add_argument('--map2',
                    required=True,
                    type=str,
                    help='pathname to kappalist fits file') 

parser.add_argument('--mapout',
                    required=True,
                    type=str,
                    help='pathname to output fits file')

args, unknown = parser.parse_known_args()

map1   = args.map1
map2   = args.map2
mapout = args.mapout

##- Read maps (kappa and weight)
khead = fits.open(map1)[1].header
k1 = fits.open(map1)[1].data.kappa
w1 = fits.open(map1)[1].data.wkappa
k2 = fits.open(map2)[1].data.kappa
w2 = fits.open(map2)[1].data.wkappa

##- Calculate weighted average
kap = np.zeros(w1.size)
wkap = np.zeros(w1.size)

for i in range(w1.size):
    if (w1[i]==0):
        kap[i] = 0.
        wkap[i] = 0.
    else:
        kap[i] = (k1[i]*w1[i] + k2[i]*w2[i])/(w1[i] + w2[i])
        wkap[i] = w1[i] + w2[i]

##- Write map
out = fitsio.FITS(mapout, 'rw', clobber=True)
head = {}
head['RPMIN']=khead['RPMIN']
head['RPMAX']=khead['RPMAX']
head['RTMIN']=khead['RTMIN']
head['RTMAX']=khead['RTMAX']
head['NT']=khead['NT']
head['NP']=khead['NP']
head['NSIDE']=khead['NSIDE']
out.write([kap, wkap], names=['kappa', 'wkappa'], header=head)
out.close()

