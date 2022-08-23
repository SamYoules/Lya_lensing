#!/usr/bin/env python

'''Write power spectra for auto-cf, cross-cf and input-auto-cfs for 100 
realisations of estimated kappa, to be used for plotting errors. The 'cross' 
refers to being crossed with input map, not lya-qso.'''

from astropy.io import fits
import healpy as hp
import numpy as np
import sys

def get_correlations(filenumber, maptype, n):
    '''Open the estimated and input kappa files and get the auto and cross
       power spectrum'''

    print (f'{filenumber} of {n}')
    ##- Open maps
    filename = f'2.0-1000/Maps/Output/Gaussian/{maptype}{filenumber}.fits.gz'
    ifilename = f'Input_Maps/kappa_input{filenumber}.fits'
    kest = fits.open(filename)[1].data.kappa
    wkest = fits.open(filename)[1].data.wkappa
    kinput = fits.open(ifilename)[1].data.I

    ##- Get resolution of map
    nside=int(hp.npix2nside(kest.size))

    ##- Reset resolution of input maps to match estimated maps
    alm_in = hp.map2alm(kinput, lmax=3*nside)
    kinput = hp.alm2map(alm_in, nside=nside)

    ##- Mask area outside footprint
    mask = wkest!=0
    mask &= (wkest>np.percentile(wkest[mask], 0.1)) & \
                (wkest<np.percentile(wkest[mask], 99.9))
    kinput = kinput*(mask)+hp.UNSEEN*(~mask) 
    kinput[~mask]=hp.UNSEEN
    kest[~mask]=hp.UNSEEN

    ##- Get the power spectra from the maps
    Cl_auto = hp.sphtfunc.anafast(kest, lmax=3*nside-1)
    Cl_cross = hp.sphtfunc.anafast(kest, kinput, lmax=3*nside-1)

    Cl_input = hp.sphtfunc.anafast(kinput, lmax=3*nside-1)
    return Cl_auto, Cl_cross, Cl_input

##- Input (e.g.  auto, cross)
maptype = sys.argv[1]
N_files = sys.argv[2]

C_ells = []

for i in range(int(N_files)):
    j = str(i+1)
    C_ells.append(get_correlations(j, maptype, int(N_files)))

Cl_array = np.asarray(C_ells)
Cl_autos = Cl_array[:,0]
Cl_crosses = Cl_array[:,1]
Cl_inputs = Cl_array[:,2]
np.savetxt(f'2.0-1000/Maps/Cls/{maptype}/Cls.txt', Cl_autos)
np.savetxt(f'2.0-1000/Maps/Cls/{maptype}/Cls_crossed_with_input.txt', Cl_crosses)
np.savetxt('2.0-1000/Maps/Cls/Cl_inputs.txt', Cl_inputs)



