#!/usr/bin/env python

#-- Plot 2d residuals (drt_in - drt_out) as function of number of realisations
#-- Make 2 forests with N deltas, keep rt constant, vary rp for each delta by
#-- 1 Mpc/h, to model a range of r, build a covariance matrix covering each
#-- pair, use Cholesky decomposition to correlate the deltas, keep variance
#-- constant. Calculate drt_out for single realisation. Rinse & repeat.
#-- Write output to fits file.

import numpy as np
from numpy import linalg
import scipy as sp
from scipy import interpolate
import fitsio
import h5py
import argparse
import pylab as plt
import matplotlib.mlab as mlab
from matplotlib.offsetbox import AnchoredText
from astropy.table import Table
from multiprocessing import Pool,Process,Lock,Manager,cpu_count,Value
import psutil

def load_model(file_xi, file_fit, nbins=50) :
    '''Read xi file and fit file to make an interpolation object for the
       correlation between deltas, in bins of rp and rt.'''

    h = fitsio.FITS(file_xi)
    ff = h5py.File(file_fit, 'r')
    base = file_xi 
    fit = ff[base+'/fit'][...] 
    data_rp = h[1]['RP'][:]
    data_rt = h[1]['RT'][:]
    h.close()
    ff.close()

    rpmin = data_rp.reshape(50, 50)[0].max()
    rpmax = data_rp.reshape(50, 50)[-1].min()
    rtmin = data_rt.reshape(50, 50)[:, 0].max()
    rtmax = data_rt.reshape(50, 50)[:, -1].min()

    #-- create the regular grid for griddata
    rp = np.linspace(rpmin, rpmax, nbins)
    rt = np.linspace(rtmin, rtmax, nbins)

    xim = sp.interpolate.griddata((data_rt, data_rp), fit,
                (np.outer(np.ones(rp.size), rt).ravel(),
                np.outer(rp, np.ones(rt.size)).ravel()), method='cubic')

    #-- create interpolator object
    xi2d = sp.interpolate.RectBivariateSpline(rt, rp, \
               xim.reshape((50, 50)).T )

    return xi2d

def calc_drt(p):
    '''Calculate drt_out for p pairs of forests'''

    drt_out = np.zeros(nreal)

    #-- Construct Gaussian random variables to represent 2 forests
    g = np.random.randn(2*npix)

    #-- Apply transformation to correlate variables
    delta = L.dot(g)

    #-- Add pixel noise to deltas (with variance = var2)
    delta += np.random.randn(2*npix) * np.sqrt(var2)

    delta = np.reshape(delta, (2, npix))

    num = 0
    den = 0

    #-- create matrix of rp separations between each pixel pair (we assume
    #-- separation between each forest pixel is 1 Mpc/h)
    rp_i = np.arange(npix)
    rp_j = np.arange(npix)
    rp = abs(rp_i - rp_j[:, None])

    #-- Create matrix of each delta pair from forest0 and forest1
    d1d2 = delta[0] * delta[1][:, None] 
 
    #-- Mask off any pairs with separations >=100
    w = rp < 100
    rp = rp[w]
    d1d2 = d1d2[w]

    #-- Get correlation function and derivative for lensed separations
    xi_lens = xi2d(rt_lens, rp, grid=False)
    xip_lens = xi2d(rt_lens, rp, dx=1, grid=False)
 
    #-- Estimate drt, subtracting product of each pixel pair from xi
    num = np.sum( (xi_lens - d1d2) * xip_lens)
    den = np.sum( xip_lens**2 ) 
    drt_out = num / den

    return drt_out


if __name__=='__main__':

    parser = argparse.ArgumentParser(formatter_class =
               argparse.ArgumentDefaultsHelpFormatter,
               description='Check kappa estimator for lya-lya.')

    parser.add_argument('--drt_in', required=True, type=float, \
               help='drt value to model')
    parser.add_argument('--rt', required=True, type=float, \
               help='rt value to model')
    parser.add_argument('--nreal', required=True, type=int, default=1000000, \
               help='number of realisations')
    parser.add_argument('--npix', required=False, type=int, default=500, \
               help='number of pixels in forests')
    parser.add_argument('--var1', required=True, type=float, default=0.05, \
               help='lss variance of deltas (~ 0.04 at z=2 to 0.15 at z=3)')
    parser.add_argument('--var2', required=True, type=float, default=0.01, \
               help='pixel variance')
    parser.add_argument('--xi', required=True, \
               help='text file containing model')
    parser.add_argument('--fit', required=True, \
               help='text file containing correlation function')
    parser.add_argument('--nproc', required=False, type=int, default=1, \
               help='number of procs used in calculation')

    args, unknown = parser.parse_known_args()

    drt_in = args.drt_in
    rt_unl = args.rt
    nreal = args.nreal
    npix = args.npix
    var1 = args.var1
    var2 = args.var2
    xi_file = args.xi
    fit_file = args.fit
    nprocs = args.nproc

    #-- Apply lensing
    rt_lens = rt_unl + drt_in

    #-- Get xi_model interpolator object
    xi2d = load_model(xi_file, fit_file)

    #-- Construct a covariance matrix for a pair of forests of length npix.
    #-- We assume a distance of 1 Mpc/h in rp between each pixel
    #-- rt is zero for pixels in same forest, or rt_unl for pairs from
    #-- different forests
    
    C = np.zeros((2*npix, 2*npix), dtype=float)
    for i in range(2*npix):
        for j in range(i, 2*npix):
            if i==j:
                C[i,j] = var1
            else:
                rt_i = (i // npix)*rt_unl
                rt_j = (j // npix)*rt_unl
                rp_i = i % npix
                rp_j = j % npix 
                rt = abs(rt_i-rt_j)
                rp = abs(rp_i-rp_j)
                
                if rt > 100 or rp > 100:
                    C[i, j] = 0.
                    C[j, i] = 0.
                else:
                    xi = xi2d(rt, rp, grid=False)
                    C[i, j] = xi
                    C[j, i] = xi
                    
    #-- Cholesky decomposition, used to correlate variables, see 2.1, 1108.5606
    L = np.linalg.cholesky(C)

    #-- Use multiple processors    
    cpu_data = {}
    for p in range(nreal):
        cpu_data[p] = [p]

    pool = Pool(processes=nprocs)
    results = pool.map(calc_drt, cpu_data.values())
    pool.close()

    #-- Get mean and error on the mean             
    drt_out_mean  = np.mean(results)
    drt_out_error = np.std(results)/np.sqrt(nreal)
    resids = (np.array(results) - drt_in) / drt_in
    std = np.std(results)                                                        

    #-- Add to table of all results
    t = Table.read('drt/drtsv2.fits')
    t.add_row([drt_in, rt_unl, nreal, npix, var1, var2, drt_out_mean,
               drt_out_error])
    t.write('drt/drtsv2.fits', format='fits', overwrite=True)

    #-- Save results
    t1 = Table([results])
    t1.write(f'drt/drt_resultsv2_{drt_in}_rt{rt_unl}_var1{var1}_var2{var2}.fits',
               format='fits', overwrite=True)

    #-- Plot histogram with Gaussian fit
    fig, ax = plt.subplots(nrows=1, ncols=1)
    hist_res = plt.hist(results, bins=100)
    x = np.linspace(min(results), max(results), 1000)
    dx = hist_res[1][1] - hist_res[1][0]
    scale = len(results)*dx
    ax.plot(x, mlab.normpdf(x, drt_out_mean, std)*scale)

    anchored_text = AnchoredText(f"{nreal} realisations\n{npix} pixels\nstd={std:.2}\nmean = {drt_out_mean:.2}\nerror on mean={drt_out_error:.2}", loc=1,borderpad=0.,frameon=False)
    ax.add_artist(anchored_text)

    ax.set_xlabel('Estimated drt')                                                
    plt.title(f'Estimated drt: drt_in = {drt_in}, rt={rt_unl}')                            
    plt.savefig(f'plots/mc/hist_drt{drt_in}_rt{rt_unl}_npix{npix}_nr{nreal}.png')

