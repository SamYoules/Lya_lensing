#!/usr/bin/env python

"""Compute estimated kappa for midpoint between qsos using estimator for
   lya x lya, and write to file"""

import numpy as np
import scipy as sp
import fitsio, glob, pickle, time
from astropy.io import fits
import argparse
import healpy as hp
import sys
from scipy import interpolate
import copy
from picca import constants
from picca.data_lens import Delta
from picca import pickle_deltas as pkl
from multiprocessing import Pool,Process,Lock,Manager,cpu_count,Value
import h5py
import psutil

class kappa:

    nside = 256
    nside_data = 32
    rot = hp.Rotator(coord=['C', 'G'])
    lambda_abs = 1215.67

    fid_Om=0.315
    cosmo = constants.Cosmo(fid_Om)
    h = 1.
    rt_min=0.
    rt_max=40.
    rp_min=0.
    rp_max=10.
    nt = 1
    np = 1

    xi2d = None

    angmax = None
    z_min_pix = None

    true_corr = False

    counter = Value('i',0)
    lock = Lock()

    data={}
    ndata=0


    @staticmethod
    def load_model(file_xi, file_fit, nbins=50) :
        '''Read in the fit of the correlation function, xi (from picca_cf.py).
           Create a 2D interpolation object for xi by rt and rp.'''

        h = fitsio.FITS(file_xi)
        print(file_fit)
        ff = h5py.File(file_fit, 'r')
        #base = file_xi 
        base = "lyalya_lyalya"
        fit = ff[base+'/fit'][...] 
        data_rp = h[1]['RP'][:]
        data_rt = h[1]['RT'][:]
        hh = h[1].read_header()
        rpmin = hh['RPMIN']
        rpmax = hh['RPMAX']
        rtmin = 0 
        rtmax = hh['RTMAX']
        h.close()
        ff.close()

        rpmin = data_rp.reshape(50, 50)[0].max()
        rpmax = data_rp.reshape(50, 50)[-1].min()
        rtmin = data_rt.reshape(50, 50)[:, 0].max()
        rtmax = data_rt.reshape(50, 50)[:, -1].min()

        #-- create the regular grid for griddata
        rp = np.linspace(rpmin, rpmax, nbins)
        rt = np.linspace(rtmin, rtmax, nbins)
        xim = sp.interpolate.griddata((data_rt, data_rp), fit, \
                    (rt[:, None], rp[None, :]), method='cubic')

        #-- create interpolator object
        xi2d = sp.interpolate.RectBivariateSpline(rt, rp, xim)

        kappa.xi2d = xi2d
        return xi2d


    @staticmethod
    def read_deltas(filename):
    #def read_deltas(in_dir, nspec=None):
    #    data = {}
    #    ndata = 0
    #    dels = []

    #    fi = glob.glob(in_dir+"/*.fits.gz")
    #    for i,f in enumerate(fi):
    #        print(f"Read {i} of {len(fi)} {ndata}")
    #        hdus = fitsio.FITS(f)
    #        dels += [Delta.from_fitsio(h) for h in hdus[1:]]
    #        ndata+=len(hdus[1:])
    #        hdus.close()
    #        if nspec and ndata > nspec:
    #            break

    #    t1 = time.time()
    #    print(f'Time reading deltas : {(t1-t0)/60:.3f} minutes')

        dels = pkl.load_deltas_from_pickle(filename)
        t1 = time.time()
        print(f'Time reading the pickle file : {(t1-t0)/60:.3f} minutes')
        ndata = len(dels)
        data = {}

        phi = [d.ra for d in dels]
        th = [np.pi/2.-d.dec for d in dels]
        pix = hp.ang2pix(kappa.nside_data, th, phi)

        z_min_pix = 10**dels[0].log_lambda[0]/kappa.lambda_abs-1

        for d,p in zip(dels,pix):
            if p not in data:
                data[p]=[]

            data[p].append(d)
            z = 10**d.log_lambda/kappa.lambda_abs-1.
            z_min_pix = np.amin( np.append([z_min_pix], z) )
            d.z = z
            d.r_comov = kappa.cosmo.get_r_comov(z)
            d.weights *= ((1.+z)/(1.+2.25))**(2.9-1.)
            d.project()

        kappa.z_min_pix = z_min_pix
        kappa.angmax = 2.*\
                np.arcsin(kappa.rt_max/(2.*kappa.cosmo.get_r_comov(z_min_pix)))
        kappa.ndata = ndata
        kappa.data = data

        return data

    @staticmethod
    def fill_neighs():
        '''Find all pairs of deltas (forests) within given separation distance:
        Read in data (dictionary) of deltas (key is the Healpixel). For each
        delta in each Healpixel, get a list of all Healpixels within a given
        radius (angmax) of the delta. Reject any that have no deltas in the
        data. List all neighbouring deltas within selected Healpixels. Find
        angles between delta and its neighbours. Ensure each angle is less than
        angmax. Avoid double counting of pairs (d1.ra > d.ra). Include list of
        neighbours in data dictionary.'''

        data = kappa.data
        print('\n Filling neighbors')
        for healpix in data.keys():
            for delta in data[healpix]:
                healpix_neighbours = hp.query_disc(
                    kappa.nside_data,
                    [delta.x_cart, delta.y_cart, delta.z_cart],
                    kappa.angmax, inclusive = True)
                healpix_neighbours = [
                    other_healpix for other_healpix in healpix_neighbours
                    if other_healpix in data
                ]
                neighbours = [
                    other_delta for other_healpix in healpix_neighbours
                    for other_delta in data[other_healpix]
                ]
                ang = delta.get_angle_between(neighbours)
                w = ang < kappa.angmax
                neighbours = np.array(neighbours)[w]
                delta.neighbours = [
                    other_delta for other_delta in neighbours
                    if delta.ra > other_delta.ra]

    @staticmethod
    def get_kappa(pixels):

        id1 = []
        id2 = []
        skappa = []
        wkappa = []

        for healpix in pixels:
            for i,d1 in enumerate(kappa.data[healpix]):
                kcounter = round(kappa.counter.value*100./kappa.ndata,2)
                if (kcounter%1==0):
                    sys.stderr.write("\rcomputing kappa: {}%".format(kcounter))
                with kappa.lock:
                    kappa.counter.value += 1
                for d2 in d1.neighbours:

                    #-- angle between skewers
                    ang = d1.get_angle_between(d2)
                    if kappa.true_corr:
                        ang_delensed = d1.delensed_angle(d2)

                    if kappa.true_corr:
                        sk, wk = kappa.fast_kappa_true(
                                d1.z, d1.r_comov,
                                d2.z, d2.r_comov,
                                ang_delensed, ang)
                    else:
                        sk, wk = kappa.fast_kappa(
                                d1.z, d1.r_comov, d1.weights, d1.delta,
                                d2.z, d2.r_comov, d2.weights, d2.delta,
                                ang)

                    if wk != 0:
                        id1.append(d1.thingid)
                        id2.append(d2.thingid)
                        skappa.append(sk)
                        wkappa.append(wk)


                setattr(d1, "neighs", None)

        return id1, id2, skappa, wkappa


    @staticmethod
    def fast_kappa(z1,r1,w1,d1,z2,r2,w2,d2,ang):
        '''Estimates kappa using quadratic estimator.'''

        rp = abs(r1-r2[:,None])*np.cos(ang/2)
        rt = (r1+r2[:,None])*np.sin(ang/2)
        d12 = d1*d2[:, None]
        w12 = w1*w2[:,None]

        w = (rp>=kappa.rp_min) & (rp<=kappa.rp_max) & (rt<=kappa.rt_max) & (rt>=kappa.rt_min)

        rp = rp[w]
        rt = rt[w]
        w12 = w12[w]
        d12 = d12[w]

        #-- getting model and first derivative
        xi_model = kappa.xi2d(rt, rp, grid=False)
        xip_model = kappa.xi2d(rt, rp, dx=1, grid=False)

        #-- weight of estimator
        R = -1/(xip_model*rt)

        ska = np.sum( ((d12 - xi_model)*w12 /R))
        wka = np.sum( w12/R**2 )

        return ska, wka

    @staticmethod
    def fast_kappa_true(z1, r1, z2, r2, ang, ang_lens):
        '''Computes kappa using the true correlation function instead of
           estimating it from the product of forest deltas.'''

        rp      = (r1-r2[:,None])*np.cos(ang/2)
        rt      = (r1+r2[:,None])*np.sin(ang/2)
        rp_lens = (r1-r2[:,None])*np.cos(ang_lens/2)
        rt_lens = (r1+r2[:,None])*np.sin(ang_lens/2)

        #z = (z1+z2[:,None])/2

        w = (rp>=kappa.rp_min) & (rp<=kappa.rp_max) & (rt<=kappa.rt_max) & (rt>=kappa.rt_min)

        rp = rp[w]
        rt = rt[w]
        rp_lens = rp_lens[w]
        rt_lens = rt_lens[w]
        #z  = z[w]

        #-- getting model and first derivative
        xi_model  = kappa.xi2d(rt,      rp,       grid=False)
        xi_lens   = kappa.xi2d(rt_lens, rp_lens,  grid=False)
        xip_model = kappa.xi2d(rt,      rp, dx=1, grid=False)
        R = -1/(xip_model*rt)

        ska = np.sum( (xi_lens - xi_model)/R)
        wka = np.sum( 1/R**2  )

        return ska, wka

def compute_kappa(p):
    id1, id2, skappa, wkappa = kappa.get_kappa(p)
    return id1, id2, skappa, wkappa

if __name__=='__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Compute the convergence (kappa) between pairs of deltas.')

    parser.add_argument('--in-dir', required=True, type=str, \
               help='folder containing deltas in pix format')
    parser.add_argument('--xi', required=True, \
               help='fits file containing model')
    parser.add_argument('--fit', required=True, \
               help='h5 file containing fit')
    parser.add_argument('--out', required=True, \
               help='output h5 file with kappa values')
    parser.add_argument('--nproc', required=False, type=int, default=1, \
               help='number of procs used in calculation')
    parser.add_argument('--nspec', required=False, type=int, default=None, \
               help='number of spectra to process')
    parser.add_argument('--rt_min', required=False, type=float, default=3., \
               help='minimum transverse separation')
    parser.add_argument('--rp_min', required=False, type=float, default=3., \
               help='minimum radial separation')
    parser.add_argument('--rt_max', required=False, type=float, default=40., \
               help='maximum transverse separation')
    parser.add_argument('--rp_max', required=False, type=float, default=10., \
               help='maximum radial separation')
    parser.add_argument('--true_corr', required=False, default=False,\
               action='store_true', help='use actual lensed correlation')
    parser.add_argument('--nside', required=False, type=float, default=256, \
               help='resolution of map')
    args, unknown = parser.parse_known_args()

    t0 = time.time()

    kappa.true_corr = args.true_corr
    kappa.rt_min = args.rt_min
    kappa.rp_min = args.rp_min
    kappa.rt_max = args.rt_max
    kappa.rp_max = args.rp_max
    kappa.nside  = args.nside
    kappa.load_model(args.xi, args.fit)
    #kappa.read_deltas(args.deltas)
    kappa.read_deltas(args.in_dir, nspec=args.nspec)
    t1 = time.time()
    kappa.fill_neighs()
    t2 = time.time()
    print(f'kappa.py - Time filling neighbours: {(t2-t1)/60:.3f} minutes')

    cpu_data = {}
    for p in kappa.data.keys():
        cpu_data[p] = [p]
    print(' ', len(kappa.data.keys()), 'pixels with data')

    #-- Memory usage before multiprocessing
    mem_dict = dict(psutil.virtual_memory()._asdict())
    for p in mem_dict:
        print(p, mem_dict[p])
    nprocs_ideal = int(np.floor( mem_dict['available']/mem_dict['used']))+1
    nprocs_avail = cpu_count()
    nprocs = nprocs_avail if nprocs_avail < nprocs_ideal else nprocs_ideal
    if nprocs > 1:
        nprocs -= 1
    print('Will use', nprocs, 'processors')

    pool = Pool(processes=nprocs)
    results = pool.map(compute_kappa, cpu_data.values())
    pool.close()

    t3 = time.time()
    print(f'kappa.py - Time computing kappa: {(t3-t2)/60:.3f} minutes')

    id1 = np.empty(0, dtype=int)
    id2 = np.empty(0, dtype=int)
    skappa = np.empty(0)
    wkappa = np.empty(0)
    
    for r in results:
        id1 = np.append(id1, np.array(r[0]).astype(int))
        id2 = np.append(id2, np.array(r[1]).astype(int))
        skappa = np.append(skappa, np.array(r[2]))
        wkappa = np.append(wkappa, np.array(r[3]))

    #-- Write file with estimated kappa and weights for each pair of forests
    f = h5py.File(args.out, "w")
    thid1 = f.create_dataset("THID1", data=id1)
    thid2 = f.create_dataset("THID2", data=id2)
    kap   = f.create_dataset("SKAPPA", data=skappa)
    wkap  = f.create_dataset("WKAPPA", data=wkappa)

    f.attrs['RPMIN'] = kappa.rp_min
    f.attrs['RPMAX']=kappa.rp_max
    f.attrs['RTMIN']=kappa.rt_min
    f.attrs['RTMAX']=kappa.rt_max
    f.attrs['NT']=kappa.nt
    f.attrs['NP']=kappa.np
    f.attrs['NSIDE']=kappa.nside
    f.close()

    t4 = time.time()
    print(f'kappa.py - Time total : {(t4-t0)/60:.3f} minutes')




