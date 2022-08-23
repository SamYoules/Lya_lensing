#!/usr/bin/env python

"""Compute estimated kappa for midpoint between qsos using estimator for
   lya x qso, and write to file"""

import numpy as np
import scipy as sp
import fitsio, glob, pickle, time
from astropy.io import fits
import argparse
import healpy as hp
import sys
from scipy import interpolate
import copy
from picca import constants, xcf, io_lens
from picca.data_lens import Delta
from multiprocessing import Pool,Process,Lock,Manager,cpu_count,Value
import h5py
import psutil

class kappa:

    #nside = 256
    nside_data = 32
    rot = hp.Rotator(coord=['C', 'G'])
    lambda_abs = 1215.67

    fid_Om=0.31
    #cosmo = constants.Cosmo(Om=args.fid_Om,
    #                        Or=args.fid_Or,
    #                        Ok=args.fid_Ok,
    #                        wl=args.fid_wl,
    #                        blinding=blinding)
    #cosmo = constants.Cosmo(fid_Om)
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

    counter = Value('i',0)
    lock = Lock()

    data={}
    ndata=0


    @staticmethod
    def load_model(file_xi, file_fit, nbins=50) :

        h = fitsio.FITS(file_xi)
        ff = h5py.File(file_fit, 'r')
        base = "lyalya_qso" 
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

        rpmin = data_rp.reshape(100, 50)[0].max()
        rpmax = data_rp.reshape(100, 50)[-1].min()
        rtmin = data_rt.reshape(100, 50)[:, 0].max()
        rtmax = data_rt.reshape(100, 50)[:, -1].min()

        #-- create the regular grid for griddata
        rp = np.linspace(rpmin, rpmax, nbins*2)
        rt = np.linspace(rtmin, rtmax, nbins)

        #xim = sp.interpolate.griddata((data_rt, data_rp), fit,
        #            (rt[:, None], rp[None, :]), method='cubic')

        xim = sp.interpolate.griddata((data_rt, data_rp), fit,
                    (np.outer(np.ones(rp.size), rt).ravel(),
                    np.outer(rp, np.ones(rt.size)).ravel()), method='cubic')

        #-- create interpolator object
        xi2d = sp.interpolate.RectBivariateSpline(rt, rp,
                   xim.reshape((100, 50)).T )

        kappa.xi2d = xi2d
        return xi2d

    @staticmethod
    #def read_deltas(filename):
    def read_deltas(in_dir, nspec=None):
        data = {}
        ndata = 0
        dels = []

        fi = glob.glob(in_dir+"/*.fits.gz")
        for i,f in enumerate(fi):
            print(f"Read {i} of {len(fi)} {ndata}")
            hdus = fitsio.FITS(f)
            dels += [Delta.from_fitsio(h) for h in hdus[1:]]
            ndata+=len(hdus[1:])
            hdus.close()
            if nspec and ndata > nspec:
                break

        t2 = time.time()
        print(f'Time reading deltas : {(t2-t1)/60:.3f} minutes')

        #dels = pickle.load(open(filename, 'rb'))
        #t2 = time.time()
        #print(f'Time reading the pickle file : {(t2-t1)/60:.3f} minutes')
        #ndata = len(dels)
        #data = {}

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
            #d.r_comov = kappa.cosmo.r_comov(z)
            d.r_comov = cosmo.get_r_comov(z)
            d.weights *= ((1.+z)/(1.+2.25))**(2.9-1.)
            d.project()

        kappa.z_min_pix = z_min_pix
        #r_comov = kappa.cosmo.r_comov(z_min_pix)
        r_comov = cosmo.get_r_comov(z_min_pix)
        kappa.angmax = 2.* np.arcsin(kappa.rt_max/(2.*r_comov))
        kappa.ndata = ndata
        kappa.data = data

        return data

    @staticmethod
    def fill_neighs():
        data = kappa.data
        objs = kappa.objs
        print('\n Filling neighbors')
        for healpix in data.keys():
            for delta in data[healpix]:
                healpix_neighbours = hp.query_disc(
                    kappa.nside_data,
                    [delta.x_cart, delta.y_cart, delta.z_cart],
                    kappa.angmax, inclusive = True)
                healpix_neighbours = [
                    other_healpix for other_healpix in healpix_neighbours
                    if other_healpix in objs
                ]
                neighbours = [
                    quasar for other_healpix in healpix_neighbours
                    for quasar in objs[other_healpix]
                    if quasar.thingid != delta.thingid]
                ang = delta.get_angle_between(neighbours)
                w = ang < kappa.angmax
                neighbours = np.array(neighbours)[w]
                delta.qneighbours = [quasar for quasar in neighbours]

        
    @staticmethod
    def get_kappa(pixels):

        id1 = []
        id2 = []
        skappa = [] 
        wkappa = []

        for healpix in pixels:
            for i,d in enumerate(kappa.data[healpix]):
                sys.stderr.write("\rcomputing kappa: {}%".format(
                       round(kappa.counter.value*100./kappa.ndata,2)))
                with kappa.lock:
                    kappa.counter.value += 1
                for q in d.qneighbours:

                    #-- angle between skewers
                    ang = d.get_angle_between(q)
                    if kappa.true_corr:
                        ang_delensed = d.delensed_angle(q)
 
                    if kappa.true_corr:
                        sk, wk = kappa.fast_kappa_true(
                                d.z, d.r_comov,
                                q.z_qso, q.r_comov,
                                ang_delensed, ang)
                    else:
                        #print(q.z_qso, q.r_comov, q.weights, ang)
                        sk, wk = kappa.fast_kappa(
                                d.z, d.r_comov, d.weights, d.delta,
                                q.z_qso, q.r_comov, q.weights, ang)
                    #print(wk)
                    if wk != 0:
                        id1.append(q.thingid)
                        id2.append(d.thingid)
                        skappa.append(sk)
                        wkappa.append(wk)

                setattr(d, "neighs", None)
                
        return id1, id2, skappa, wkappa

    @staticmethod
    def fast_kappa(z1,r1,w1,d1,zq,rq,wq,ang):
        rp = (r1[:,None]-rq)*np.cos(ang/2)
        rt = (r1[:,None]+rq)*np.sin(ang/2)
        de = d1[:,None]
        we = w1[:,None]*wq

        #print(np.min(rp), np.min(rt)) ## problem numbers too big ~2000, ~90
        w = (rp>=kappa.rp_min) & (rp<=kappa.rp_max) & (rt<=kappa.rt_max) & (rt>=kappa.rt_min)

        #print(np.sum(w)) ## problem. = 0
        rp = rp[w]
        rt = rt[w]
        we = we[w]
        de = de[w]
 
        #-- getting model and first derivative
        xi_model = kappa.xi2d(rt, rp, grid=False)
        xip_model = kappa.xi2d(rt, rp, dx=1, grid=False)

        #-- weight of estimator
        R = -1/(xip_model*rt)

        ska = np.sum( ((de - xi_model)/R*we))
        wka = np.sum( we/R**2 ) 

        return ska, wka

    @staticmethod
    def fast_kappa_true(z1, r1, zq, rq, ang, ang_lens):
        
        rp      = (r1[:,None]-rq)*np.cos(ang/2)
        rt      = (r1[:,None]+rq)*np.sin(ang/2)
        rp_lens = (r1[:,None]-rq)*np.cos(ang_lens/2)
        rt_lens = (r1[:,None]+rq)*np.sin(ang_lens/2)

        w = (rp>=kappa.rp_min) & (rp<=kappa.rp_max) & (rt<=kappa.rt_max) & (rt>=kappa.rt_min)

        rp = rp[w]
        rt = rt[w]
        rp_lens = rp_lens[w]
        rt_lens = rt_lens[w]

        #-- getting model and first derivative
        xi_model  = kappa.xi2d(rt,      rp,       grid=False)
        xi_lens   = kappa.xi2d(rt_lens, rp_lens,  grid=False)
        xip_model = kappa.xi2d(rt,      rp, dx=1, grid=False)
        R = -1/(xip_model*rt)

        ska = np.sum( ((xi_lens - xi_model)/R))
        wka = np.sum( 1/R**2 )

        return ska, wka


def compute_kappa(p):
    id1, id2, skappa, wkappa = kappa.get_kappa(p)
    #print(f"size skappa = {len(skappa)}")
    return id1, id2, skappa, wkappa
            

if __name__=='__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Compute the convergence (kappa) between qsos and deltas.')

    parser.add_argument('--in-dir', required=True, type=str,
               help='pathname to deltas folder')
    parser.add_argument('--xi', required=True,
               help='text file containing model')
    parser.add_argument('--fit', required=True,
               help='text file containing fit')
    parser.add_argument('--out', required=True,
               help='output H5PY file with kappa map')
    parser.add_argument('--nproc', required=False, type=int, default=1,
               help='number of procs used in calculation')
    parser.add_argument('--nspec', required=False, type=int, default=None,
               help='number of spectra to process')
    parser.add_argument('--drq', type=str, default=None,
               required=True, help='Catalog of objects in DRQ format')
    parser.add_argument('--z_min_obj', type=float, default=None,
               required=False, help='Min redshift for object field')
    parser.add_argument('--z_max_obj', type=float, default=None,
               required=False, help='Max redshift for object field')
    parser.add_argument('--z_cut_min', type=float, default=0.,
               required=False,
               help='Use only pairs of forest x object with the mean of the last absorber redshift and the object redshift larger than z_cut_min')
    parser.add_argument('--z_cut_max', type=float, default=10., required=False,
              help='Use only pairs of forest x object with the mean of the last absorber redshift and the object redshift smaller than z_cut_max')
    parser.add_argument('--z_evol_obj', type=float, default=1.,
               required=False, help='Exponent of the redshift evolution of  the object field')
    parser.add_argument('--z_ref', type=float, default=2.25, required=False,
               help='Reference redshift')
    parser.add_argument(
        '--fid-Om',
        type=float,
        default=0.315,
        required=False,
        help='Omega_matter(z=0) of fiducial LambdaCDM cosmology')

    parser.add_argument(
        '--fid-Or',
        type=float,
        default=0.,
        required=False,
        help='Omega_radiation(z=0) of fiducial LambdaCDM cosmology')

    parser.add_argument('--fid-Ok',
                        type=float,
                        default=0.,
                        required=False,
                        help='Omega_k(z=0) of fiducial LambdaCDM cosmology')

    parser.add_argument(
        '--fid-wl',
        type=float,
        default=-1.,
        required=False,
        help='Equation of state of dark energy of fiducial LambdaCDM cosmology')

    parser.add_argument('--rt_min', required=False, type=float, default=3.,
               help='minimum transverse separation')
    parser.add_argument('--rp_min', required=False, type=float, default=-70.,
               help='minimum radial separation')
    parser.add_argument('--rt_max', required=False, type=float, default=70.,
               help='maximum transverse separation')
    parser.add_argument('--rp_max', required=False, type=float, default=70.,
               help='maximum radial separation')
    parser.add_argument('--true_corr', required=False, default=False,
               action='store_true', help='use actual lensed correlation')
    args, unknown = parser.parse_known_args()

    # read blinding keyword
    blinding = io_lens.read_blinding(args.in_dir)

    ### Read objects
    cosmo = constants.Cosmo(Om=args.fid_Om,
                            Or=args.fid_Or,
                            Ok=args.fid_Ok,
                            wl=args.fid_wl,
                            blinding=blinding)

    t0 = time.time()

    # Find the redshift range
    z_min = 1.7
    z_max = 5.0
    if args.z_min_obj is None:
        r_comov_min = cosmo.get_r_comov(z_min)
        r_comov_min = max(0., r_comov_min)
        args.z_min_obj = cosmo.distance_to_redshift(r_comov_min)
        print("z_min_obj = {}".format(args.z_min_obj), end="")
    if args.z_max_obj is None:
        r_comov_max = cosmo.get_r_comov(z_max)
        r_comov_max = max(0., r_comov_max)
        args.z_max_obj = cosmo.distance_to_redshift(r_comov_max)
        print("z_max_obj = {}".format(args.z_max_obj), end="")

    objs,zmin_obj = io_lens.read_objects(args.drq, kappa.nside_data,
                                   args.z_min_obj, args.z_max_obj,
                                   args.z_evol_obj, args.z_ref,cosmo)

    t1 = time.time()
    print(f'picca_kappa.py - Time reading qsos: {(t1-t0)/60:.3f} minutes')

    kappa.objs = objs
    kappa.true_corr = args.true_corr
    kappa.rt_min = args.rt_min
    kappa.rp_min = args.rp_min
    kappa.rt_max = args.rt_max
    kappa.rp_max = args.rp_max
    kappa.z_cut_max = args.z_cut_max
    kappa.z_cut_min = args.z_cut_min
    kappa.load_model(args.xi, args.fit)
    #kappa.read_deltas(args.deltas)
    kappa.read_deltas(args.in_dir, nspec=args.nspec)
    t2 = time.time()
    kappa.fill_neighs()
    t3 = time.time()
    print(f'picca_xkappa.py - Time filling neighbours : {(t3-t2)/60:.3f} minutes')

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

    pool = Pool(processes=args.nproc)
    results = pool.map(compute_kappa, cpu_data.values())
    pool.close()

    t4 = time.time()
    print(f'picca_kappa.py - Time computing kappa: {(t4-t3)/60:.3f} minutes')

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
    thidq = f.create_dataset("THIDQ", data=id1)
    thid  = f.create_dataset("THID", data=id2)
    kap   = f.create_dataset("SKAPPA", data=skappa)
    wkap  = f.create_dataset("WKAPPA", data=wkappa)

    f.attrs['RPMIN'] = kappa.rp_min
    f.attrs['RPMAX']=kappa.rp_max
    f.attrs['RTMIN']=kappa.rt_min
    f.attrs['RTMAX']=kappa.rt_max
    f.attrs['NT']=kappa.nt
    f.attrs['NP']=kappa.np
    f.close()

    t5 = time.time()
    print(f'picca_xkappa.py - Time total : {(t5-t0)/60:.3f} minutes')


