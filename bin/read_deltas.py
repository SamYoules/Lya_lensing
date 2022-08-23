import fitsio, glob, pickle, time
import numpy as np
from picca import io_lens, data_lens, constants
from picca.data_lens import Delta
import sys

delta_path = sys.argv[1]
filename   = sys.argv[2]

#delta_path = '/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Deltas/Blob/'
#filename = '/global/cfs/projectdirs/desi/users/syoules/Lensing/2.0-1000/Deltas/blob-deltas.pkl'

all_pix = glob.glob(delta_path+'/*.fits.gz')
print(len(all_pix))

def read_deltas(delta_path, max_num_spec=None):
    ''' This is just a wrapper of picca.io.read_deltas
        It reads fits files, and computes comoving distances
        and updates weights. Also does the projection.
    '''
    # load fiducial cosmology
    cosmo = constants.Cosmo(Om=0.31,
                            Or=0.,
                            Ok=0.,
                            wl=-1.)

    t0 = time.time()
    data, num_data, z_min, z_max = io_lens.read_deltas(delta_path,
                    16,
                    1216.,
                    1.,
                    2.25,
                    cosmo,
                    max_num_spec=max_num_spec,
                    no_project=True)
    t1 = time.time()
    print(f'Time reading fits files : {(t1-t0)/60:.3f} minutes')
    return data

def read_deltas_simple(delta_path, max_num_spec=None):
    ''' Simply reads the fits files into a list. 
        It does not compute distances as read_deltas does.
    '''
    t0 = time.time()
    files = np.sort(glob.glob(delta_path+'/*.fits.gz'))
    deltas = []
    num_data = 0
    index = 0
    for filename in files:
        print(f"Read {index} of {len(files)} {num_data}")
        hdul = fitsio.FITS(filename)
        deltas += [Delta.from_fitsio(hdu) for hdu in hdul[1:]]
        index += 1
        hdul.close()
        num_data = len(deltas)
        if max_num_spec is not None and num_data > max_num_spec:
            break
    t1 = time.time()
    print(f'Time reading fits files : {(t1-t0)/60:.3f} minutes')
    
    return deltas

def save_deltas_to_pickle(deltas, filename):
    pickle.dump(deltas, open(filename, 'wb'))

def load_deltas_from_pickle(filename):
    t0 = time.time()
    deltas = pickle.load(open(filename, 'rb'))
    t1 = time.time()
    print(f'Time reading the pickle file : {(t1-t0)/60:.3f} minutes')
    
    return deltas

#deltas = read_deltas(delta_path, max_num_spec=10000)
#save_deltas_to_pickle(deltas, 'deltas.pkl')
#deltas_pickle = load_deltas_from_pickle('deltas.pkl')

#deltas_simple = read_deltas_simple(delta_path, max_num_spec=10000)
#save_deltas_to_pickle(deltas_simple, 'deltas_simple.pkl')
#deltas_simple_pickle = load_deltas_from_pickle('deltas_simple.pkl')

#len(deltas_simple), len(deltas_simple_pickle)

#for k in deltas.keys():
#    print(f'Pixel #{k} {len(deltas[k])} {len(deltas_pickle[k])}')

#-- Pick a few random deltas and compare
#keys = list(deltas.keys())
#for k in keys[::3]:
#    d0 = deltas[k][0]
#    d1 = deltas_pickle[k][0]
#    plt.figure()
#    plt.plot(d0.delta)
#    plt.plot(d1.delta)

