## SY 7/6/18
## Create a text file with a table of all unlensed and lensed RA & DEC for delta files

import fitsio
import numpy as N
import sys
import glob

# input directory name (lensed deltas)
dir = sys.argv[1]

# output file name
fout = open(sys.argv[2], 'w')

print('RA, DEC, RA0, DEC0, Z, THINGID', file = fout)

alldeltas =  glob.glob(dir+'/*.fits.gz')
ndel = len(alldeltas)
i=0
for filename in alldeltas:
    hdus = fitsio.FITS(filename)
    print(i, ndel)
    i+=1

    for hdu in hdus[1:]:
        header = hdu.read_header()
        ra = header['RA']
        dec = header['DEC']
        ra0 = header['RA0']
        dec0 = header['DEC0']
        z = header['Z']
        thingid = header['THING_ID']
        print (ra, dec, ra0, dec0, z, thingid, file = fout)

fout.close()
