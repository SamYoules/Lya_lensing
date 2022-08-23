#-- SY 6/5/18
#-- Plots displacements of Mock QSOs when lensed by a giant blob placed in centre of NGC quasars

import numpy as np
from matplotlib import pylab as plt
from astropy.table import Table
import sys
#-- For removing duplicate labels in plot
from collections import OrderedDict

#-- input pathname (e.g. outputs/table-blob.txt)
#fin1 = sys.argv[1]
fin1 = "Misc/table-blob.txt"

#-- Create a table of unlensed and lensed mocks qsos
t = Table(np.loadtxt(fin1, skiprows=2), names=('ra', 'dec', 'ra0', 'dec0', 'z', 'thing_id'))

#-- Pick random numbers to use as indices for plotting a thinned out sample of the data
rand_ind = np.random.choice(np.arange(len(t)), replace=False, size=1000)

plt.rcParams.update({'font.size':14})

plt.figure(figsize=(7, 7))
for i in rand_ind:
    x = t['ra0'][i]
    y = t['dec0'][i]
    dx = t['ra'][i] - t['ra0'][i]
    dy = t['dec'][i] - t['dec0'][i]
    plt.arrow(x, y, dx, dy, color='#141651', width=0.001)
plt.xlim((2.6,3.6))
plt.ylim((0.0,1.0))
plt.xlabel('RA [Radians]')
plt.ylabel('DEC [Radians]')
plt.title('Displacements in Mock QSOs')
plt.savefig('Plots/displaced_qsos.png')
plt.savefig('Plots/displaced_qsos.pdf')

