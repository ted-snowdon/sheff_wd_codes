from astropy.io import fits
from math import ceil
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from scipy.constants import speed_of_light

c = speed_of_light/1000.

file_names = [f for f in listdir('VIS') if '.fits' in f]
file_names.sort()

data = []

for file in file_names:
    hdu = fits.getdata('VIS/'+file, ext=1)[0]
    wvl = hdu['WAVE']
    x = c * (wvl/818.904 - 1.)
    y = hdu['FLUX']
    mindex, maxdex = np.searchsorted(x, (-500,500))
    xx, yy = x[mindex:maxdex], y[mindex:maxdex]
    yy_norm = yy / np.median(yy)
    data.append(yy_norm)

axlab = ['-500', '0', '+500']
axlocs = [0, ceil(len(data[0])/2), len(data[2])-1]

fig, ax = plt.subplots()
im = ax.imshow(data, aspect=2)
ax.invert_yaxis()
ax.set_xticks(axlocs, labels=axlab)
ax.set_xlabel('Velocity / (km/s)')
ax.set_ylabel('Spectrum no.')
ax.set_title('ZTFJ0618 Na 818/819')
cbar = fig.colorbar(im)
plt.show()