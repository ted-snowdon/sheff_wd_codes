# Use pre-supplied ephemerides to determine phase of .fits observations

from astropy.io import fits
from os import listdir

# t0 = JD time of observation corresponding to t=0 (epoch time)
# p = period

# ZTFJ0126+1210
#t0 = 2458717.9237281
#p = 0.07087968994

# ZTFJ0406+0958
#t0 = 2458790.986815
#p = 0.1147833978

# ZTFJ0618-0919
#t0 = 2459137.006
#p = 0.292698

# ZTFJ1922+1038
#t0 = 2458206.226915
#p = 0.17245542

# ZTFJ2220+0721
#t0 = 2458676.57946
#p = 0.14041767

files = [f for f in listdir() if '.fits' in f]

for file in files:
    hdul = fits.open(file, mode='update')
    bjd = hdul[0].header['BJD-OBS']
    phi = (bjd-t0)/p % 1
    phi_card = ('PHI', phi, 'Orbital phase')
    hdul[0].header.insert(16, phi_card)
    hdul.close()