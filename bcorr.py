# loads .fits files in the same directory and applies a barycentric correction to wavelength arrays
# also applies barycentric time correction
# intended to be used for XSHOOTER data following telluric correction

from astropy import time, coordinates as coord, units as u
from astropy.io import fits
import os
from scipy.constants import speed_of_light
from time import sleep

def create_filename(input_filename):
    # molecfit outputs pretty generic filenames
    # creates individual filename from header info
    hdr = fits.getheader(input_filename, ext=0)
    object = hdr['OBJECT'][:10]
    date = hdr['DATE-OBS'][2:10].replace('-','')
    time = hdr['DATE-OBS'][11:19].replace(':','')
    return('_'.join([object,date,time])+'_tb.fits')

def time_correction(jd, skycoords, site):
    # calculate barycentric time correction
    # if not using data from Paranal, change site on line 46
    ra = str(skycoords[0])
    dec = str(skycoords[1])
    ip_peg = coord.SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
    site_coords = coord.EarthLocation.of_site(site)
    times = time.Time(jd, format='jd', scale='utc', location=site_coords)

    return(times.light_travel_time(ip_peg))

c = speed_of_light/1000.

# assuming we're in a directory with folders for each tell. corr. with 1 file each
folders = [f for f in os.listdir() if os.path.isdir(f)]
folders.sort()

for i in range(len(folders)):
    folder = folders[i]
    file = os.listdir(folder)[0]
    path = folder+'/'+file
    os.system(f'cp {path} .')
    hdul = fits.open(file, mode='update')
    # barycentric time correction
    skycoords = [hdul[0].header['RA'], hdul[0].header['DEC']]
    jd1 = hdul[0].header['MJD-OBS'] + 2400000.5
    jd2 = hdul[0].header['MJD-END'] + 2400000.5
    bary_time1 = time_correction(jd1, skycoords, site='Paranal Observatory (ESO)')
    bary_time_card1 = ('BJD-OBS', float(jd1)+bary_time1.value, '[d] Start of observations (days)')
    bary_time2 = time_correction(jd2, skycoords, site='Paranal Observatory (ESO)')
    bary_time_card2 = ('BJD-END', float(jd2)+bary_time2.value, '[d] End of observations (days)')
    hdul[0].header.insert(14, bary_time_card2)
    hdul[0].header.insert(14, bary_time_card1)
    # barycentric velocity correction
    v_b = hdul[0].header['ESO QC VRAD BARYCOR']
    print(f'BARYCOR VELOCITY = {v_b}')
    data = hdul[1].data
    data['WAVE'][0] *= (1. + v_b/c)
    new_name = create_filename(file)
    print(f'CORRECTION APPLIED, SAVED TO {new_name}\n')
    hdul.close()
    os.system(f'mv {file} {new_name}')
    sleep(0.1) # needed to ensure we actually get all the files we want