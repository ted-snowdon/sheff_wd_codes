#Measures shifts in bright lines in VIS arm 2D sky plots
#DOES NOT APPLY CORRECTION - writes it to header instead
from astropy.io import fits
import numpy as np
from os import listdir
from scipy.constants import speed_of_light
from scipy.optimize import curve_fit
from tkinter.filedialog import askdirectory

def gauss(x, *p):
    #Needed for fitting line positions later
    a, b, c = p
    return a*np.exp(-(x-b)**2/(2.*c**2))

### SKY LINE REFERENCE WAVELENGTHS ###
#Taken from Osterbrock et al. 1996, add/change if desired
ref_wvls = [630.0304, 791.3708, 834.4602, 882.7096]
#Wavelength margin around line position to use (nm)
margin = 0.1

print('SELECT TARGET FOLDER:')
#Observations should be organised as /STAR/VIS
source = askdirectory()
dirlist = [f for f in listdir(source) if f != 'README']

sky_files = []

for dir in dirlist:
    folder = source+'/'+dir
    filenames = [f for f in listdir(folder) if 'SKY_SLIT_MERGE1D_VIS.fits' in f]
    for filename in filenames:
        sky_files.append(folder+'/'+filename)

print(f'FOUND {len(sky_files)} SKY FILES')

for filename in sky_files:
    print('#======================================#')
    print(filename)
    with fits.open(filename, mode='update') as hdul:
        hdr = hdul[0].header
        cdelta = hdr['CDELT1']
        crval1 = hdr['CRVAL1']
        flx = hdul[0].data
        wvl = (np.arange(len(flx))*cdelta)+crval1
        vel_shifts = 0

        for line in ref_wvls:
            a = np.argwhere(wvl >= line-margin)
            b = np.argwhere(wvl <= line+margin)
            lindex = np.intersect1d(a,b)
            line_flx = flx[lindex]
            line_wvl = wvl[lindex]
            #Fit gaussian to find line centre
            p0 = [max(line_flx),line,0.01]
            coeff, var_matrix = curve_fit(gauss,line_wvl,line_flx,p0=p0)
            #Calculate shifts
            peak_pos = np.round(coeff[1], 4)
            wvl_shift = line-peak_pos
            vel_shift = (wvl_shift*speed_of_light)/line/1000.
            vel_shifts += vel_shift
            print(f'PEAK POSITION:\t{peak_pos}')
            print(f'REF. POSITION:\t{line}')
            print(f'DV = {vel_shift} KM/S')

        mean_shift = np.round(vel_shifts/len(ref_wvls), 4)
        print(f'MEAN SHIFT:\t{mean_shift}')
        #Write to header
        hdr['SKYSHIFT'] = (mean_shift, 'sky line velocity along LOS')
    






