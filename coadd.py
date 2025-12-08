from astropy.io import fits
from itertools import islice
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from scipy.constants import speed_of_light
import spectres

ref_lines = [656.279,   #H alpha
             818.326,   # Na doublet
             819.482,
             766.490,   # K doublet
             769.896,
             1252.214   # NIR K line
]

save_flag = True        # If true, write coadded spec to file
line_margin = 10        # Margin in nm around line to plot

c = speed_of_light/1000.

for i in range(len(ref_lines)):
    print(f'[{i}] - {ref_lines[i]} nm')

line_selection = int(input(f'\nEnter index of line to select [0-{len(ref_lines)-1}]: '))
selected_line = ref_lines[line_selection]

lines_file = [f for f in listdir() if '_lines.csv' in f][0]

filenames = []
vrads = []

with open(lines_file, 'r') as file:
    for line in islice(file, 1, None):
        splitstr = line.split(',')
        if float(splitstr[4]) == selected_line:
            filenames.append(splitstr[0])
            vrads.append(float(splitstr[-2]))

wvls = []
flxs = []
errs = []
snrs = []

obj_name = None

if selected_line >= 1000:
    folder = 'NIR/'
else:
    folder = 'VIS/'

for i in range(len(filenames)):
    hdu = fits.open(folder+filenames[i])

    if obj_name == None:
        obj_name = hdu[0].header['OBJECT']

    wvl = hdu[1].data['WAVE'][0]
    flx = hdu[1].data['FLUX'][0]
    err = hdu[1].data['ERR'][0]
    
    wvl *= (1. - vrads[i]/c)
    mindex, maxdex = np.searchsorted(wvl, 
                                    (selected_line-line_margin,
                                    selected_line+line_margin))
    wvl = wvl[mindex:maxdex]
    flx = flx[mindex:maxdex]
    err = err[mindex:maxdex]

    snr = np.mean(flx/err)

    med = np.median(flx)
    sig = np.std(flx)

    wvl_ctm = wvl.copy()
    flx_ctm = flx.copy()
    
    # Take all data within 1 sigma of median as continuum
    for i in range(len(flx_ctm)):
        if not med-sig <= flx_ctm[i] <= med+sig:
            wvl_ctm[i] = None
            flx_ctm[i] = None

    # Fit parabola to continuum
    wvl_ctm = wvl_ctm[~np.isnan(wvl_ctm)]
    flx_ctm = flx_ctm[~np.isnan(flx_ctm)]
    continuum = np.polyfit(wvl_ctm, flx_ctm, 2)

    # Divide flux by continuum parabola to normalise
    for i in range(len(flx)):
        y_curve = (continuum[0]*(wvl[i]**2))+(continuum[1]*wvl[i])+continuum[2]
        flx[i] /= y_curve

    wvls.append(wvl)
    flxs.append(flx)
    errs.append(err)
    snrs.append(snr)

wvl0 = wvls[0]                   # Base wvl array to rebin to
coadd_flux = flxs[0]*snrs[0]

# Rebin spectra to common wavelength axis
for i in range(1, len(wvls)):
    flxs[i], errs[i] = spectres.spectres(new_wavs = wvl0,
                                spec_wavs = wvls[i],
                                spec_fluxes = flxs[i],
                                spec_errs = errs[i],
                                fill=None,
                                verbose=False)
    coadd_flux += (flxs[i]*snrs[i])
    plt.plot(wvl0, flxs[i], alpha=0.3)

coadd_flux /= np.sum(snrs)

plt.plot(wvl0, coadd_flux, color='k')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Normalised Flux')
plt.title(f'{obj_name} / {str(selected_line)} nm')
plt.show()
