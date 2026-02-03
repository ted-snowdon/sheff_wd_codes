from astropy.io import fits
from itertools import islice
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import numpy as np
from os import listdir
from scipy.constants import speed_of_light
import spectres

ref_lines = [656.279,   #H alpha
             818.326,   # Na doublet
             819.482,
             766.490,   # K doublet
             769.896,
             1252.214,  # NIR K line
             'NIR',
             'VIS',
             'UVB'
]

save_flag = True        # If true, write coadded spec to file
p0 = [-217.0765173,      # Guess params for sine curve fit
      2*np.pi, 
      3, 
      45.00533485]
line_margin = 10        # Margin in nm around line to plot

def select_continuum(eclick, erelease):
    
    x1, x2 = eclick.xdata, erelease.xdata

    cont_as.append(x1)
    cont_bs.append(x2)
    mindex, maxdex = np.searchsorted(wvl, (x1, x2))
    ax.plot(wvl[mindex:maxdex], flx[mindex:maxdex], 'r')
    ax.axvspan(x1, x2, ymin=0, ymax=1, color='y', alpha=0.3)
    plt.draw()

def model_sin(x, A, w, p, off):
    return(A * np.sin(w * x + p) + off)

c = speed_of_light/1000.

for i in range(len(ref_lines)):
    print(f'[{i}] - {ref_lines[i]}')

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

if type(selected_line) == str:
    folder = selected_line+'/'
elif selected_line >= 1000:
    folder = 'NIR/'
else:
    folder = 'VIS/'

if len(filenames) == 0:
    filenames = [f for f in listdir(folder) if '.fits' in f]
    filenames.sort()

for i in range(len(filenames)):
    hdu = fits.open(folder+filenames[i])

    if len(vrads) == 0:
        phi = hdu[0].header['PHI']
        vrad = model_sin(phi, p0[0], p0[1], p0[2], p0[3])
    else:
        vrad = vrads[i]

    if obj_name == None:
        obj_name = hdu[0].header['OBJECT']

    wvl = hdu[1].data['WAVE'][0]
    flx = hdu[1].data['FLUX'][0]
    err = hdu[1].data['ERR'][0]
    
    wvl *= (1. - vrad/c)

    if type(selected_line) == float:
        mindex, maxdex = np.searchsorted(wvl, 
                                        (selected_line-line_margin,
                                        selected_line+line_margin))
        wvl = wvl[mindex:maxdex]
        flx = flx[mindex:maxdex]
        err = err[mindex:maxdex]

    snr = np.mean(flx/err)

    med = np.median(flx)
    sig = np.std(flx)

    # Select continuum and normalise
    # Uses same continuum regions for all spectra
    if i == 0: # Only want to do this once
        cont_as = []
        cont_bs = []
        fig, ax = plt.subplots()
        plt.plot(wvl, flx, 'k')
        plt.title('Select continuum regions')
        plt.xlim(min(wvl), max(wvl))
        plt.ylim(-1.1*max(flx), 1.1*max(flx))
        selectors = (RectangleSelector(
                        ax, select_continuum,
                        useblit=True,
                        button=[1, 3],
                        minspanx=2, minspany=2,
                        spancoords='pixels',
                        interactive=True))
        plt.show()

    wvl_ctm = np.asarray([])
    flx_ctm = np.asarray([])

    if len(cont_as) == 0:
        print('No continuum regions selected \n Closing...')
        quit()

    for j in range(len(cont_as)):
        contdex = np.searchsorted(wvl, (cont_as[j], cont_bs[j]))
        region_wvl = wvl[contdex[0]:contdex[1]]
        region_flx = flx[contdex[0]:contdex[1]]
        wvl_ctm = np.concatenate((wvl_ctm, region_wvl), axis=None)
        flx_ctm = np.concatenate((flx_ctm, region_flx), axis=None)

    continuum = np.polyfit(wvl_ctm, flx_ctm, 4)
    p = np.poly1d(continuum)
    y_curve = p(wvl)

    if i == 0:
        plt.plot(wvl, flx, 'k')
        plt.plot(wvl, y_curve, 'b')
        plt.title('Normalisation Curve')
        plt.show()

    flx /= y_curve
    err /= y_curve

    if i == 0:
        plt.plot(wvl, flx, 'k')
        plt.show()

    wvls.append(wvl)
    flxs.append(flx)
    errs.append(err)
    snrs.append(snr)

wvl0 = wvls[0]                   # Base wvl array to rebin to
coadd_flux = flxs[0]*snrs[0]
coadd_errs = errs[0]*snrs[0]

# Rebin spectra to common wavelength axis
for i in range(1, len(wvls)):
    flxs[i], errs[i] = spectres.spectres(new_wavs = wvl0,
                                spec_wavs = wvls[i],
                                spec_fluxes = flxs[i],
                                spec_errs = errs[i],
                                fill=None,
                                verbose=False)
    coadd_flux += (flxs[i]*snrs[i])
    coadd_errs += (errs[i]*snrs[i])
    plt.plot(wvl0, flxs[i], alpha=0.3)

coadd_flux /= np.sum(snrs)
coadd_errs /= np.sum(snrs)

plt.plot(wvl0, coadd_flux, color='k')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Normalised Flux')
plt.title(f'{obj_name} / {str(selected_line)} nm')
plt.show()

if save_flag:
    with open(f'{obj_name}_coadded.dat', 'w') as file:
        for i in range(1, len(wvl0)-1):
            writestr = ' '.join((str(wvl0[i]),str(coadd_flux[i]),str(coadd_errs[i]),'\n'))
            file.write(writestr)
