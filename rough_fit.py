from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import numpy as np
from os import listdir
import spectres

def parse_filename(filename):
    teff = ''.join([x for x in filename.split('_')[0] if x.isdigit()])
    logg = ''.join([x for x in filename.split('_')[1] if x.isdigit()])
    return(int(teff)/1000, float(logg)/100.)

def window_crop(wvl, flx, window):
    mindex, maxdex = np.searchsorted(wvl, window)
    wvl = wvl[mindex:maxdex]
    flx = flx[mindex:maxdex]
    return(wvl, flx)

def read_model(filename):
    wvl = []
    flx = []
    with open(filename, 'r') as file:
        for line in file:
            if line[0] != '#':
                splitstr = [float(x) for x in line[:-2].split(' ') if x]
                wvl.append(splitstr[0]/10)
                flx.append(splitstr[1])
    wvl_out = np.asarray(wvl)
    flx_out = np.asarray(flx)
    return(wvl_out, flx_out)

def lazychisq(o, e):
    chisq = []

    for i in range(len(o)):
        chisq.append((o[i]-e[i])**2/e[i])

    return(sum(chisq[10:-10]))

###### CONTROL STUFF ######
wvl_window = (400, 550)
instrument_convolve = False
# Line wavelength
lmbd = 434.047
# Resolution
R = 8000
dlmbd = lmbd/R
###########################

files = [f for f in listdir('koester2/') if '.dk.dat' in f]
files.sort()

# DEFINING CONTINUUM REGIONS

fig, ax = plt.subplots()

cont_as = []
cont_bs = []
wvl_test, flx_test = read_model(f'koester2/{files[0]}')

mindex, maxdex = np.searchsorted(wvl_test, wvl_window)
wvl_test, flx_test = window_crop(wvl_test, flx_test, wvl_window)
ax.plot(wvl_test,flx_test,'k')
ax.set_xlim(wvl_window)
ax.set_title('Select continuum regions')

def select_continuum(eclick, erelease):
    
    x1, x2 = eclick.xdata, erelease.xdata

    cont_as.append(x1)
    cont_bs.append(x2)
    mindex, maxdex = np.searchsorted(wvl_test, (x1, x2))
    ax.plot(wvl_test[mindex:maxdex], flx_test[mindex:maxdex], 'r')
    ax.axvspan(x1, x2, ymin=0, ymax=1, color='y', alpha=0.3)
    plt.draw()

selectors = (RectangleSelector(
        ax, select_continuum,
        useblit=True,
        button=[1, 3],
        minspanx=2, minspany=2,
        spancoords='pixels',
        interactive=True))
plt.show()

# LOADING OBSERVED SPECTRUM

wvl0, flx0 = read_model('ZTFJ0406_UVB_coadd.dat')
wvl0 *= 10
wvl0, flx0 = window_crop(wvl0, flx0, wvl_window)

# CREATE EVENLY-SPACED WVL AXIS AND REBIN

new_axis = np.linspace(min(wvl0), max(wvl0), len(wvl0))
flx0 = spectres.spectres(new_wavs = new_axis,
                        spec_wavs = wvl0,
                        spec_fluxes = flx0,
                        spec_errs = np.zeros(len(flx0)),
                        fill=1.0,
                        verbose=False)[0]

flx0s = convolve(flx0, Box1DKernel(10))

chisqs = []

# LOAD MODELS AND NORMALISE

for i in range(len(files)):
    file = files[i]
    print(f'\rREADING FILE {i+1}/{len(files)} {file}...', end='')
    wvl1, flx1 = read_model(f'koester2/{file}')
    wvl1, flx1 = window_crop(wvl1, flx1, wvl_window)
    cont_wvl = np.asarray([])
    cont_flx = np.asarray([])
    for i in range(len(cont_as)):
        region_wvl, region_flx = window_crop(wvl1, flx1, (cont_as[i], cont_bs[i]))
        cont_wvl = np.concatenate((cont_wvl, region_wvl), axis=None)
        cont_flx = np.concatenate((cont_flx, region_flx), axis=None)
    fit = np.polyfit(cont_wvl, cont_flx, 4)
    p = np.poly1d(fit)
    cont_curve = p(wvl1)
    flx1 /= cont_curve
    # REBIN TO OBSERVATION
    flx2 = spectres.spectres(new_wavs = new_axis,
                            spec_wavs = wvl1,
                            spec_fluxes = flx1,
                            spec_errs = np.zeros(len(flx1)),
                            fill=1.0,
                            verbose=False)[0]
    # INSTRUMENTAL BROADENING
    if instrument_convolve:
        pix_size = new_axis[1]-new_axis[0]
        fwhm = pix_size * dlmbd
        flx2 = convolve(flx2, Gaussian1DKernel(fwhm))

    chisq = lazychisq(flx0, flx2)
    chisqs.append(chisq)

best = np.argmin(chisqs)
print(f'\nLOWEST CHI SQUARED = {np.min(chisqs)}')
print(f'FOR MODEL {files[best]}')

# PLOT BEST FIT

fig2, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios':[2,1]})
fig2.tight_layout()

wvl3, flx3 = read_model(f'koester2/{files[best]}')
wvl3, flx3 = window_crop(wvl3, flx3, wvl_window)
cont_wvl = np.asarray([])
cont_flx = np.asarray([])
for i in range(len(cont_as)):
    region_wvl, region_flx = window_crop(wvl3, flx3, (cont_as[i], cont_bs[i]))
    cont_wvl = np.concatenate((cont_wvl, region_wvl), axis=None)
    cont_flx = np.concatenate((cont_flx, region_flx), axis=None)
fit = np.polyfit(cont_wvl, cont_flx, 4)
p = np.poly1d(fit)
cont_curve = p(wvl3)
flx3 /= cont_curve

ax[0].plot(new_axis, flx0, 'gray')
ax[0].plot(new_axis, flx0s, 'k')
ax[0].plot(wvl3, flx3, 'r')
ax[0].set_xlabel('Wavelength (nm)')
ax[0].set_ylabel('Normalised Flux')

# CHISQ SURFACE PLOT
teffs = []
loggs = []

for file in files:
    teff, logg = parse_filename(file)
    teffs.append(teff)
    loggs.append(logg)

teffs = list(set(teffs))
teffs.sort()
tefflocs = np.arange(0, len(teffs))
loggs = list(set(loggs))
loggs.sort()
logglocs = np.arange(0, len(loggs))

best_teff, best_logg = parse_filename(files[best])
best_teff = teffs.index(best_teff)
best_logg = loggs.index(best_logg)

chisqs_plot = np.asarray(chisqs).reshape(len(teffs),len(loggs))

im = ax[1].imshow(chisqs_plot, aspect='auto')
ax[1].scatter(best_logg, best_teff, color='r', marker='x')
ax[1].invert_yaxis
ax[1].set_yticks(tefflocs, labels=teffs)
ax[1].set_ylabel('$T_eff$ (kK)')
ax[1].set_xticks(logglocs, labels=loggs, rotation=90)
ax[1].set_xlabel('log g')
plt.show()
