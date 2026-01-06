# Takes spectrum in .dat format and computes normalisation curve
# Continuum regions defined by user

import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import numpy as np
from tkinter.filedialog import askopenfilename

fig, ax = plt.subplots()

def read_model(teff, logg):
    filename = askopenfilename()
    wvl = []
    flx = []
    nmify = input('Convert flux to nm? (y/n): ')
    with open(filename, 'r') as file:
        for line in file:
            if line[0] != '#':
                splitstr = [float(x) for x in line[:-2].split(' ') if x]
                if nmify == 'y':
                    wvl.append(splitstr[0]/10)
                else:
                    wvl.append(splitstr[0])
                flx.append(splitstr[1])
    return(wvl, flx)

cont_as = []
cont_bs = []
wvl, flx = read_model(20000, 800)
wvl_window = (400, 550)
mindex, maxdex = np.searchsorted(wvl, wvl_window)
wvl = wvl[mindex:maxdex]
flx = flx[mindex:maxdex]
ax.plot(wvl,flx,'k')
ax.set_xlim(wvl_window)
ax.set_title('Select continuum regions')

def select_continuum(eclick, erelease):
    
    x1, x2 = eclick.xdata, erelease.xdata

    cont_as.append(x1)
    cont_bs.append(x2)
    mindex, maxdex = np.searchsorted(wvl, (x1, x2))
    ax.plot(wvl[mindex:maxdex], flx[mindex:maxdex], 'r')
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

cont_wvl = np.asarray([])
cont_flx = np.asarray([])

if len(cont_as) == 0:
    print('No continuum regions selected \n Closing...')
    quit()

for i in range(len(cont_as)):
    contdex = np.searchsorted(wvl, (cont_as[i], cont_bs[i]))
    region_wvl = wvl[contdex[0]:contdex[1]]
    region_flx = flx[contdex[0]:contdex[1]]
    plt.plot(region_wvl, region_flx, 'r')
    cont_wvl = np.concatenate((cont_wvl, region_wvl), axis=None)
    cont_flx = np.concatenate((cont_flx, region_flx), axis=None)

fit = np.polyfit(cont_wvl, cont_flx, 4)
p = np.poly1d(fit)
cont_curve = p(wvl)

plt.plot(wvl, flx, 'k')
plt.plot(wvl, cont_curve, 'b')
plt.title('Normalisation Curve')
plt.show()

norm_flx = flx / cont_curve

plt.plot(wvl, norm_flx, 'k')
plt.title('Normalised Spectrum')
plt.show()

save_yn = input('Save to file? (y/n): ')

if save_yn == 'y':
    savename = input('Enter filename to save as: ')
    with open(savename, 'w') as file:
        for i in range(len(wvl)):
            writestr = ' '.join([str(wvl[i]), str(norm_flx[i]), '\n'])
            file.write(writestr)
