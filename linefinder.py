# fit Gaussians to locations of pre-determined spectral features
# determines radial velocity shifts relative to reference wavelength
# intended to be run in the corrected observations folder with UVB/VIS/NIR subfolders

from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel
import csv
from itertools import islice
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from scipy.constants import speed_of_light
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

def gauss(x, y0, a, x0, sigma):
    return y0+a*np.exp(-(x-x0)**2/(2*sigma**2))

def model_sin(x, A, w, p, off):
    return(A * np.sin(w * x + p) + off)

box_kernel = Box1DKernel(5)

c = speed_of_light/1000.

### REFERENCE WAVELENGTHS ###
H_a = 656.279

vis_lines = [818.326,   # Na doublet
             819.482,
             766.490,   # K doublet
             769.896,]

nir_lines = [1252.214   # K line
            ]

Ha_margin = 0.3
vis_margin = 0.4            # wavelength margin around line pos. to use (nm)
nir_margin = 0.5

Ha_plot_flag = False
vis_plot_flag = False
nir_plot_flag = False

#uvb_files = [f for f in listdir('UVB/') if '.fits' in f]
vis_files = [f for f in listdir('VIS/') if '.fits' in f]
nir_files = [f for f in listdir('NIR/') if '.fits' in f]
#uvb_files.sort()
vis_files.sort()
nir_files.sort()

obj_name = (fits.getheader('VIS/'+vis_files[0], ext=0)['OBJECT']).split('.')[0]

with open(f'{obj_name}_lines.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',')
    csvwriter.writerow(['Filename', 'Date', 'Time', 'Phase', 'Line', 'Peak', 'Vel. shift'])
    for file in vis_files:
        hdu = fits.open('VIS/'+file)
        date, time = (hdu[0].header['DATE-OBS']).split('T')
        phi = hdu[0].header['PHI']
        wvl = hdu[1].data['WAVE'][0]
        flx = hdu[1].data['FLUX'][0]
        flx_sm = convolve(flx, box_kernel)

        H_min, H_max = np.searchsorted(wvl, (H_a-3, H_a+3))
        x_Ha, y_Ha = wvl[H_min:H_max], flx_sm[H_min:H_max]
        Ha_peaks = find_peaks(y_Ha)[0]
        Ha_peak_flxs = y_Ha[Ha_peaks]
        Ha_peak_dex = Ha_peaks[np.argmax(Ha_peak_flxs)]
        x_Ha_pk, y_Ha_pk = x_Ha[Ha_peak_dex], y_Ha[Ha_peak_dex]
        Ha_min, Ha_max = np.searchsorted(wvl, (x_Ha_pk-Ha_margin, x_Ha_pk+Ha_margin))
        xx_Ha, yy_Ha = wvl[Ha_min:Ha_max], flx_sm[Ha_min:Ha_max]
        Ha_mean = sum(xx_Ha * yy_Ha)/sum(yy_Ha)
        Ha_sigma = np.sqrt(sum(yy_Ha*(xx_Ha-Ha_mean)**2)/sum(yy_Ha))
        p0_Ha = [np.median(y_Ha), y_Ha_pk, Ha_mean, Ha_sigma]
        coeff, var_matrix = curve_fit(gauss, xx_Ha, yy_Ha, p0=p0_Ha)
        curve_Ha = coeff[0] + coeff[1]*np.exp(-(xx_Ha-coeff[2])**2/(2.*coeff[3]**2))
        Ha_peak_pos = np.round(coeff[2], 4)
        Ha_dwvl = Ha_peak_pos-H_a
        Ha_dvrd = (Ha_dwvl*c)/H_a

        csvwriter.writerow([file, date, time, phi, H_a, Ha_peak_pos, Ha_dvrd])
        if Ha_plot_flag:
            plt.plot(x_Ha, flx[H_min:H_max], color='gray')
            plt.plot(xx_Ha, yy_Ha, color='k')
            plt.plot(xx_Ha, curve_Ha, color='r')
            plt.plot([H_a, H_a],[min(y_Ha), max(y_Ha)], 'g--')
            plt.plot([min(xx_Ha), max(xx_Ha)], [np.median(y_Ha), np.median(y_Ha)], 'r|--')
            plt.plot([Ha_peak_pos, Ha_peak_pos],[np.median(y_Ha), max(curve_Ha)], 'rx--')
            plt.show()

        for line in vis_lines:
            line_0 = line
            line *= (1. + Ha_dvrd/c)

            x, y = wvl, flx
            mindex, maxdex = np.searchsorted(x, (line-1, line+1))
            xr, yr = x[mindex:maxdex], y[mindex:maxdex]
            xrs = xr[2:-2]
            yrs = convolve(yr, box_kernel)[2:-2]
            peaks = find_peaks(-yrs)[0]
            peak_ys = yrs[peaks]
            rough_dex = peaks[np.argmin(peak_ys)]
            rough_loc = xrs[rough_dex]
            rough_hgt = yrs[rough_dex]
            mindex, maxdex = np.searchsorted(xrs, (rough_loc-vis_margin, rough_loc+vis_margin))
            xx, yy = xrs[mindex:maxdex], yrs[mindex:maxdex]

            mean = sum(xx * yy)/sum(yy)
            sigma = np.sqrt(sum(yy*(xx-mean)**2)/sum(yy))
            p0 = [np.median(yrs), rough_hgt, mean, sigma]
            coeff, var_matrix = curve_fit(gauss, xx, yy, p0=p0)
            curve = coeff[0] + coeff[1]*np.exp(-(xx-coeff[2])**2/(2.*coeff[3]**2))
            peak_pos = np.round(coeff[2], 4)
            wvl_shift = peak_pos-line_0
            vel_shift = (wvl_shift*c)/line_0

            csvwriter.writerow([file, date, time, phi, line_0, peak_pos, vel_shift])
            if vis_plot_flag:
                plt.plot(xr,yr, color='gray')
                plt.plot(xrs, yrs, 'k--')
                plt.plot(xx, yy, 'k')
                plt.plot([min(xrs),max(xrs)],[np.median(yrs), np.median(yrs)],'r--')
                plt.plot([peak_pos, peak_pos],[min(yr), max(yr)],'r--')
                plt.plot(xx,curve, color='r')
                plt.plot([line_0,line_0],[min(yr),max(yr)],'b--')
                plt.plot([line,line],[min(yr),max(yr)],'g--')
                plt.scatter(rough_loc, rough_hgt)
                plt.xlabel(str(line))
                plt.title(file)
                plt.show()

    for file in nir_files:
        hdu = fits.open('NIR/'+file)
        date, time = (hdu[0].header['DATE-OBS']).split('T')
        phi = hdu[0].header['PHI']
        wvl = hdu[1].data['WAVE'][0]
        flx = hdu[1].data['FLUX'][0]

        for line in nir_lines:
                x, y = wvl, flx
                mindex, maxdex = np.searchsorted(x, (line-2, line+2))
                xr, yr = x[mindex:maxdex], y[mindex:maxdex]
                xrs = xr[2:-2]
                yrs = convolve(yr, box_kernel)[2:-2]
                peaks = find_peaks(-yrs)[0]
                peak_ys = yrs[peaks]
                rough_dex = peaks[np.argmin(peak_ys)]
                rough_loc = xrs[rough_dex]
                rough_hgt = yrs[rough_dex]
                mindex, maxdex = np.searchsorted(xrs, (rough_loc-nir_margin, rough_loc+nir_margin))
                xx, yy = xrs[mindex:maxdex], yrs[mindex:maxdex]

                mean = sum(xx * yy)/sum(yy)
                sigma = np.sqrt(sum(yy*(xx-mean)**2)/sum(yy))
                p0 = [np.median(yrs), rough_hgt, mean, sigma]
                coeff, var_matrix = curve_fit(gauss, xx, yy, p0=p0)
                curve = coeff[0] + coeff[1]*np.exp(-(xx-coeff[2])**2/(2.*coeff[3]**2))
                peak_pos = np.round(coeff[2], 4)
                wvl_shift = peak_pos-line
                vel_shift = (wvl_shift*c)/line

                csvwriter.writerow([file, date, time, phi, line, peak_pos, vel_shift])

                if nir_plot_flag:
                    plt.plot(xr,yr, color='gray')
                    plt.plot(xrs, yrs, 'k--')
                    plt.plot(xx, yy, 'k')
                    plt.plot([min(xrs),max(xrs)],[np.median(yrs), np.median(yrs)],'r--')
                    plt.plot([peak_pos, peak_pos],[min(yr), max(yr)],'r--')
                    plt.plot(xx,curve, color='r')
                    plt.plot([line,line],[min(yr),max(yr)],'g--')
                    plt.scatter(rough_loc, rough_hgt)
                    plt.show()

vrads = []
phis = []

with open(f'{obj_name}_lines.csv', 'r') as file:
    for line in islice(file, 1, None): 
        splitstr = line.split(',')
        vrads.append(float(splitstr[-1]))
        phis.append(float(splitstr[3]))

p0_sin = [np.max(vrads)-np.min(vrads), 2*np.pi, 3, np.max(vrads)]

poptsin, pcovsin = curve_fit(model_sin, phis, vrads, p0=p0_sin, maxfev=5000)

sinx = np.linspace(0.,1.,1000)
siny = model_sin(sinx, poptsin[0], poptsin[1], poptsin[2], poptsin[3])

print(poptsin)

semi_amp = abs(poptsin[0])
v_sys = np.median(siny)

print('=====VELOCITY CURVE=====')
print(f'SEMI-AMPLITUDE      = {semi_amp} km/s')
print(f'SYSTEMIC VELOCITY   = {v_sys} km/s')

fig = plt.figure()
ax = fig.add_subplot()

print(np.max(siny))

ax.plot([0,1],[0,0],'k--', alpha=0.5)
ax.plot([0,1],[v_sys,v_sys],'r--', alpha=0.5)
ax.scatter(phis, vrads, color='k', marker='x')
ax.plot(sinx, siny, 'r')
ax.set_title(f'{obj_name} velocity curve')
ax.set_ylim(-1.1*np.max(siny), 1.1*np.max(siny))
ax.set_xlabel('Phase')
ax.set_ylabel('v_rad / km/s')
ax.text(0.95,0.95, f'SEMI-AMPLITUDE: {np.round(semi_amp, 2)} km/s',
        va='top', ha='right', transform=ax.transAxes)
ax.text(0.95,0.9, f'SYSTEMIC VELOCITY: {np.round(v_sys, 2)} km/s',
        va='top', ha='right', transform=ax.transAxes)
plt.show()
