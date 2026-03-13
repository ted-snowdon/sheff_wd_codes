import hipercam as hcam
import matplotlib.pyplot as plt
from tkinter.filedialog import askopenfilename

# For quick plotting of ultracam reductions

#.log file column definitions
#0, 1 - CCD, nframe
#2, 3 - MJD, MJDok 
#4 - Exptime 
#5, 6 - mfwhm, mbeta 
#7, 8, 9, 10 - x_1, xe_1, y_1, ye_1 
#11, 12, 13, 14 - fwhm_1, fwhme_1, beta_1, betae_1 
#15, 16 - counts_1, countse_1 
#17, 18, 19 - sky_1, skye_1, nsky_1 
#20, 21, 22 - nrej_1, cmax_1, flag_1 
#23, 24, 25, 26 - x_2, xe_2, y_2, ye_2 
#27, 28, 29, 30 - fwhm_2, fwhme_2, beta_2, betae_2 
#31, 32 - counts_2 countse_2
#33, 34, 35 - sky_2, skye_2, nsky_2 
#36, 37, 38 - nrej_2, cmax_2, flag_2 

print('Select .log file: ')
log_file = askopenfilename()
hlog = hcam.hlog.Hlog.rascii(log_file)

t0 = 2458790.986815
p = 0.1147833978

def read_ccd(ccd):
    ccd_out = hcam.hlog.Tseries(
        t = hlog[ccd]['MJD'],
        y = hlog[ccd]['counts_1'],
        ye = hlog[ccd]['countse_1'],
        te = hlog[ccd]['Exptim'],
        cpy = True
    )
    
    ccd_out.mjd2tdb(
        position = '4:06:35.00 +10:00:59.0',
        telescope = 'Thai National Observatory',
        inplace = True
    )

    ccd_out.flag_outliers(
        sigma = 3.0,
        inplace = True)

    ccd_out.to_mag(
        inplace = True
    )

    return(ccd_out)

fig, ax = plt.subplots(3,2)
fig.tight_layout()

ccd1 = read_ccd('1')
ccd2 = read_ccd('2')
ccd3 = read_ccd('3')

ccd1p = ccd1.phase(
    t0 = t0,
    period = p,
    fold = True,
    inplace = False,
    sort = True,
)

ccd2p = ccd2.phase(
    t0 = t0,
    period = p,
    fold = True,
    inplace = False,
    sort = True,
)

ccd3p = ccd3.phase(
    t0 = t0,
    period = p,
    fold = True,
    inplace = False,
    sort = True,
)

ccd1.mplot(axes = ax[0,0], color = 'r')
ccd1p.mplot(axes = ax[0,1], color = 'r')
ccd2.mplot(axes = ax[1,0], color = 'g')
ccd2p.mplot(axes = ax[1,1], color = 'g')
ccd3.mplot(axes = ax[2,0], color = 'b')
ccd3p.mplot(axes = ax[2,1], color = 'b')

for axis in ax.flatten():
    axis.invert_yaxis()

plt.show()