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
print(log_file)

if not '.log' in log_file:
    print('Not a valid .log output file')
    print('Closing...')
    quit

CCD1 = []
CCD2 = []
CCD3 = []

with open(log_file, 'r') as file:
    for _ in range(308):
        next(file)
    for line in file:
        if not '#' in line:
            splitline = line.split(' ')
            match splitline[0]:
                case '1':
                    CCD1.append(splitline)
                case '2':
                    CCD2.append(splitline)
                case '3':
                    CCD3.append(splitline)

fig, ax = plt.subplots(3,1)
fig.tight_layout()

mjd1 = []
counts1 = []
dcounts1 = []

for line in CCD1:
    mjd1.append(float(line[2]))
    counts1.append(float(line[15]))
    dcounts1.append(float(line[16]))

mjd2 = []
counts2 = []
dcounts2 = []

for line in CCD2:
    mjd2.append(float(line[2]))
    counts2.append(float(line[15]))
    dcounts2.append(float(line[16]))

mjd3 = []
counts3 = []
dcounts3 = []

for line in CCD3:
    mjd3.append(float(line[2]))
    counts3.append(float(line[15]))
    dcounts3.append(float(line[16]))

ax[0].errorbar(mjd1, counts1, yerr=dcounts1, color='k', ecolor='gray')
ax[0].set_title('CCD1')
ax[0].set_xticks([])
ax[1].errorbar(mjd2, counts2, yerr=dcounts2, color='k', ecolor='gray', label='CCD2')
ax[1].set_title('CCD2')
ax[1].set_ylabel('Counts')
ax[1].set_xticks([])
ax[2].errorbar(mjd3, counts3, yerr=dcounts3, color='k', ecolor='gray', label='CCD3')
ax[2].set_title('CCD3')
ax[2].set_xlabel('MJD')

plt.legend()
plt.show()