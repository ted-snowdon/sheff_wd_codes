# Just a couple of functions for reading .dat format spectra
# Not to be run on its own, just so I can copy these into other scripts

def read_spec(filename):
  # For dM model spectra as well as general use
    wvl = []
    flx = []
    err = []
    with open(filename, 'r') as file:
        for line in file:
            splitstr = [float(x) for x in line[:-2].split(' ') if x]
            wvl.append(splitstr[0])
            flx.append(splitstr[1])
            err.append(np.abs(splitstr[2]))
    wvl_out = np.asarray(wvl)
    flx_out = np.asarray(flx)
    err_out = np.asarray(err)
    return(wvl_out, flx_out, err_out)

def read_model(teff, logg):
  # For Koester WD model spectra
    filename = f'da{teff}_{logg}.dk.dat.txt'
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

