import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from scipy.stats import chisquare
import spectres

class _Gdl:

    def __init__(self, vsini, epsilon):
        """
        Calculate the broadening profile.

        Parameters
        ----------
        vsini : float
            Projected rotation speed of the star [km/s]
        epsilon : float
            Linear limb-darkening coefficient
        """
        self.vc = vsini / 299792.458
        self.eps = epsilon

    def gdl(self, dl, refwvl, dwl):
        """
        Calculates the broadening profile.

        Parameters
        ----------
        dl : array
            'Delta wavelength': The distance to the reference point in
            wavelength space [A].
        refwvl : array
            The reference wavelength [A].
        dwl : float
            The wavelength bin size [A].

        Returns
        -------
        Broadening profile : array
            The broadening profile according to Gray. 
        """
        self.dlmax = self.vc * refwvl
        self.c1 = 2.*(1. - self.eps) / \
            (np.pi * self.dlmax * (1. - self.eps/3.))
        self.c2 = self.eps / (2. * self.dlmax * (1. - self.eps/3.))
        result = np.zeros(len(dl))
        x = dl/self.dlmax
        indi = np.where(np.abs(x) < 1.0)[0]
        result[indi] = self.c1 * \
            np.sqrt(1. - x[indi]**2) + self.c2*(1. - x[indi]**2)
        # Correct the normalization for numeric accuracy
        # The integral of the function is normalized, however, especially in the case
        # of mild broadening (compared to the wavelength resolution), the discrete
        # broadening profile may no longer be normalized, which leads to a shift of
        # the output spectrum, if not accounted for.
        result /= (np.sum(result) * dwl)
        return result

def read_spec(filename):
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

def fastRotBroad(wvl, flux, epsilon, vsini, effWvl=None):
    # Wavelength binsize
    dwl = wvl[1] - wvl[0]

    if effWvl is None:
        effWvl = np.mean(wvl)

    gdl = _Gdl(vsini, epsilon)

    # The number of bins needed to create the broadening kernel
    binnHalf = int(np.floor(((vsini / 299792.458) * effWvl / dwl))) + 1
    gwvl = (np.arange(4*binnHalf) - 2*binnHalf) * dwl + effWvl
    # Create the broadening kernel
    dl = gwvl - effWvl
    g = gdl.gdl(dl, effWvl, dwl)
    # Remove the zero entries
    indi = np.where(g > 0.0)[0]
    g = g[indi]

    result = np.convolve(flux, g, mode="same") * dwl
    return result

def lazychisq(o, e):
    # Scipy chisq is being weird so i'm gonna do it myself
    chisq = []

    for i in range(len(o)):
        chisq.append((o[i]-e[i])**2/e[i])

    return(sum(chisq[10:-10]))

wvl0, flx0, err0 = read_spec('ZTFJ040627.23+095826.97_coadded.dat')

standards = [f for f in listdir() if '_VIS.dat' in f]
standards.sort()

vsinis = np.arange(100,131,1)

for file in standards[6:10]:

    wvl1, flx1, err1 = read_spec(file)
    wvl1 /= 10.

    cutmin, cutmax = np.searchsorted(wvl1, (min(wvl0), max(wvl0)))
    wvl1 = wvl1[cutmin:cutmax]
    flx1 = flx1[cutmin:cutmax]
    err1 = err1[cutmin:cutmax]

    contmin, contmax = np.searchsorted(wvl1, (808.5, 816.5))
    cont_flx = np.median(flx1[contmin:contmax])

    flx1 /= cont_flx

    ruffmin, ruffmax = np.searchsorted(wvl1, (818.0,818.6))
    ruff_wvl, ruff_flx = wvl1[ruffmin:ruffmax], flx1[ruffmin:ruffmax]
    ruff_spot = ruff_wvl[np.argmin(ruff_flx)]

    dwvl = ruff_spot-818.326
    dvrd = (dwvl*299792.458)/818.326

    print(dvrd)

    wvl1 *= (1.-dvrd/299792.458)

    profmin, profmax = np.searchsorted(wvl1, (ruff_spot-1.0, ruff_spot+1.0))
    prof_wvl, prof_flx = wvl1[profmin:profmax], flx1[profmin:profmax]

    chisqs = []
    chisq_clip = np.searchsorted(wvl0, (817.75, 820.25))
    for vsini in vsinis:
        broad_prof = fastRotBroad(wvl1, flx1, 0.5, vsini)

        # Rebin broadened profile to observed wvl scale
        broad_prof_rebin = spectres.spectres(new_wavs = wvl0,
                                             spec_wavs = wvl1,
                                             spec_fluxes = broad_prof,
                                             spec_errs = np.zeros(len(broad_prof)),
                                             fill=1.0,
                                             verbose=False)[0]
        obs = broad_prof_rebin[chisq_clip[0]:chisq_clip[1]]
        exp = flx0[chisq_clip[0]:chisq_clip[1]]
        phleb = wvl0[chisq_clip[0]:chisq_clip[1]]
        
        #chisqs.append(lazychisq(obs, exp))
        chisqs.append(lazychisq(broad_prof_rebin, flx0))
    
    best_vsini = np.argmin(chisqs)
    print(file, vsinis[best_vsini], chisqs[best_vsini])

    plot_broad_prof = fastRotBroad(wvl1, flx1, 0.5, vsinis[best_vsini])
    
    fig, ax = plt.subplots()
    ax.plot(wvl0, flx0, 'k')
    ax.plot(wvl1, flx1, 'r', alpha=0.3)
    ax.plot(wvl1, plot_broad_prof, 'r')
    ax.text(0.95,0.95,str(vsinis[best_vsini]), va='top', ha='right', transform=ax.transAxes)
    ax.set_title(file)
    #ax.set_xlim(817,821)
    plt.show()