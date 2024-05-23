import os
import numpy as np

from scipy.optimize import curve_fit
from scipy.special import erf

basepath = "%s/data/" % '/global/common/software/des/jesteves/y3_cluster_cpp/'
fname_z = "int_lc_lt_des_t_zt_bins.txt"
fname_lt = 'int_lc_lt_des_t_lt_bins.txt'
fname_bin = "int_lc_lt_des_t_lbd_bin{}.txt"
zpiv = 0.35

def readIntBin(i):
    return np.array(np.genfromtxt(basepath+fname_bin.format(i)))

## Fit interpolator
def errFunction(x, xmin, xmax, sigma):
    base = np.sqrt(2) * sigma
    return 1*(erf((xmax - x) / base) - erf((xmin - x) / base)) / 2.0

compute_rms = lambda x, y: np.nanstd(y-x)#/np.nanmax(x)

def fitErrFunc(lbd_vec, int_lbd, verbose=False, p0=[25.844, 40.0429, 5.]):
    popt_erf, pcov_erf = curve_fit(errFunction, lbd_vec, int_lbd, maxfev = 5000, 
                                   p0=p0)
    
    ypred = errFunction(lbd_vec, *popt_erf)
    rms = compute_rms(int_lbd, ypred)#/np.nanmax(int_lbd)
    if verbose:
        a, b, c = np.round(popt_erf,3)
        print(f"Window Func fit parameters: xmin = {a}, xmax = {b}, sigma = {c}")
    return popt_erf, rms

def get_fitted_pars(bin_id):
    ## get paramaters for a given lambda bin
    intVec = np.array(int_bins[bin_id]).reshape(zt.size, lt.size)
    err = np.zeros(zt.size, dtype=float)
    pars = np.zeros((zt.size,3), dtype=float)

    for i,zm in enumerate(zt):
        yt = intVec[np.argmin(np.abs(zt-zm)),:]
        p, _ = fitErrFunc(lbd_vec, yt)
        pars[i] = p
    return pars

def fit_pars_to_poly(pars, deg=2, nth=1):
    pm = np.poly1d(np.polyfit(zt-zpiv, pars, deg))(zt-zpiv)
    outliers = np.abs(pars-pm) > nth*np.nanstd(pars-pm)
    coefs = np.polyfit(zt[~outliers]-zpiv, pars[~outliers], deg)
    return coefs

def interpolate_window_pars(pars, deg=2):
    return [np.flip(fit_pars_to_poly(pars[:,ix],deg)) for ix in range(3)]

def measure_interpolation_error(bin_id, coefs):
    models = [np.poly1d(np.flip(coefs[j])) for j in range(3)]
    intVec = np.array(int_bins[bin_id]).reshape(zt.size, lt.size)

    err = np.zeros(zt.size, dtype=float)
    for i,zm in enumerate(zt):
        guess = [m(zm-zpiv) for m in models]
        yt = intVec[np.argmin(np.abs(zt-zm)),:]
        ypre = errFunction(lbd_vec, *guess)
        err[i] = compute_rms(yt, ypre)
    return err

def print_pars_values(bin_id, coefs, err):
    coefs = np.round(coefs, 2)
    print(2*"--------------------")
    print(f"Lambda Bin {bin_id}")
    print(f"At pivot value z={zpiv}")
    print(f"lmin, lmax, sigma: {coefs[0][0], coefs[1][0], coefs[2][0]}")
    print(f"RMS Error: {np.interp(zpiv, zt, err)}")
    print(2*"--------------------")
    print("")

def print_pars(label,arr):
    l1 = "std::vector<std::vector<double>> %s = {"%label
    print(l1)
    sp = len(l1)*" "
    for row in arr:
        print(sp+"{", end="")
        for i in range(len(row)):
            print("{:.8f}".format(row[i]), end="")
            if i < len(row) - 1:
                print(", ", end="")
        print("},")
    print(sp+"}")
    print("")

## Read each lambda obs bin file
zt = np.array(np.genfromtxt(basepath+fname_z))
lt = np.array(np.genfromtxt(basepath+fname_lt))
int_bins = [readIntBin(i) for i in range(4)]

## Fit the window function to each lambda bin
lbd_vec = lt
lminV = []
lmaxV = []
sigmaV = []

for bin_id in range(4):
    pars = get_fitted_pars(bin_id)
    coefs = interpolate_window_pars(pars, deg=2)
    err = measure_interpolation_error(bin_id, coefs)
    print_pars_values(bin_id, coefs, err)
    lminV.append(coefs[0])
    lmaxV.append(coefs[1])
    sigmaV.append(coefs[2])

## print the parameters as cpp vectors
print_pars("lminV",lminV)
print_pars("lmaxV",lmaxV)
print_pars("sigmaV",sigmaV)