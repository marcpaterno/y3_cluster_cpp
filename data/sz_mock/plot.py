#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jackson O'Donnell
#   jacksonhodonnell@gmail.com

from fitsio import FITS
import matplotlib.pyplot as plt
import pandas as pd

fname_fmt = 'mass_ge{:.0e}_lt{:.0e}_z_ge{:g}_lt{:g}.fits'

if __name__ == '__main__':
    df = {'mbin': [],
          'zbin': [],
          'R_norm': [],
          'y': [],
          'yerr_jacknife': []}

    Ms = [2e14, 5e14, 8e14, 1e15]
    zs = [0.2, 0.25, 0.3, 0.35, 0.4]
    for im, (mlo, mhi) in enumerate(zip(Ms[:], Ms[1:])):
        for iz, (zlo, zhi) in enumerate(zip(zs[:], zs[1:])):
            fname = fname_fmt.format(mlo, mhi, zlo, zhi)
            block = FITS(fname)[1]
            for name in df.keys():
                if name == 'mbin':
                    df[name].extend([im]*18)
                elif name == 'zbin':
                    df[name].extend([iz]*18)
                else:
                    df[name].extend(block[name].read())

    df = pd.DataFrame(df)
    for im in range(len(Ms) - 1):
        for iz in range(len(zs) - 1):
            print(df[(df.mbin == im) & (df.zbin == iz)].R_norm)

    fig, axs = plt.subplots(ncols=4, figsize=(12, 4), sharey=True)
    axs[0].set_ylabel('Compton-$y$', fontsize=20)
    for zbin, ax in zip(range(len(zs) - 1), axs):
        ax.loglog()
        ax.set_title('{:g} < z < {:g}'.format(zs[zbin], zs[zbin+1]))
        ax.set_xlabel('Radius (arcmin)')
        for mbin in range(len(Ms) - 1):
            bin_ = df[(df['mbin'] == mbin) & (df['zbin'] == zbin)]
            ax.errorbar(bin_['R_norm'], bin_['y'], yerr=bin_['yerr_jacknife'],
                        label='M={:.0e}-{:.0e}'.format(Ms[mbin], Ms[mbin+1]))
        ax.legend()
    plt.show()
