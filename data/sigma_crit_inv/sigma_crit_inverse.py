#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jackson O'Donnell
#   jacksonhodonnell@gmail.com

from __future__ import division, print_function
from astropy.constants import c, G
from astropy.cosmology import FlatLambdaCDM
import numpy as np
import os

default_cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.7, Ob0=0.06)


def sig_crit_inv(zs, zl, cosmo=default_cosmo):
    """
    Computes the inverse of the weak lensing critical surface density, for a
    source at redshift zs and lens at redshift zl.

    Returns values in units of Mpc^2 / M_{sol}.
    """
    zs = np.array(zs)
    zl = np.array(zl)
    const = 4 * np.pi * (G / (c*c)).to('Mpc/solMass')
    dl = cosmo.angular_diameter_distance(zl)
    ds = cosmo.angular_diameter_distance(zs)
    dls = cosmo.angular_diameter_distance_z1z2(zl, zs)
    result = const * dl * dls / ds
    result[result < 0] = 0
    return result


if __name__ == '__main__':
    zl = np.arange(0.01, 1.005, 0.01)
    zs = np.copy(zl)
    folder = os.path.dirname(__file__)
    np.savetxt(os.path.join(folder, 'zl.txt'), zl,
               header='Lens redshifts for sigma_crit_inv test')
    np.savetxt(os.path.join(folder, 'zs.txt'), zs,
               header='Source redshifts for sigma_crit_inv test')

    zl, zs = np.meshgrid(zl, zs)
    msg = '\n'.join(['Inverse of critical surface density',
                     'Units: Mpc2 / M_{sol}'])
    np.savetxt(os.path.join(folder, 'sigma_crit_inverse.txt'),
               sig_crit_inv(zs, zl),
               header=msg)
    # TODO: try with other cosmologies?
