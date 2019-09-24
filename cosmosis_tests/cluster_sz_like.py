#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jackson O'Donnell
#   jacksonhodonnell@gmail.com

from __future__ import division, print_function
from cosmosis.datablock import names, option_section
from fitsio import FITS
import numpy as np
import os


def get_data(data_dir, mmin, mmax, zmin, zmax):
    '''
    Gets the data for a given mass-redshift bin.

    Returns 1-tuple of vectors, (ys, y errors).

    NB: This simply assumes the pipeline's y values are computed for the same
    angular bins. Fix this to be config-error proof!
    '''
    fname_fmt = os.path.join(data_dir,
                             'mass_ge{:.0e}_lt{:.0e}_z_ge{:g}_lt{:g}.fits')
    tbl = FITS(fname_fmt.format(mmin, mmax, zmin, zmax))[1]
    return tbl['y'].read(), tbl['yerr_jacknife'].read()


def setup(options):
    data_dir = options[option_section, 'data_dir']

    Mlows = options[option_section, 'M_low']
    Mhighs = options[option_section, 'M_high']
    zlows = options[option_section, 'z_low']
    zhighs = options[option_section, 'z_high']

    ys = []
    yerrs = []
    for Mlo, Mhi, zlo, zhi in zip(Mlows, Mhighs, zlows, zhighs):
        these_ys, these_yerrs = get_data(data_dir, Mlo, Mhi, zlo, zhi)
        ys.extend(these_ys)
        yerrs.extend(these_yerrs)

    return {'y': np.array(ys),
            'yerr': np.array(yerrs)}


def execute(block, config):
    results = block['compton_y_sims', 'compton_y_vals'].flatten()
    diff = results - config['y']

    # Covariancee matrix is diagonal in this case
    # (I'm making it square to easily extend to non-diagonal covariance,
    #  or to add other covariances)
    cov = np.diag(config['yerr']**2)
    # The precision matrix is the inverse of the covariance matrix
    precision = np.linalg.inv(cov)

    cov_mult = np.matmul(diff, np.matmul(precision, diff))
    # The normalization of a multivariate gaussian is:
    #   1 / \sqrt{ (2\pi)^k |Cov| }
    # Where k is the dimensions of the data vector

    # The following line is mathematically correct, but the determinant causes
    # an error because the covariances are too small
    # (The determinant == 0 so the np.log fails)
    # normalization = np.log(2 * np.pi * np.linalg.det(cov)) * (diff.size / 2)
    ln_det = np.log(np.diag(cov)).sum()
    normalization = (np.log(2 * np.pi) + ln_det) * (diff.size / 2)
    log_P = (- cov_mult / 2) - normalization

    block[names.likelihoods, 'compton_y_sims_like'] = log_P
    return 0


def cleanup(config):
    # nothing to do here
    pass
