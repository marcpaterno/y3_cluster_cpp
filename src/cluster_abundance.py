#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jackson O'Donnell
#   jacksonhodonnell@gmail.com

from __future__ import division, print_function
from cosmosis.datablock import option_section
from cython import Integrand
import numpy as np
import traceback


def process_ranges(block, name):
    edges = block[option_section, name + '_bins']
    return list(zip(edges[:-1], edges[1:]))


def setup(block):
    config = {}
    config['verbosity'] = block[option_section, 'verbosity']
    config['lo_bins'] = process_ranges(block, 'lo')
    config['zo_bins'] = process_ranges(block, 'zo')
    config['radii_bins'] = block[option_section, 'radii_bins']
    return config


def execute(sample, config):
    print('About to initialize integrand...')
    igrnd = Integrand(sample, config['radii_bins'],
                      config['lo_bins'], config['zo_bins'])
    print('Initialized integrand.')

    ctrd_res = igrnd.integrate_centered()
    if config['verbosity'] > 0:
        print('Centered:')
        print(ctrd_res)

    misctrd_res = igrnd.integrate_miscentered()

    if config['verbosity'] > 0:
        print('Miscentered:')
        print(misctrd_res)
    else:
        print('Centered and miscentered integration')

    # Store all of the results in the integrand in the datablock
    # Centered and miscentered values stored separately
    try:
        for name, res in [('centered', ctrd_res),
                          ('miscentered', misctrd_res)]:
            results = {'gamma_ts': [],
                       'cluster_counts': [],
                       'cluster_count_errors': [],
                       'cluster_biased_counts': []}

            for bin_ in res:
                for gt in bin_.gamma_ts:
                    results['gamma_ts'].append(gt)
                results['cluster_counts'].append(bin_.N)
                results['cluster_count_errors'].append(bin_.N_error)
                results['cluster_biased_counts'].append(bin_.Nb)

            for vname, value in results.items():
                sample['cluster_abundance',
                       name + '_' + vname] = np.array(value)
    except:
        traceback.print_exc()
        raise

    del igrnd
    return 0


def cleanup(config):
    pass
