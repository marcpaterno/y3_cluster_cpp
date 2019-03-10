#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jackson O'Donnell
#   jacksonhodonnell@gmail.com

from __future__ import division, print_function
from cosmosis.datablock import option_section
from cython import Integrand


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

    print('About to integrate...')
    ctrd_res = igrnd.integrate_centered()
    if config['verbosity'] > 0:
        print('Centered:')
        print(ctrd_res)

    misctrd_res = igrnd.integrate_miscentered()
    print('Finished integration.')

    if config['verbosity'] > 0:
        print('Miscentered:')
        print(misctrd_res)
    else:
        print('Centered and miscentered integration')

    del igrnd
    return 0


def cleanup(config):
    pass


if __name__ == '__main__':
    pass
