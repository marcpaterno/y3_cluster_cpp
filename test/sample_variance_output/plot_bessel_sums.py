#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt

bessel_sums = pd.read_csv('bessel_sums.csv')

plt.figure('Bessel Sums')

sums = ['sum_of_bessels_z{}_z{}'.format(i, j)
        for (i, j) in [(0, 0), (1, 1), (2, 2),
                       (0, 1), (0, 2), (1, 2)]]
for name in sums:
    bessel_sums.plot('k_h', name, ax=plt.gca())

plt.xscale('log')
plt.yscale('log')

plt.ylabel(r'$\sum_l R_{il}(k) R_{jl}(k) \Theta_l^2$')
plt.xlabel(r'$k_h (Mpc^{-1} h)$')

plt.legend()
plt.show()
