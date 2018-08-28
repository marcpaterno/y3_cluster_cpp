#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt

bessel_sums = pd.read_csv('bessel_sums.csv')
bessel_sums['residual'] = (bessel_sums['sum_of_bessels'] - bessel_sums['manual_sum_of_bessels']) / bessel_sums['manual_sum_of_bessels']

fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
ax1.set_title('Bessel Sums')

bessel_sums.plot('k_h', 'sum_of_bessels', ax=ax1)
bessel_sums.plot('k_h', 'manual_sum_of_bessels', ax=ax1)
bessel_sums.plot('k_h', 'residual', ax=ax2)

ax1.set_xscale('log')
ax1.set_yscale('log')
ax2.set_yscale('log')

plt.legend()
plt.show()
