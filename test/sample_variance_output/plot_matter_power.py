#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt

matter_power = pd.read_csv('matter_power.csv')

z1 = matter_power[matter_power['z'] == 0.1].dropna()

fig = plt.figure('Matter Power')
ax = fig.gca()
#z1.plot('k_h', 'p_k', ax=ax)
z1.plot('k_h', 'p_k_k_cubed', ax=ax)

plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
