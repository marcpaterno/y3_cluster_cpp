import numpy as np
import matplotlib.pyplot as plt

nc=np.genfromtxt('comp_nc.txt', delimiter=',')
ma=np.genfromtxt('comp_mass.txt', delimiter=',' )
ma=np.log10(ma/nc)
nc_m=np.genfromtxt('comp_nc_matteo.txt', delimiter=',')
ma_m=np.genfromtxt('comp_mass_matteo.txt', delimiter=',')
print nc.shape, ma.shape, nc_m.shape, ma_m.shape
print nc

fig, axs = plt.subplots(4, 1,sharex=True, figsize=[8, 15])
xx=np.array([20, 30, 45, 60])
axs[0].plot(xx, nc[0, :], 'bo', label='low redshift')
axs[0].plot(xx, nc[1, :], 'rs', label='medium redshift')
axs[0].plot(xx, nc[2, :], 'k*', label='highest redshift')
axs[0].plot(xx, nc_m[0, :], 'b', label='Matteo')
axs[0].plot(xx, nc_m[1, :], 'r')
axs[0].plot(xx, nc_m[2, :], 'k')
axs[1].plot(xx, nc[0, :]/nc_m[0, :], 'bo')
axs[1].plot(xx, nc[1, :]/nc_m[1, :], 'rs')
axs[1].plot(xx, nc[2, :]/nc_m[2, :], 'k*')


axs[2].plot(xx,ma[0, :], 'bo', label='low redshift')
axs[2].plot(xx,ma[1, :], 'rs', label='medium redshift')
axs[2].plot(xx,ma[2, :], 'k*', label='highest redshift')
axs[2].plot(xx,ma_m[0, :], 'b', label='Matteo')
axs[2].plot(xx,ma_m[1, :], 'r')
axs[2].plot(xx,ma_m[2, :], 'k')
axs[3].plot(xx,ma[0, :]-ma_m[0, :], 'bo')
axs[3].plot(xx,ma[1, :]-ma_m[1, :], 'rs')
axs[3].plot(xx,ma[2, :]-ma_m[2, :], 'k*')

axs[3].set_xlabel('richness bins')
axs[0].set_ylabel('NC')
axs[1].set_ylabel('NC Diff (Y/matteo)')
axs[2].set_ylabel('log M')
axs[3].set_ylabel('Log M Diff (Y-Matteo)')

axs[0].legend(loc=0)
axs[3].set_xlim(10, 80)
axs[3].set_ylim(-0.015, 0.015)
fig.tight_layout()
plt.savefig('temp.png')
plt.show()
