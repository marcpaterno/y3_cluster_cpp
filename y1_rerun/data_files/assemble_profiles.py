import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import imshow

nvec = 12
ndim = 15
rad_vectors = np.zeros([nvec, ndim]) 
data_vectors = np.zeros([nvec, ndim]) 
data_vectors_err = np.zeros([nvec, ndim]) 
cov_vectors = np.zeros([nvec*ndim, nvec*ndim]) 


for z in range(3):
    for l in range(4):
        lam = l+3
        dat=np.genfromtxt('wl_data_files/full-unblind-v2-mcal-zmix_y1subtr_l%i_z%i_profile.dat'%(lam, z));
        print(dat[:, 0])
        cov = np.genfromtxt('SAC_files/SAC_z%i_l%i.txt'%(z, lam))
        start = (z*4 + l)
        stop = (z*4 + l)+1
        print(dat.shape, cov.shape, start, stop,  np.shape(rad_vectors[start:stop]))
        rad_vectors[start, :]=dat[:, 0]
        data_vectors[start, :]=dat[:, 1]
        data_vectors_err[start, :]=dat[:, 2]
        cov_vectors[(start*ndim):(stop*ndim), (start*ndim):(stop*ndim)]=cov
        if z == 2:
           plt.errorbar(rad_vectors[start, :], data_vectors[start, :], yerr=data_vectors_err[start, :])
#plt.errorbar(rad_vectors, data_vectors, yerr=np.diag(cov_vectors),c='k')
#plt.plot(rad_vectors, data_vectors_err, c='r')
#plt.plot(rad_vectors, data_vectors_err - np.sqrt(np.diag(cov_vectors)), '.', c='r')
#plt.plot(rad_vectors, np.sqrt(np.diag(cov_vectors)),c='k')
plt.xscale('log')
plt.yscale('symlog')
plt.show()

np.savetxt('wl_data_vector.txt', data_vectors)
np.savetxt('wl_data_vector_err.txt', data_vectors_err)
np.savetxt('wl_data_vector_radius.txt', rad_vectors)
np.savetxt('wl_cov.txt', cov_vectors)
imshow(np.log(cov_vectors))
plt.show()

