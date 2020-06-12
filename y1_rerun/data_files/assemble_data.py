import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import imshow

nvec = 12
ndim = 15
rad_vectors = np.zeros(nvec*ndim + ndim) 
data_vectors = np.zeros(nvec*ndim + ndim) 
data_vectors_err = np.zeros(nvec*ndim + ndim) 
cov_vectors = np.zeros([nvec*ndim+ndim, nvec*ndim+ndim]) 
ncs = np.genfromtxt('nc.txt')
ncs_err = np.genfromtxt('nc_err.txt')
print(ncs.shape)


for l in range(4):
    for z in range(3):
        data_vectors[z+l*3] = ncs[z, l]
        cov_vectors[z+l*3, z+l*3] = ncs_err[z, l]**2
        lam = l+3
        dat=np.genfromtxt('wl_data_files/full-unblind-v2-mcal-zmix_y1subtr_l%i_z%i_profile.dat'%(lam, z));
        cov = np.genfromtxt('SAC_files/SAC_z%i_l%i.txt'%(z, lam))
        start = (z + l*3)*ndim + ndim
        stop = (z + l * 3)*ndim + ndim + ndim
        print(dat.shape, cov.shape, start, stop,  np.shape(rad_vectors[start:stop]))
        rad_vectors[start:stop]=dat[:, 0]
        data_vectors[start:stop]=dat[:, 1]
        data_vectors_err[start:stop]=dat[:, 2]
        cov_vectors[start:stop, start:stop]=cov

#plt.errorbar(rad_vectors, data_vectors, yerr=data_vectors_err,c='r')
#plt.errorbar(rad_vectors, data_vectors, yerr=np.diag(cov_vectors),c='k')
plt.plot(rad_vectors, data_vectors_err, c='r')
plt.plot(rad_vectors, np.sqrt(np.diag(cov_vectors)),c='k')
plt.xscale('log')
plt.yscale('log')
plt.show()

np.savetxt('wl_data_vector.txt', data_vectors)
np.savetxt('wl_data_vector_radius.txt', rad_vectors)
np.savetxt('wl_cov.txt', cov_vectors)
imshow(np.log(cov_vectors))
plt.show()

