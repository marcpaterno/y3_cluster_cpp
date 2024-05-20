import cluster_toolkit as clusterwl
import numpy as np
import os

def sigma(R, M, c, om, klin, Plin, Rp, xi_mm):
    xi_nfw   = clusterwl.xi.xi_nfw_at_r(R, M, c, om)
    bias = 1
    xi_2halo = clusterwl.xi.xi_2halo(bias, xi_mm)
    xi_hm    = clusterwl.xi.xi_hm(xi_nfw, xi_2halo)
    Sigma_mm  = clusterwl.deltasigma.Sigma_at_R(Rp, R, xi_2halo, M, c, om)
    Sigma_nfw = clusterwl.deltasigma.Sigma_nfw_at_R(Rp, M, c, om)
    return Sigma_nfw, Sigma_mm, xi_nfw, xi_2halo

if __name__ == "__main__":
  # All of the test files have been commited into the git repo. This code is kept here for record.
  # If you'd like to run it yourself, run cosmosis ../deltasigma.ini first to generate the matter power spectrum files
  # test are done at M=3.199267137797384375e+14, z=0.2010101, c=5.0, om=0.3
  dirname = os.path.expandvars('${Y3_CLUSTER_CPP_DIR}/deltasigma/OUTPUT/matter_power_nl')
  k = np.loadtxt(os.path.join(dirname, "k_h.txt"))
  P = np.loadtxt(os.path.join(dirname, "p_k.txt"))
  zz = np.loadtxt(os.path.join(dirname, "z.txt"))
  ind, =np.where((zz > 0.2) & (zz < 0.202))
  P=P[ind, :]
  c = 5
  om = 0.3
  M=3.199267137797384375e+14

  NR = 1000
  R = np.logspace(-2, 3, NR, base=10) #Xi_hm MUST be evaluated to higher than BAO
  Rp = np.logspace(-2, 1.7, NR, base=10)
  xi_mm    = clusterwl.xi.xi_mm_at_r(R, k, P)

  SigmaNFW, SigmaMM, xiNFW, xiMM=sigma(R, M, c, om, k, P, Rp, xi_mm)
  
  dirname = os.path.expandvars('${Y3_CLUSTER_CPP_DIR}/deltasigma/test')
  ff=open(os.path.join(dirname, 'Sigma_mm_2halo.txt'), 'w')
  for ii in range(len(Rp)):
     ff.write('%f, %.12e \n'%(Rp[ii], SigmaMM[ii]))
  ff.close 

  ff=open(os.path.join(dirname, 'Sigma_nfwonly.txt'), 'w')
  for ii in range(len(Rp)):
     ff.write('%f, %.12e \n'%(Rp[ii], SigmaNFW[ii]))
  ff.close

  ff=open(os.path.join(dirname, 'xi_mm_2halo.txt'), 'w')
  for ii in range(len(R)):
      ff.write('%f,  %.12e \n'%(R[ii], xiMM[ii]))
  ff.close
  
  ff=open(os.path.join(dirname, 'xi_nfwonly.txt'), 'w')     
  for ii in range(len(R)):
      ff.write('%f,  %.12e \n'%(R[ii], xiNFW[ii]))
  ff.close

  #generate bias test files using linear power spectrum
  dirname = os.path.expandvars('${Y3_CLUSTER_CPP_DIR}/deltasigma/OUTPUT/matter_power_lin')
  klin = np.loadtxt(os.path.join(dirname, "k_h.txt"))
  Plin = np.loadtxt(os.path.join(dirname, "p_k.txt"))
  zz = np.loadtxt(os.path.join(dirname, "z.txt"))
  ff=open(os.path.expandvars('${Y3_CLUSTER_CPP_DIR}/deltasigma/test/bias.txt'), 'w')
  ff.write('## Mass(solar mass/h), redshift, bias  \n')
  M=10.0**(np.arange(12.0, 16.0, 0.2))
  ind, =np.where((zz > 0.2) & (zz < 0.202))
  P=Plin[ind, :]
  bias=np.zeros(len(M))
  for ii in range(len(M)):
      bias[ii]=clusterwl.bias.bias_at_M(M[ii], klin, P, om)
      ff.write('%f, %f, %.12e \n'%(M[ii], zz[ind], bias[ii]))
 
  ind, =np.where(zz < 0.101)
  P=Plin[ind, :]
  bias=np.zeros(len(M))
  for ii in range(len(M)):
      bias[ii]=clusterwl.bias.bias_at_M(M[ii], klin, P, om)
      ff.write('%f, %f, %.12e\n'%(M[ii], zz[ind], bias[ii]))
  ff.close
