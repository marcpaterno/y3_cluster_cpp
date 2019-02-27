import cluster_toolkit as clusterwl
import numpy as np

def delta_sigma(R, M, c, om, klin, Plin, Rp, xi_mm, Rmis, fmis):
    xi_nfw   = clusterwl.xi.xi_nfw_at_R(R, M, c, om)
    bias = clusterwl.bias.bias_at_M(M, klin, Plin, om)
    print bias
    xi_2halo = clusterwl.xi.xi_2halo(bias, xi_mm)
    xi_hm    = clusterwl.xi.xi_hm(xi_nfw, xi_2halo)
    Sigma  = clusterwl.deltasigma.Sigma_at_R(Rp, R, xi_hm, M, c, om)
    Sigma_nfw = clusterwl.deltasigma.Sigma_nfw_at_R(Rp, M, c, om)
    #DeltaSigma = clusterwl.deltasigma.DeltaSigma_at_R(Rp, Rp, Sigma, M, c, om)
    #Sigma_exp  = clusterwl.miscentering.Sigma_mis_at_R(Rp, Rp, Sigma, M, c, om, Rmis, kernel="exponential")
    #DeltaSigma_exp = clusterwl.miscentering.DeltaSigma_mis_at_R(Rp, Rp, Sigma_exp)
    #DeltaSigma_comb = DeltaSigma#*(1-fmis)+DeltaSigma_exp*fmis
    return Sigma

if __name__ == "__main__":
  k = np.loadtxt("test_data/matter_power_nl/k_h.txt")
  P = np.loadtxt("test_data/matter_power_nl/p_k.txt")
  zz = np.loadtxt("test_data/matter_power_nl/z.txt")
  print k.shape, P.shape, zz.shape
  ind, =np.where((zz > 0.2) & (zz < 0.202))
  P=P[ind, :]

  NR = 1000
  R = np.logspace(-2, 3, NR, base=10) #Xi_hm MUST be evaluated to higher than BAO
  Rp = np.logspace(-2, 2.4, NR, base=10)
  c = 5
  om = 0.3
  xi_mm    = clusterwl.xi.xi_mm_at_R(R, k, P)

  #truth
  M=10.0**(14.5)
  Rmis = 0.15 #Mpc/h
  fmis = 0.32
  DeltaSigma_truth=delta_sigma(R, M, c, om, k, P, Rp, xi_mm, Rmis, fmis)
  
  #fmis uncertainty
  ff=open('Sigma_14.5_tot.txt', 'w')
  ff.write('## log_M (truth 14.5),  C(truth 5) redshift %f \n## R comoving (Mpc/h) Sigma  \n'%(zz[ind]))
  for ii in range(len(R)):
     ff.write('%f, %f \n'%(Rp[ii], DeltaSigma_truth[ii]))
  ff.close 
