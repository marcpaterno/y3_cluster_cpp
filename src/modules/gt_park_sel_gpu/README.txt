In this folder we apply the Park et al. 2023 model to \Sigma:
    \tilde{\Sigma}(R) = \Pi(R) \Sigma(R)

In this module we compute the average \tilde{\Sigma(R)}
for a richness, redshift bin. 

The gamma_t is obtained in post-processing with the following:
    \Delta \Sigma = \bar{\Sigma}(<R) - \Sigma(R)

The \bar{\Sigma}(<R) is the cumulative \Sigma profile. 
Here we compute:
Iij = \int_Ri^Rj \Pi x \Sigma(R) R dR

which can be used later as:
\bar{\Sigma}(<Rk) = 2/Rk^2 \sum_{i=1}^{k} i I0i(k)

The description of the modules of this folder are:
avgCent: centered model of the \Sigma profile.
avgMiscApprox2: drop integration ovxer \lambda_c, Rmis and \theta,
                uses an interoplation table for delta sigma.

The \Pi x \Sigma module names are:
avgCentSigmaPark
avgMiscApprox2SigmaPark

The cumulative \bar{\Sigma} models are defined as:
\bar{\Sigma} = \int_0^Rp  \Pi x \Sigma(Rp) Rp dRp 

with module names:
avgCentSigmaCumPark
avgMiscApprox2SigmaCumPark