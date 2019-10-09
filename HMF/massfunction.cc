#include "massfunction.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_sf.h"
#include "gsl/gsl_spline.h"
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#define rhocrit 2.77533742639e+11
#define del 1e-6
#define M_PIl 3.141592653589793238462643383279502884L /* pi */

#define d 1.97
#define e 1.0
#define f 0.51
#define g 1.228

#define ABSERR 0.0
#define RELERR 1e-5
#define delta_c 1.686 // Critical collapse density
#define rhocrit 2.77533742639e+11
// 1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are Msun h^2/Mpc^3
#define workspace_size 8000
#define KEY 4 // Used for GSL QAG function

typedef struct integrand_params {
  gsl_spline* spline;
  gsl_interp_accel* acc;
  double r;
  double* kp; // pointer to wavenumbers
  double* Pp; // pointer to P(k) array
  int Nk;     // length of k and P arrays
} integrand_params;

double
integrand(double lk, void* params)
{
  integrand_params pars = *(integrand_params*)params;
  double k = exp(lk);
  double x = k * pars.r;
  double P = gsl_spline_eval(pars.spline, k, pars.acc);
  double w = (sin(x) - x * cos(x)) * 3.0 / (x * x * x); // Window function
  return k * k * k * P * w * w;
}

double
M_to_R(double M, double Omega_m)
{
  // Lagrangian radius Mpc/h
  return pow(M / (1.33333333333 * M_PIl * rhocrit * Omega_m), 0.3333333333);
}

int
sigma2_at_R_arr(double R, double* k, double* P, int Nk, double* s2)
{
  // Initialize GSL things and the integrand structure.
  gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, Nk);
  gsl_interp_accel* acc = gsl_interp_accel_alloc();
  gsl_integration_workspace* workspace =
    gsl_integration_workspace_alloc(workspace_size);
  gsl_function F;
  integrand_params* params =
    (integrand_params*)malloc(sizeof(integrand_params));
  double lkmin = log(k[0]);
  double lkmax = log(k[Nk - 1]);
  double result, abserr;
  int i;
  gsl_spline_init(spline, k, P, Nk);
  params->spline = spline;
  params->acc = acc;
  params->kp = k;
  params->Pp = P;
  params->Nk = Nk;
  F.function = &integrand;
  F.params = params;
  for (i = 0; i < 1; i++) {
    params->r = R;
    gsl_integration_qag(&F,
                        lkmin,
                        lkmax,
                        ABSERR,
                        RELERR,
                        workspace_size,
                        KEY,
                        workspace,
                        &result,
                        &abserr);
    s2[i] = result / (2 * M_PIl * M_PIl);
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return 0;
}

int
sigma2_at_M_arr(double M, double* k, double* P, int Nk, double om, double* s2)
{
  double R = M_to_R(M, om);
  sigma2_at_R_arr(R, k, P, Nk, s2);
  return 0;
}

int
G_at_sigma_arr(double* sigma, double* G)
{
  // Compute the prefactor B
  double d2 = 0.5 * d;
  double gamma_d2 = gsl_sf_gamma(d2);
  double f2 = 0.5 * f;
  double gamma_f2 = gsl_sf_gamma(f2);
  double B = 2. / (pow(e, d) * pow(g, -d2) * gamma_d2 + pow(g, -f2) * gamma_f2);
  int i;
  for (i = 0; i < 1; i++) {
    G[i] = B * exp(-g / (sigma[i] * sigma[i])) *
           (pow(sigma[i] / e, -d) + pow(sigma[i], -f));
  }
  return 0;
}

int
dndM_sigma2_precomputed(double M,
                        double* sigma2,
                        double* sigma2_top,
                        double* sigma2_bot,
                        double Omega_m,
                        double* dndM)
{
  // This function exists to make emulator tests faster with sigma^2(M)
  // precomputed.
  double dndMconst = Omega_m * rhocrit * 0.5 / del; // normalization coefficient
  double* sigma = (double*)malloc(sizeof(double));
  double* Gsigma = (double*)malloc(sizeof(double));
  int i;
  for (i = 0; i < 1; i++) {
    sigma[i] = sqrt(sigma2[i]);
  }
  G_at_sigma_arr(sigma, Gsigma);
  for (i = 0; i < 1; i++) {
    dndM[i] =
      dndMconst * Gsigma[i] * log(sigma2_top[i] / sigma2_bot[i]) / (M * M);
  }
  free(sigma);
  free(Gsigma);
  return 0;
}

int
dndM_at_M_arr(double M, double* k, double* P, int Nk, double om, double* dndM)
{
  double* sigma2 = (double*)malloc(sizeof(double));
  double* sigma2_top = (double*)malloc(sizeof(double));
  double* sigma2_bot = (double*)malloc(sizeof(double));
  double M_top = M * (1 - del * 0.5);
  double M_bot = M * (1 + del * 0.5);
  sigma2_at_M_arr(M, k, P, Nk, om, sigma2);
  sigma2_at_M_arr(M_top, k, P, Nk, om, sigma2_top);
  sigma2_at_M_arr(M_bot, k, P, Nk, om, sigma2_bot);
  dndM_sigma2_precomputed(M, sigma2, sigma2_top, sigma2_bot, om, dndM);
  free(sigma2);
  free(sigma2_top);
  free(sigma2_bot);
  return 0;
}

double
dndM_at_M(double M, double* k, double* P, int Nk, double om)
{
  double* dndM = (double*)malloc(sizeof(double));
  double result;
  dndM_at_M_arr(M, k, P, Nk, om, dndM);
  result = dndM[0];
  free(dndM);
  return result;
}
