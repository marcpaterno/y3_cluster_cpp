// GSL-driven variant of ShearPrj: outer integral over z done by
// gsl_integration_qagp (adaptive piecewise Gauss-Kronrod with an
// explicit breakpoint at z=zob to handle the xi_nl(Δchi) cusp).
// Inner (lnM, θ) still use fixed GL.  Writes to sections
// sigma_prj_gsl / dsigma_prj_gsl / shear_prj_gsl so it can run
// alongside the default fixed-GL evaluator.
#include "models/sigma_prj_t.hh"
#include "utils/module_macros.hh"

using ShearPrjGsl = y3_cluster::ShearPrjGsl;

DEFINE_COSMOSIS_SCALAR_EVALUATOR_MODULE(ShearPrjGsl)
