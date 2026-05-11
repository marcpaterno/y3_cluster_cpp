// Single-pass evaluator for the 2-halo projection tangential shear.
// Co-computes Sigma_prj, DSigma_prj, and gamma_t^prj on one outer
// (zt, lnM, theta) GL loop, sharing every (z, M, theta)-dependent
// factor between the two NFW profiles.  Replaces the old trio of
// sigma_prj + dsigma_prj + red_shear_prj.py modules.
//
// Module label: "shear_prj".  Emits datablock sections sigma_prj,
// dsigma_prj, shear_prj — each with {vals, rnd, cl} subfields.
#include "models/sigma_prj_t.hh"
#include "utils/module_macros.hh"

using ShearPrjEvaluator = y3_cluster::ShearPrjEvaluator;

DEFINE_COSMOSIS_SCALAR_EVALUATOR_MODULE(ShearPrjEvaluator)
