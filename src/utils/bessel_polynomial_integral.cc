#include "cubacpp/cubacpp.hh"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_bessel.h"
#include <cmath>

#include "bessel_polynomial_integral.hh"
#include "primitives.hh"
#include "sin_cos_polynomial_integral.hh"

using namespace y3_cluster;

// For computing A_{nl} and B^n_{li}, defined above
static double
product_l_n_k(int l, int n, int max)
{
  double output = 1;
  for (auto k = 0; k < max; k++) output *= l + n - 2 * k - 1;
  return output;
}

static std::vector<double>
bessels_array(const unsigned l, const double x)
{
  std::vector<double> output(l + 1);
  // WARNING: These values are not exact! See GSL docs:
  // https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_bessel_jl_array

  // Contrary to GSL docs, Steed's method seems to be faster. Maybe our maxl is
  // too small to see the benefit of the recursive relations.
  auto retval = gsl_sf_bessel_jl_steed_array(l, x, &output[0]);
  if (retval != GSL_SUCCESS) throw std::runtime_error("GSL error!");
  return output;
}

static double
bessel_polynomial_integral_quadrature(const int n,
                                      const unsigned l,
                                      const double k,
                                      const double range_min,
                                      const double range_max)
{
  cubacpp::QNG integrator(range_min, range_max);
  auto f = [=](double r) {
    return integer_pow(r, n) * gsl_sf_bessel_jl(l, k * r);
  };
  const auto res = integrator.integrate(f, 1e-5, 1e-18);
  if (res.status != 0)
    throw std::runtime_error("Bessel quadrature did not converge!");
  return res.value;
}

// Performs the integral:
//      \int dx x^n j_l(kx)
// On the range x \in [range_min, range_max], where j_l(x) has already been
// computed at k*range_min and k*range_max for l = 0, 1, 2, ... l - 1,
// and \Delta X_{n - l - 1} has been computed on the same range.
static double
bessel_polynomial_integral_precomputed(const int n,
                                       const unsigned l,
                                       const double k,
                                       const double range_min,
                                       const double range_max,
                                       const std::vector<double>& bessels_min,
                                       const std::vector<double>& bessels_max,
                                       const double Xn)
{
  if ((bessels_min.size() != bessels_max.size()) || (bessels_min.size() < l))
    throw std::runtime_error(
      "bessel_polynomial_integral_precomputed: Bad jl vector size");

  // The analytic integration is unstable for ranges very close to
  // zero, punt to a quadrature version for integrations within
  // the first zero of the Bessel function. See Bloomfield et al.,
  // Equation 15
  const auto diff = std::abs(k * (range_max - range_min));
  if (((l < 10) && (diff < (pi() + l))) || (diff < (4.75 + 1.05 * l)))
    return bessel_polynomial_integral_quadrature(n, l, k, range_min, range_max);

  // Compute the sum over B^n_{li}
  double running_sum = 0;
  for (auto i = 0u; i < l; i++) {
    const auto exponent = n - l + 1 + i;
    const double B = product_l_n_k(l, n, l - i - 1);
    running_sum -= B * (integer_pow(k * range_max, exponent) * bessels_max[i] -
                        integer_pow(k * range_min, exponent) * bessels_min[i]);
  }

  // Calculate \Delta X_{n - l - 1}
  const double Anl = (l > 0) ? product_l_n_k(l, n, l) : 1.0;
  return integer_pow(k, -n - 1) * (Anl * Xn + running_sum);
}

// Performs the integral:
//      \int dx x^n j_l(kx)
// On the range x \in [range_min, range_max]
double
y3_cluster::bessel_polynomial_integral(const int n,
                                       const unsigned l,
                                       const double k,
                                       const double range_min,
                                       const double range_max)
{
  // First, compute bessels up front
  std::vector<double> bessels_min, bessels_max;

  if (l > 0) {
    bessels_min = bessels_array(l - 1, k * range_min);
    bessels_max = bessels_array(l - 1, k * range_max);
  }

  // Then pass off to other function
  return bessel_polynomial_integral_precomputed(
    n,
    l,
    k,
    range_min,
    range_max,
    bessels_min,
    bessels_max,
    sinusoid_polynomial_integral(n - l - 1, k * range_min, k * range_max)
      .first);
}

// Computes I^n_l(x; k) on the range x = [range_min, range_max] for all l
// from 0 to l, inclusive, and returns a vector of the results.
// i.e., output[l] = I_n^l(x; k)
std::vector<double>
y3_cluster::bessel_polynomial_integrals(const int n,
                                        const unsigned maxl,
                                        const double k,
                                        const double range_min,
                                        const double range_max)
{
  // First, make output vector, and compute bessels up front
  std::vector<double> output(maxl + 1), bessels_min, bessels_max;

  if (maxl > 0) {
    bessels_min = bessels_array(maxl - 1, k * range_min);
    bessels_max = bessels_array(maxl - 1, k * range_max);
  }

  const auto sin_integral =
    sinusoid_polynomial_integrals(n - maxl - 1, n, k * range_min, k * range_max)
      .first;

  for (auto l = 0u; l <= maxl; l++)
    output[l] = bessel_polynomial_integral_precomputed(n,
                                                       l,
                                                       k,
                                                       range_min,
                                                       range_max,
                                                       bessels_min,
                                                       bessels_max,
                                                       sin_integral[maxl - l]);

  return output;
}
