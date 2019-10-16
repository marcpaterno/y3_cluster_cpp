#ifndef _Y3_CLUSTER_SIN_COS_POLYNOMIAL_HH
#define _Y3_CLUSTER_SIN_COS_POLYNOMIAL_HH

#include <utility>
#include <vector>

/* Contains routines for analytic integrations of the form:
 *
 *  X_n = \int dx x^n sin(x)
 *  Y_n = \int dx x^n cos(x)
 *
 * Using recursion relations from Bloomfield et al.
 * (https://arxiv.org/pdf/1703.06428.pdf), Section 3
 */

namespace y3_cluster {
  // Expands sine & cosine to 12th order around x=0, and analytically
  // integrates with a polynomial
  std::pair<double, double> sinusoid_polynomial_taylor_approx(
    int n,
    const double range_min,
    const double range_max);

  // Returns (\Delta X_n, \Delta Y_n)
  std::pair<double, double> sinusoid_polynomial_integral(
    int n,
    const double range_min,
    const double range_max);

  // Computes \Delta X_n and \Delta Y_n on the range
  // x = [range_min, range_max] for all n from minn to maxn, inclusive.
  //
  // Returns a pair of vectors, ([\Delta X_{minn} ... \Delta X_{maxn}],
  //                             [\Delta Y_{minn} ... \Delta Y_{maxn}])
  std::pair<std::vector<double>, std::vector<double>>
  sinusoid_polynomial_integrals(int minn,
                                int maxn,
                                const double range_min,
                                const double range_max);
} // namespace y3_cluster

#endif // _Y3_CLUSTER_SIN_COS_POLYNOMIAL_HH
