#ifndef _Y3_BESSEL_POLYNOMIAL_HH
#define _Y3_BESSEL_POLYNOMIAL_HH

#include <vector>

/* Contains routines for analytic integrations of the form:
 *
 *  I^n_l = \int dx x^n j_l(kx)
 *
 * Wher j_l is the `l`th spherical Bessel function of the first kind.
 *
 * Using recursion relations from Bloomfield et al. (https://arxiv.org/pdf/1703.06428.pdf),
 * Section 4
 *
 *  I^n_l = A_{nl} X_{n - l - 1} + \sum_{i = 0}^{l - 1} B^n_{li} x^{n - l + 1 + i} j_i
 *  A_{nl} = \product_{k = 0}^{l - 1} (l + n - 2k - 1)
 *  B^n_{li} = \product_{k = 0}^{l - i - 1} (l + n - 2k - 1)
 *
 * See `sin_cos_polynomial_integral.hh` for an explanation of X_n.
 */

namespace y3_cluster {
    // Performs the integral:
    //      \int dx x^n j_l(kx)
    // On the range x \in [range_min, range_max]
    double
    bessel_polynomial_integral(int n, unsigned l, double k,
                               double range_min, double range_max);

    // Computes I^n_l(x; k) on the range x = [range_min, range_max] for all l
    // from 0 to l, inclusive, and returns a vector of the results.
    // i.e., output[l] = I_n^l(x; k)
    std::vector<double>
    bessel_polynomial_integrals(int n, unsigned maxl, double k,
                                double range_min, double range_max);
}

#endif
