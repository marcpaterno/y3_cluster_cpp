#ifndef _Y3_BESSEL_POLYNOMIAL_HH
#define _Y3_BESSEL_POLYNOMIAL_HH

#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>

#include "sin_cos_polynomial_integral.hh"
#include "primitives.hh"

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
    namespace {
        // For computing A_{nl} and B^n_{li}, defined above
        double
        product_l_n_k(int l, int n, int max)
        {
            double output = 1;
            for (auto k = 0; k < max; k++)
                output *= l + n - 2*k - 1;
            return output;
        }
    }

    // Performs the integral:
    //      \int dx x^n j_l(kx)
    // On the range x \in [range_min, range_max]
    double
    bessel_polynomial_integral(const int n, const unsigned l, const double k,
                               const double range_min, const double range_max)
    {
        // First, compute bessels up front
        std::vector<double> bessels_min(l), bessels_max(l);

        if (l > 0) {
            // WARNING: These values are not exact! See GSL docs:
            // https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_bessel_jl_array
            auto retval = gsl_sf_bessel_jl_array(l - 1, k * range_min, &bessels_min[0]);
            if (retval != GSL_SUCCESS)
              throw std::runtime_error("GSL error!");

            retval = gsl_sf_bessel_jl_array(l - 1, k * range_max, &bessels_max[0]);
            if (retval != GSL_SUCCESS)
              throw std::runtime_error("GSL error!");
        }

        // Compute the sum over B^n_{li}
        double running_sum = 0;
        for (auto i = 0u; i < l; i++) {
            const auto exponent = n - l + 1 + i;
            const double B = product_l_n_k(l, n, l - i - 1);
            running_sum -= B * (integer_pow(k*range_max, exponent) * bessels_max[i]
                              - integer_pow(k*range_min, exponent) * bessels_min[i]);
        }

        // Calculate \Delta X_{n - l - 1}
        const double xdiff = sinusoid_polynomial_integral(n - l - 1, k * range_min, k * range_max).first,
                     Anl = (l > 0) ? product_l_n_k(l, n, l) : 1.0;
        return integer_pow(k, -n - 1) * (Anl * xdiff + running_sum);
    }
}

#endif
