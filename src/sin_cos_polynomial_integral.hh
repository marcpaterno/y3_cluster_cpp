#ifndef _Y3_CLUSTER_SIN_COS_POLYNOMIAL_HH
#define _Y3_CLUSTER_SIN_COS_POLYNOMIAL_HH

#include <cmath>
#include <utility>
#include <gsl/gsl_sf_expint.h>

#include "primitives.hh"

/* Contains routines for analytic integrations of the form:
 *
 *  X_n = \int dx x^n sin(x)
 *  Y_n = \int dx x^n cos(x)
 *
 * Using recursion relations from Bloomfield et al. (https://arxiv.org/pdf/1703.06428.pdf),
 * Section 3
 */

namespace y3_cluster {
    namespace {
        // Helper function for recursion relation
        std::pair<double, double>
        sinusoid_polynomial_iteration(const std::pair<double, double> curstep, const int n,
                                               const double x, const double cosx, const double sinx)
        {
            if (n > 0)
                return {n * curstep.second - integer_pow(x, n) * cosx,
                        -n * curstep.first + integer_pow(x, n) * sinx};

            if (n < -1)
                return {(integer_pow(x, n + 1) * sinx - curstep.second) / (n + 1),
                        (integer_pow(x, n + 1) * cosx + curstep.first) / (n + 1)};

            throw std::runtime_error("n shouldn't be a base case");
        }
    }

    // Returns (\Delta X_n, \Delta Y_n)
    std::pair<double, double>
    sinusoid_polynomial_integral(int n, const double range_min, const double range_max)
    {
        const double sin_min = std::sin(range_min),
                     sin_max = std::sin(range_max),
                     cos_min = std::cos(range_min),
                     cos_max = std::cos(range_max);

        std::pair<double, double> vals_min{-cos_min, sin_min},
                                  vals_max{-cos_max, sin_max};
        if (n >= 0) {
            for (auto i = 1; i <= n; i++) {
                vals_min = sinusoid_polynomial_iteration(vals_min, i, range_min, cos_min, sin_min);
                vals_max = sinusoid_polynomial_iteration(vals_max, i, range_max, cos_max, sin_max);
            }
        } else if (n < 0) {
            // Set X, Y for base case (n = -1)
            // Ci = - \int_x^{+inf} cos(t)/t dt
            // Ci = - \int (odd function) dt
            //  therefore, Ci is _even_
            vals_min = {gsl_sf_Si(range_min), gsl_sf_Ci(std::abs(range_min))};
            vals_max = {gsl_sf_Si(range_max), gsl_sf_Ci(std::abs(range_max))};

            for (auto i = -2; i >= n; i--) {
                vals_min = sinusoid_polynomial_iteration(vals_min, i, range_min, cos_min, sin_min);
                vals_max = sinusoid_polynomial_iteration(vals_max, i, range_max, cos_max, sin_max);
            }
        }

        return {vals_max.first - vals_min.first, vals_max.second - vals_min.second};
    }

    // Computes \Delta X_n and \Delta Y_n on the range
    // x = [range_min, range_max] for all n from minn to maxn, inclusive.
    //
    // Returns a pair of vectors, ([\Delta X_{minn} ... \Delta X_{maxn}],
    //                             [\Delta Y_{minn} ... \Delta Y_{maxn}])
    std::pair<std::vector<double>, std::vector<double>>
    sinusoid_polynomial_integrals(int minn, int maxn, const double range_min, const double range_max)
    {
        if (minn > maxn)
            throw std::runtime_error("sinusoid_polynomial_integrals: Bad range");

        int diff = maxn - minn, pivot = -minn;
        std::pair<std::vector<double>, std::vector<double>>
            output = {std::vector<double>(diff + 1),
                      std::vector<double>(diff + 1)};

        const double sin_min = std::sin(range_min),
                     sin_max = std::sin(range_max),
                     cos_min = std::cos(range_min),
                     cos_max = std::cos(range_max);

        std::pair<double, double> vals_min{-cos_min, sin_min},
                                  vals_max{-cos_max, sin_max};
        if (maxn >= 0) {
            if (minn <= 0) {
                output.first[pivot] = vals_max.first - vals_min.first;
                output.second[pivot] = vals_max.second - vals_min.second;
            }

            for (auto i = 1; i <= maxn; i++) {
                vals_min = sinusoid_polynomial_iteration(vals_min, i, range_min, cos_min, sin_min);
                vals_max = sinusoid_polynomial_iteration(vals_max, i, range_max, cos_max, sin_max);

                if (i >= minn) {
                    output.first[pivot + i] = vals_max.first - vals_min.first;
                    output.second[pivot + i] = vals_max.second - vals_min.second;
                }
            }
        }

        if (minn < 0) {
            // Set X, Y for base case (n = -1)
            // Ci = - \int_x^{+inf} cos(t)/t dt
            // Ci = - \int (odd function) dt
            //  therefore, Ci is _even_
            vals_min = {gsl_sf_Si(range_min), gsl_sf_Ci(std::abs(range_min))};
            vals_max = {gsl_sf_Si(range_max), gsl_sf_Ci(std::abs(range_max))};

            if (maxn >= 0) {
                output.first[pivot - 1] = vals_max.first - vals_min.first;
                output.second[pivot - 1] = vals_max.second - vals_min.second;
            }

            for (auto i = -2; i >= minn; i--) {
                vals_min = sinusoid_polynomial_iteration(vals_min, i, range_min, cos_min, sin_min);
                vals_max = sinusoid_polynomial_iteration(vals_max, i, range_max, cos_max, sin_max);

                if (i <= maxn) {
                    output.first[pivot + i] = vals_max.first - vals_min.first;
                    output.second[pivot + i] = vals_max.second - vals_min.second;
                }
            }
        }

        return output;
    }
}

#endif // _Y3_CLUSTER_SIN_COS_POLYNOMIAL_HH
