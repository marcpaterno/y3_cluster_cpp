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

    // Expands sine & cosine to 12th order around x=0, and analytically
    // integrates with a polynomial
    std::pair<double, double>
    sinusoid_polynomial_taylor_approx(int n, const double range_min, const double range_max)
    {
        // n + 1 power because we are integrating
        double powmin = integer_pow(range_min, n + 1),
               powmax = integer_pow(range_max, n + 1),
               sinint_min_sum = 0,
               cosint_min_sum = 0,
               sinint_max_sum = 0,
               cosint_max_sum = 0;

        int factorial = 1;

        for (auto i = 0; i <= 12; i++) {
            factorial *= (i == 0) ? 1 : i;
            double sign = ((i / 2) % 2) ? -1 : 1,
                   min_item = sign * (powmin / (n + i + 1)) / factorial,
                   max_item = sign * (powmax / (n + i + 1)) / factorial;

            if ((n + i + 1) == 0) {
                min_item = sign * std::log(std::abs(range_min)) / factorial;
                max_item = sign * std::log(std::abs(range_max)) / factorial;
            }

            if (i % 2) {
                // Odd, sine
                sinint_min_sum += min_item;
                sinint_max_sum += max_item;
            } else {
                // Even, cosine
                cosint_min_sum += min_item;
                cosint_max_sum += max_item;
            }

            powmin *= range_min;
            powmax *= range_max;
        }

        return {sinint_max_sum - sinint_min_sum, cosint_max_sum - cosint_min_sum};
    }

    // Returns (\Delta X_n, \Delta Y_n)
    std::pair<double, double>
    sinusoid_polynomial_integral(int n, const double range_min, const double range_max)
    {
        if ((std::abs(range_max - range_min) < 1)
                && (range_max < 2)
                && (range_min > -2)) {
            return sinusoid_polynomial_taylor_approx(n, range_min, range_max);
        }

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

        // On range around zero, defer to taylor expansion
        if ((std::abs(range_max - range_min) < 2)
                && (range_max < 2)
                && (range_min > -2)) {
            for (auto i = 0; i < diff; i++) {
                const auto [sin, cos] = sinusoid_polynomial_taylor_approx(minn + i, range_min, range_max);
                output.first[i] = sin;
                output.second[i] = cos;
            }
            return output;
        }

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
