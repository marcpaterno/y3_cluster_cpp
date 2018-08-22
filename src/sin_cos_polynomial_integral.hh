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

    // Returns (\Delta X_n, \Delta Y_n)
    std::pair<double, double>
    sinusoid_polynomial_integral(int n, const double range_min, const double range_max)
    {
        const double sin_min = std::sin(range_min),
                     sin_max = std::sin(range_max),
                     cos_min = std::cos(range_min),
                     cos_max = std::cos(range_max);

        double X_min = 0, X_max = 0, Y_min = 0, Y_max = 0;
        if (n >= 0) {
            // Set X, Y for base case
            X_min = -cos_min;
            X_max = -cos_max;
            Y_min = sin_min;
            Y_max = sin_max;

            for (auto i = 1; i <= n; i++) {
                const double X_min_new = i * Y_min - integer_pow(range_min, i) * cos_min,
                             X_max_new = i * Y_max - integer_pow(range_max, i) * cos_max,
                             Y_min_new = -i * X_min + integer_pow(range_min, i) * sin_min,
                             Y_max_new = -i * X_max + integer_pow(range_max, i) * sin_max;

                X_min = X_min_new;
                X_max = X_max_new;
                Y_min = Y_min_new;
                Y_max = Y_max_new;
            }
        } else if (n < 0) {
            // Set X, Y for base case (n = -1)
            X_min = gsl_sf_Si(range_min);
            X_max = gsl_sf_Si(range_max);
            Y_min = gsl_sf_Ci(range_min);
            Y_max = gsl_sf_Ci(range_max);

            for (auto i = -2; i >= n; i--) {
                const double X_min_new = (integer_pow(range_min, i + 1) * sin_min - Y_min) / (i + 1),
                             X_max_new = (integer_pow(range_max, i + 1) * sin_max - Y_max) / (i + 1),
                             Y_min_new = (integer_pow(range_min, i + 1) * cos_min + X_min) / (i + 1),
                             Y_max_new = (integer_pow(range_max, i + 1) * cos_max + X_max) / (i + 1);

                X_min = X_min_new;
                X_max = X_max_new;
                Y_min = Y_min_new;
                Y_max = Y_max_new;
            }
        }

        return {X_max - X_min, Y_max - Y_min};
    }
}

#endif // _Y3_CLUSTER_SIN_COS_POLYNOMIAL_HH
