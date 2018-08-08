#ifndef _Y3_CLUSTER_CPP_SUB_INTEGRAL_MEMOIZER_HH
#define _Y3_CLUSTER_CPP_SUB_INTEGRAL_MEMOIZER_HH

#include <gsl/gsl_errno.h>
// For 1D integration routines
#include <gsl/gsl_integration.h>
#include <memory>

namespace y3_cluster {
    // The below is basically a duplication of cubacpp for GSL - this is OK
    // because we need to integrate deterministically over 1 dimension,
    // and cuba does not do that.
    namespace {
        template<typename F>
        double
        gsl_integrand(double x, void *params)
        {
            F *f = static_cast<F*>(params);
            return (*f)(x);
        }

        template<typename F>
        gsl_function_struct
        make_gsl_integrand(F *f)
        {
            gsl_function out;
            out.function = gsl_integrand<F>;
            out.params = static_cast<void*>(f);
            return out;
        }

        // See constructor of `MemoizedIntegral` for use
        template<typename F>
        std::shared_ptr<const Interp2D>
        make_memoized_integral(F f,
                               const IntegrationRange arange,
                               const IntegrationRange xrange,
                               const std::size_t xsteps,
                               const IntegrationRange yrange,
                               const std::size_t ysteps,
                               const double epsrel,
                               const double epsabs)
        {
            // Max iterations is 2**n, so 20 is reasonable
            auto *integrator = gsl_integration_romberg_alloc(20);
            if (!integrator)
                throw std::runtime_error("Failure to allocate GSL romberg integrator");

            std::vector<double> xs(xsteps),
                                ys(ysteps),
                                values(xsteps * ysteps);

            for (auto i = 0u; i < xsteps; i++) {
                for (auto j = 0u; j < ysteps; j++) {
                    // Setup
                    const double x = xrange.transform(((double) i) / (xsteps - 1)),
                                 y = yrange.transform(((double) j) / (ysteps - 1));
                    double result;
                    std::size_t neval;

                    // Perform integration
                    auto fn = [&f, x, y](double a) { return f(a, x, y); };
                    auto integrand = make_gsl_integrand(&fn);
                    auto retval = gsl_integration_romberg(&integrand,
                                                          arange.transform(0.0),
                                                          arange.transform(1.0),
                                                          epsabs,
                                                          epsrel,
                                                          &result,
                                                          &neval,
                                                          integrator);

                    // Did something go wrong?
                    if (retval != GSL_SUCCESS)
                        throw std::runtime_error("GSL Romberg Integration failed in MemoizedIntegration constructor");

                    // Record the result of this integration
                    xs.push_back(x);
                    ys.push_back(y);
                    values.push_back(result);
                }
            }

            // Need to free our integrator, as this is C
            gsl_integration_romberg_free(integrator);
            return std::make_shared<const Interp2D>(xs, ys, values);
        }
    }

    /* Class to memoize a 1D integration.
     *
     * For integrals of the form:
     *
     *   A(x, y) = \int_{\Delta a} F(a, x, y) da
     *
     * This class will evaluate A at a grid of (x, y) points, and at future
     * calls via `operator()` will interpolate over the results.
     */
    class MemoizedIntegration {
    private:
        std::shared_ptr<const Interp2D> results;
    public:
        MemoizedIntegration() = delete;

        template<typename F>
        MemoizedIntegration(F f,
                            const IntegrationRange arange,
                            const IntegrationRange xrange,
                            const std::size_t xsteps,
                            const IntegrationRange yrange,
                            const std::size_t ysteps,
                            const double epsrel=1e-5,
                            const double epsabs=1e-12)
            : results(make_memoized_integral(f, arange, xrange, xsteps, yrange, ysteps, epsrel, epsabs))
        {}

        double operator()(double x, double y) const
        {
            return results->eval(x, y);
        }
    };
}

#endif // _Y3_CLUSTER_CPP_SUB_INTEGRAL_MEMOIZER_HH
