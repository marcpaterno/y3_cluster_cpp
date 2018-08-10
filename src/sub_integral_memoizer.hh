#ifndef _Y3_CLUSTER_CPP_SUB_INTEGRAL_MEMOIZER_HH
#define _Y3_CLUSTER_CPP_SUB_INTEGRAL_MEMOIZER_HH

#include "integration_range.hh"
#include "interp_2d.hh"
#include <cubacpp/cubacpp.hh>
#include <memory>

namespace y3_cluster {
    // The below is basically a duplication of cubacpp for GSL - this is OK
    // because we need to integrate deterministically over 1 dimension,
    // and cuba does not do that.
    namespace {
        // See constructor of `MemoizedIntegral` for use
        template<typename F>
        std::shared_ptr<const Interp2D>
        make_memoized_integral(F const& f,
                               const IntegrationRange arange,
                               const IntegrationRange xrange,
                               const std::size_t xsteps,
                               const IntegrationRange yrange,
                               const std::size_t ysteps,
                               const double epsrel,
                               const double epsabs)
        {
            cubacpp::QAG qag(arange.transform(0.0), arange.transform(1.0),
                             GSL_INTEG_GAUSS61, 20);
            std::vector<double> xs(xsteps),
                                ys(ysteps),
                                values(xsteps * ysteps);

            for (auto i = 0u; i < xsteps; i++) {
                for (auto j = 0u; j < ysteps; j++) {
                    const double x = xrange.transform(((double) i) / (xsteps - 1)),
                                 y = yrange.transform(((double) j) / (ysteps - 1));

                    // Perform integration
                    auto result = qag.integrate([&f, x, y](double a) { return f(a, x, y); },
                                                epsabs,
                                                epsrel);

                    // Did something go wrong?
                    if (result.status != 0)
                        throw std::runtime_error("GSL Romberg Integration failed in MemoizedIntegration constructor");

                    // Record the result of this integration
                    xs.push_back(x);
                    ys.push_back(y);
                    values.push_back(result.value);
                }
            }

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
        MemoizedIntegration(F const& f,
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
