#ifndef Y3_CLUSTER_COSMOSIS_SCALAR_EVALUATOR_MODULE_HH
#define Y3_CLUSTER_COSMOSIS_SCALAR_EVALUATOR_MODULE_HH

// CosmoSISScalarEvaluatorModule<I> — non-integrating counterpart of
// CosmoSISScalarIntegrationModule<I>.  Used when the integrand class
// already performs its own quadrature (e.g. fixed Gauss-Legendre) in
// set_sample and just needs to be queried per grid point.
//
// The integrand class I must supply:
//   using grid_point_t = std::array<double, N>;
//   static constexpr std::size_t n_outputs = K;
//
//   I(cosmosis::DataBlock&);
//   void set_sample(cosmosis::DataBlock&);
//   std::array<double, n_outputs> evaluate(grid_point_t const&);
//
//   static char const* module_label();          // config section
//   static std::array<char const*, n_outputs>
//     output_sections();                        // where to write vals
//   static grid_t<N> make_grid_points(cosmosis::DataBlock&);
//
// For each of the K outputs one 1-D ndarray is written to its own
// datablock section.  By default the subfield name is "vals"; an
// integrand may optionally supply
//   static std::array<char const*, n_outputs> output_names();
// to publish multiple named subfields in the same section (e.g.
// publishing {vals, rnd, cl} for the same Sigma_prj observable).

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"

#include "utils/make_grid_points.hh"

#include <array>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace y3_cluster {

  template <typename COSMOSISINTEGRAND>
  class CosmoSISScalarEvaluatorModule {
   public:
    using IntegrandType = COSMOSISINTEGRAND;
    using grid_point_t  = typename IntegrandType::grid_point_t;
    using my_grid_t     = grid_t<std::tuple_size<grid_point_t>::value>;

    explicit CosmoSISScalarEvaluatorModule(cosmosis::DataBlock& cfg);

    void execute(cosmosis::DataBlock& sample);

   private:
    IntegrandType integrand_;
    my_grid_t     grid_points_;
  };

}  // namespace y3_cluster

template <typename I>
y3_cluster::CosmoSISScalarEvaluatorModule<I>::CosmoSISScalarEvaluatorModule(
    cosmosis::DataBlock& cfg)
try : integrand_(cfg),
      grid_points_(IntegrandType::make_grid_points(cfg)) {
}
catch (cosmosis::Exception const&) {
  std::cerr << "\nDuring construction of a CosmoSISScalarEvaluatorModule, the "
               "lookup of some parameter failed.\n";
}
catch (std::exception const& e) {
  std::cerr << "\nDuring construction of a CosmoSISScalarEvaluatorModule, an "
               "std::exception was thrown:\n"
            << e.what() << '\n';
}
catch (...) {
  std::cerr << "\nUnknown exception while constructing a "
               "CosmoSISScalarEvaluatorModule.\n";
}

namespace y3_cluster {
  namespace evaluator_detail {
    // SFINAE: does I have a static output_names() returning
    // std::array<char const*, I::n_outputs>?
    template <typename I, typename = void>
    struct has_output_names : std::false_type {};

    template <typename I>
    struct has_output_names<
        I,
        std::void_t<decltype(I::output_names())>> : std::true_type {};

    template <typename I>
    inline std::array<char const*, I::n_outputs>
    get_output_names(std::true_type)
    {
      return I::output_names();
    }

    template <typename I>
    inline std::array<char const*, I::n_outputs>
    get_output_names(std::false_type)
    {
      std::array<char const*, I::n_outputs> names{};
      for (std::size_t k = 0; k != I::n_outputs; ++k) names[k] = "vals";
      return names;
    }
  }  // namespace evaluator_detail
}  // namespace y3_cluster

template <typename I>
void
y3_cluster::CosmoSISScalarEvaluatorModule<I>::execute(
    cosmosis::DataBlock& sample)
{
  integrand_.set_sample(sample);

  constexpr std::size_t K = IntegrandType::n_outputs;
  std::size_t const n     = grid_points_.size();

  std::array<std::vector<double>, K> buffers;
  for (auto& b : buffers) b.assign(n, 0.0);

  for (std::size_t i = 0; i != n; ++i) {
    auto out = integrand_.evaluate(grid_points_.points[i]);
    for (std::size_t k = 0; k != K; ++k) buffers[k][i] = out[k];
  }

  auto const sections = IntegrandType::output_sections();
  auto const names    = evaluator_detail::get_output_names<IntegrandType>(
      evaluator_detail::has_output_names<IntegrandType>{});
  for (std::size_t k = 0; k != K; ++k) {
    cosmosis::ndarray<double> vals(buffers[k], {n});
    sample.put_val(sections[k], names[k], vals);
  }
}

#endif
