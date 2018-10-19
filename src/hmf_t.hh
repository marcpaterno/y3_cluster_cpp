#ifndef Y3_CLUSTER_CPP_HMF_T_HH
#define Y3_CLUSTER_CPP_HMF_T_HH

#include <datablock_reader.hh>
#include "interp_2d.hh"
#include <read_vector.hh>
#include <stdexcept>
#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cosmosis/datablock/ndarray.hh"
#include "/cosmosis/cosmosis/datablock/section_names.h"

#include <memory>

namespace y3_cluster {

  namespace {
      inline std::vector<double> 
      _adjust_to_log(cosmosis::DataBlock& db, const std::vector<double>& masses)
      {
          std::vector<double> output = masses;
          double omega_c = get_datablock<double>(db, "cosmological_parameters", "omega_c");
          double omega_b = get_datablock<double>(db, "cosmological_parameters", "omega_b");
          for (auto& x : output)
              x = std::log(x * (omega_c + omega_b));
          return output;
      }

      cosmosis::ndarray<double>
      log_ndarray(cosmosis::ndarray<double> input)
      {
          cosmosis::ndarray<double> output = input;
          if (output.ndims() != 2)
              throw std::runtime_error("expected 2d array");

          const auto dims = output.extents();
          for (auto i = 0u; i < dims[0]; i++)
              for (auto j = 0u; j < dims[1]; j++)
                  output(i, j) = std::log(output(i, j));
          return output;
      }
  }

  class HMF_t {
  public:
    HMF_t(std::shared_ptr<Interp2D const> nmz, double s, double q)
      : _nmz(nmz), _s(s), _q(q)
    {}

    using doubles = std::vector<double>;

    explicit HMF_t(cosmosis::DataBlock& sample)
      : _nmz(std::make_shared<Interp2D const>(
                  _adjust_to_log(sample, get_datablock<doubles>(sample, "mass_function", "m_h")),
                  get_datablock<doubles>(sample, "mass_function", "z"),
                  log_ndarray(get_datablock<cosmosis::ndarray<double>>(sample, "mass_function", "dndlnmh"))))
      , _s(get_datablock<double>(sample, "cluster_abundance", "hmf_s"))
      , _q(get_datablock<double>(sample, "cluster_abundance", "hmf_q"))
    {}

    double
    operator()(double lnM, double zt) const
    {
      return std::exp(_nmz->eval(lnM, zt)) * (_s * (lnM *0.4342944819 - 13.8124426028) + _q);
      //0.4342944819 is log(e)
    }

  private:
    std::shared_ptr<Interp2D const> _nmz;
    double _s;
    double _q;
  };
}

#endif
