#ifndef Y3_CLUSTER_SIG_T_HH
#define Y3_CLUSTER_SIG_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "ez.hh"
#include "primitives.hh"
#include "interp_2d.hh"
#include "interp_1d.hh"
#include <datablock_reader.hh>
#include <cubacpp/cubacpp.hh>
#include <math.h>
#include <memory>
#include <cmath>
#include <iostream>

namespace y3_cluster
{
  class SIG_y3 {
  private:
    std::shared_ptr<Interp2D const> _xi1;
    std::shared_ptr<Interp2D const> _xi2;
    std::shared_ptr<Interp2D const> _bias;
    double _om;

  public:
    SIG_y3(std::shared_ptr<Interp2D const> xi1, 
              std::shared_ptr<Interp2D const> xi2, 
              std::shared_ptr<Interp2D const> bias, double om)
              : _xi1(xi1), _xi2(xi2), _bias(bias), _om(om) {}

    using doubles = std::vector<double>;

    // TODO: This needs to be reading cosmosis datablock parameters
    explicit SIG_y3(cosmosis::DataBlock& sample)
      : _xi1(std::make_shared<Interp2D const>(
          get_datablock<doubles>(sample, "cluster_abundance", "m_h"),
          get_datablock<doubles>(sample, "cluster_abundance", "R_perp"),
          get_datablock<cosmosis::ndarray<double>>(sample, "cluster_abundance", "xi_nfw")))
      , _xi2(std::make_shared<Interp2D const>(
          get_datablock<doubles>(sample, "cluster_abundance", "z"),
          get_datablock<doubles>(sample, "cluster_abundance", "R_perp"),
          get_datablock<cosmosis::ndarray<double>>(sample, "cluster_abundance", "xi_nl")))
      , _bias(std::make_shared<Interp2D const>(
          get_datablock<doubles>(sample, "tinker_bias_function", "ln_mass_h"),
          get_datablock<doubles>(sample, "tinker_bias_function", "z"),
          get_datablock<cosmosis::ndarray<double>>(sample, "tinker_bias_function", "bias")))
      , _om(get_datablock<double>(sample, "cosmological_parameters", "omega_m")) 
    {}

    double
    operator()(double R, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    { 

       double rho_m = 0.277526157*_om;
       cubacpp::QAG integrator(-50, 50, GSL_INTEG_GAUSS61, 50000000);
       const auto res = integrator.integrate([=](double chi) {
		 //double chi=std::exp(lnchi);
		 double r=sqrt(R*R+chi*chi);
		 //if (r > 100.0) r = 100.0;
                 double xi1v = _xi1->eval(std::exp(lnM), r);
                 //modify the following line after testing the bias module
		 //double xi2v = _bias->eval(lnM, zt) * _xi2->eval(zt, r); 
                 double xi2v = 2.97410387975 * _xi2->eval(zt, r);		 
                 double rest = xi1v;
                 if (xi1v <= xi2v) rest = xi2v;
                 return rest;// * chi;
	         }, 1e-3, 1e-9);
      
      return res.value  * rho_m;//res_ds.value*2.0/R/R - res.value;
    }
  };
}

#endif
