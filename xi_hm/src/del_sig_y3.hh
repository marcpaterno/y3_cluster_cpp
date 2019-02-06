#ifndef Y3_CLUSTER_DEL_SIG_T_HH
#define Y3_CLUSTER_DEL_SIG_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "ez.hh"
#include "primitives.hh"
#include "interp_2d.hh"
#include <datablock_reader.hh>
#include <cubacpp/cubacpp.hh>
#include <math.h>
#include <memory>
#include <cmath>
#include <iostream>

namespace y3_cluster
{
  class DEL_SIG_y3 {
  private:
    std::shared_ptr<Interp2D const> _xi1;
    std::shared_ptr<Interp2D const> _xi2;
    std::shared_ptr<Interp2D const> _bias;

  public:
    DEL_SIG_y3(std::shared_ptr<Interp2D const> xi1, 
              std::shared_ptr<Interp2D const> xi2, 
              std::shared_ptr<Interp2D const> bias)
              : _xi1(xi1), _xi2(xi2), _bias(bias) {}

    using doubles = std::vector<double>;

    // TODO: This needs to be reading cosmosis datablock parameters
    explicit DEL_SIG_y3(cosmosis::DataBlock& sample)
      : _xi1(std::make_shared<Interp2D const>(
          get_datablock<doubles>(sample, "cluster_abundance", "R_perp"),
          get_datablock<doubles>(sample, "cluster_abundance", "lnM"),
          get_datablock<cosmosis::ndarray<double>>(sample, "cluster_abundance", "xi_nfw")))
      , _xi2(std::make_shared<Interp2D const>(
          get_datablock<doubles>(sample, "cluster_abundance", "R_perp"),
          get_datablock<doubles>(sample, "cluster_abundance", "z"),
          get_datablock<cosmosis::ndarray<double>>(sample, "cluster_abundance", "xi_nl")))
      , _bias(std::make_shared<Interp2D const>(
          get_datablock<doubles>(sample, "tinker_bias_function", "ln_mass_h"),
          get_datablock<doubles>(sample, "tinker_bias_function", "z"),
          get_datablock<cosmosis::ndarray<double>>(sample, "tinker_bias_function", "bias")))
    {}

    double
    operator()(double R, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    { 
      cubacpp::QAG integrator(-20.0, 20.0, GSL_INTEG_GAUSS61, 50000);
      /* const auto res = integrator.integrate([=](double chi) {
		 double r=sqrt(R*R+chi*chi);
                 double xi1v = _xi1->eval(r,lnM);
                 double xi2v = _bias->eval(lnM, zt) * _xi2->eval(r,zt);		 
                 double rest = xi1v;
                 if (xi1v >= xi2v) rest = xi2v;
                 return rest;
	         }, 1e-3, 1e-9);
      
      cubacpp::Cuhre integrator_ds;
      //this part need to set up jacobian
     const auto res_ds = integrator_ds.integrate([=](double chi, double r) {
			double dd=sqrt(r*r+chi*chi);
			double xi1v = _xi1->eval(dd,lnM);
			double xi2v = _bias->eval(lnM, zt) * _xi2->eval(dd,zt);
		        double rest = xi1v;
	                if (xi1v >= xi2v) rest = xi2v;
                        return rest;
			}, 1e-5, 1e-9);
      */
      return 0;//res_ds.value*2.0/R/R - res.value;
    }
  };
}

#endif
