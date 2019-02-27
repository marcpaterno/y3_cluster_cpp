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
    std::shared_ptr<Interp1D const> _xi1;

  public:
    DEL_SIG_y3(std::shared_ptr<Interp1D const> xi1) 
              : _xi1(xi1){}

    using doubles = std::vector<double>;

    // TODO: This needs to be reading cosmosis datablock parameters
    explicit DEL_SIG_y3(cosmosis::DataBlock& sample)
      : _xi1(std::make_shared<y3_cluster::Interp1D const>(
          get_datablock<doubles>(sample, "cluster_abundance", "R_perp"),
          get_datablock<doubles>(sample, "cluster_abundance", "sigma")))
    {}

    double
    operator()(double R, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    { 

       cubacpp::QAG integrator(0.002, R, GSL_INTEG_GAUSS61, 50000000);
       const auto res = integrator.integrate([=](double r) {
                 double xi1v = _xi1->eval(r);
                 double rest = xi1v * r;
                 return rest;// * chi;
	         }, 1e-3, 1e-9);
      
      return res.value *2.0/R/R - _xi1->eval(R);//res_ds.value*2.0/R/R - res.value;
    }
  };
}

#endif
