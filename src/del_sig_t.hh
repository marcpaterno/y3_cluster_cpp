#ifndef Y3_CLUSTER_DEL_SIG_T_HH
#define Y3_CLUSTER_DEL_SIG_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "ez.hh"
#include "primitives.hh"
#include "interp_2d.hh"

#include <memory>

#include <cmath>

// namespace y3_cluster
// {
//   class DEL_SIG_CEN_t {
//   public:
//     explicit DEL_SIG_CEN_t(double c) : _c(c) {}

//     explicit DEL_SIG_CEN_t(cosmosis::DataBlock& sample)
//     {
//       sample.get_val<double>("del_sig_cen_params", "c", _c);
//     }

//     double
//     operator()(double r, double lnM, double zt) const
//     /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
//     {
//       y3_cluster::EZ_sq ez_sq{0.3, 0.7, 0.};

//       double delta_c =
//         200. * _c * _c * _c / (3. * (std::log(1. + _c) - _c / (1. + _c)));
//       double rho_crit = 2.77526157E11 * ez_sq(zt);
//       /* EZ*EZ would be a little bit slower than direct definition */

//       double r_200 = std::pow(
//         3. * std::exp(lnM) / (800. * y3_cluster::pi() * rho_crit), 1. / 3.);
//       double r_s = r_200 / _c;

//       double r_ratio = r / r_s;
//       double coeff = r_s * delta_c * rho_crit;

//       if (r_ratio < 1.) {
//         return coeff *
//                (8. * std::atanh(std::sqrt((1. - r_ratio) / (1. + r_ratio))) /
//                   (r_ratio * r_ratio * std::sqrt(1. - r_ratio * r_ratio)) +
//                 4. * std::log(r_ratio / 2.) / (r_ratio * r_ratio) -
//                 2. / (r_ratio * r_ratio - 1.) +
//                 4. * std::atanh(std::sqrt((1. - r_ratio) / (1. + r_ratio))) /
//                   ((r_ratio * r_ratio - 1.) *
//                    std::sqrt(1. - r_ratio * r_ratio)));
//       } else if (r_ratio == 1.) {
//         return coeff * (10. / 3. + 4. * std::log(0.5));
//       } else {
//         return coeff *
//                (8. * std::atan(std::sqrt((r_ratio - 1.) / (r_ratio + 1.))) /
//                   (r_ratio * r_ratio * std::sqrt(r_ratio * r_ratio - 1.)) +
//                 4. * std::log(r_ratio / 2.) / (r_ratio * r_ratio) -
//                 2. / (r_ratio * r_ratio - 1.) +
//                 4. * std::atan(std::sqrt((r_ratio - 1.) / (r_ratio + 1.))) /
//                   (std::pow(r_ratio * r_ratio - 1., 1.5)));
//       }
//     }

//   private:
//     double _c;
//   };
// }

namespace y3_cluster
{
  class DEL_SIG_t {
  public:
    DEL_SIG_t(std::shared_ptr<Interp2D const> dsigma1, 
                  std::shared_ptr<Interp2D const> dsigma2, 
                  std::shared_ptr<Interp2D const> bias
                  /*,double c*/) 
                  : _dsigma1(dsigma1), _dsigma2(dsigma2), _bias(bias)/*, _c(c)*/ {}

    using doubles = std::vector<double>;

    explicit DEL_SIG_t(cosmosis::DataBlock& sample)
      : _dsigma1(std::make_shared<Interp2D const>(
          sample.view<doubles>("del_sig_params", "x1"),
          sample.view<doubles>("del_sig_params", "y1"),
          sample.view<doubles>("del_sig_params", "z1")))
      , _dsigma2(std::make_shared<Interp2D const>(
          sample.view<doubles>("del_sig_params", "x2"),
          sample.view<doubles>("del_sig_params", "y2"),
          sample.view<doubles>("del_sig_params", "z2")))
      , _bias(std::make_shared<Interp2D const>(
          sample.view<doubles>("del_sig_params", "x3"),
          sample.view<doubles>("del_sig_params", "y3"),
          sample.view<doubles>("del_sig_params", "z3")))
      // , _c(sample.view<double>("del_sig_params", "c"))
    {}

    double
    operator()(double r, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    { 
      double del_sig_1 = _dsigma1->eval(r,lnM);
      double del_sig_2 = _bias->eval(zt,lnM) * _dsigma2->eval(r,zt);
      if (del_sig_1 >= del_sig_2) {
        return (1.+zt)*(1.+zt)*(1.+zt)*del_sig_1;
      }else{
        return (1.+zt)*(1.+zt)*(1.+zt)*del_sig_2;
      }
    }

  private:
    std::shared_ptr<Interp2D const> _dsigma1;
    std::shared_ptr<Interp2D const> _dsigma2;
    std::shared_ptr<Interp2D const> _bias;
    /*double _c;*/
  };
}

#endif
