#ifndef MATTEOS_ZO_ZT_HH
#define MATTEOS_ZO_ZT_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"

#include <array>
#include <cmath>
#include <polynomial.hh>

namespace y3_cluster {
  namespace {
    static const std::array<polynomial<4>, 5>
        sigmas{{std::array<double, 4>{-1.18159413, 1.10688, -0.24906221, 0.021577},
                std::array<double, 4>{-1.22925, 1.117566, -0.25085, 0.02129},
                std::array<double, 4>{-1.26122355, 1.12986624, -0.25394517, 0.0212711},
                std::array<double, 4>{-1.16167235, 1.05896902, -0.23973177, 0.02015095},
                std::array<double, 4>{-1.20483017, 1.0778892, -0.24311016, 0.02009036}}};
  }

  class MATTEOS_ZO_ZT_t {
  public:
    MATTEOS_ZO_ZT_t() {}
    explicit MATTEOS_ZO_ZT_t(cosmosis::DataBlock&) {}

    double
    operator()(double zomin, double zomax, double zt) const
    {
      double sigma;

      /*
      if (lo < 27.9)
          sigma = sigmas[0](zt);
      else if (lo < 37.6)
          sigma = sigmas[1](zt);
      else if (lo < 50.3)
          sigma = sigmas[2](zt);
      else if (lo < 69.3)
          sigma = sigmas[3](zt);
      else
      */
          sigma = sigmas[4](zt);

      // P(zo | zt) := (1.0 / sqrt(2pi) / sigma) * exp(- (zo - zt) * (zo - zt) / (2 * sigma * sigma))
      //    (i.e., a standard gaussian)
      // So,
      //    \int P(zo|zt) d(zo), zo in [zomin, zomax]
      //     == (erf((zomax - zt) / base) - erf((zomin - zt) / base)) / 2

      // Using matteo's sigma equation
      double base = std::sqrt(2) * sigma * 1.014;
      return (std::erf((zomax - (zt + 0.002)) / base) -
              std::erf((zomin - (zt + 0.002)) / base)) / 2.0;
    }
  };
}

#endif
