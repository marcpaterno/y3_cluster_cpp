#include "catch2/catch.hpp"
#include "lc_lt_t.hh"

#include <fstream>
#include <iostream>
#include <string>

using y3_cluster::LC_LT_t;
TEST_CASE("Lc_Lt_t works")
{

  std::vector<double> ys(195);
  std::vector<double> zs(10);

  LC_LT_t lc_lt;
  for (std::size_t j = 0; j != zs.size(); ++j) {
 
      double sum = 0;
      zs[j]=j*0.05+0.1;
      for (std::size_t i = 0, sy = ys.size(); i != sy; ++i) {
         ys[i]=i+1.0;
   	 double const fz = lc_lt(5, ys[i], zs[j]);
	 sum=sum+fz;
      }
      double constexpr epsrel = 1.0e-6;
      double constexpr epsabs = 1.0e-12;
      CHECK(sum == Approx(1.0).epsilon(epsrel).margin(epsabs));

  }
}
