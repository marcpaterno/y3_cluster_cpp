#include <datablock_reader.hh>
#include "interp_2d.hh"
#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cosmosis/datablock/ndarray.hh"
#include "/cosmosis/cosmosis/datablock/section_names.h"
#include <memory>
#include <math.h>
using namespace y3_cluster;
#define PI 3.14159265

class xi_nl{
public:
  xi_nl(std::shared_ptr<Interp2D const> pkz)
      : _pkz(pkz) {}

  explicit xi_nl(cosmosis::DataBlock& sample)
	   : _pkz(std::make_shared<Interp2D const>(
		 get_datablock<std::vector<double>>(sample, "matter_power_nl", "k_h"),
		 get_datablock<std::vector<double>>(sample, "matter_power_nl", "z"),
		 get_datablock<cosmosis::ndarray<double>>(sample, "matter_power_nl", "p_k")) )
	    {}

  double
  operator()(double r, double z) const
  {
      double tmp = 0.0;
      for (std::size_t i = 0; i < 1000; i++)
      {  
         double k = 0.01/r * exp(i*0.01);
         double x = k * r;
         double w = sin(x) / x;
         tmp = tmp + w * k*k * _pkz->eval(k, z) / (2.0 * PI);
       }
       return tmp;
    }

private:
 std::shared_ptr<Interp2D const> _pkz;
};
