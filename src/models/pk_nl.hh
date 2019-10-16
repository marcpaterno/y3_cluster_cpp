#include <datablock_reader.hh>
#include "interp_2d.hh"
#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cosmosis/datablock/ndarray.hh"
#include "/cosmosis/cosmosis/datablock/section_names.h"
#include <memory>
#include <math.h>
using namespace y3_cluster;
#define PI 3.14159265
#include <cubacpp/cubacpp.hh>

class pk_nl{
public:
  pk_nl(std::shared_ptr<Interp2D const> pkz)
      : _pkz(pkz) {}

  explicit xi_nl(cosmosis::DataBlock& sample)
	   : _pkz(std::make_shared<Interp2D const>(
		 get_datablock<std::vector<double>>(sample, "matter_power_nl", "k_h"),
		 get_datablock<std::vector<double>>(sample, "matter_power_nl", "z"),
		 get_datablock<cosmosis::ndarray<double>>(sample, "matter_power_nl", "p_k")) )
	    {}

  double
  operator()(double k, double z) const
  {

      return _pkz->eval(k, z);
    }

private:
 std::shared_ptr<Interp2D const> _pkz;
};
