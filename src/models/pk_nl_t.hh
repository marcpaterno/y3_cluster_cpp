#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "cosmosis/datablock/section_names.h"
#include "utils/datablock_reader.hh"
#include "utils/interp_2d.hh"
#include <memory>
#include <math.h>
using namespace y3_cluster;
#define PI 3.14159265

class pk_nl{
public:
  pk_nl(std::shared_ptr<Interp2D const> pkz)
      : _pkz(pkz) {}

  explicit pk_nl(cosmosis::DataBlock& sample)
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
