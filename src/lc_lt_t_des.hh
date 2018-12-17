#ifndef Y3_CLUSTER_LC_LT_T_DES_HH
#define Y3_CLUSTER_LC_LT_T_DES_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "interp_2d.hh"
#include "primitives.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

namespace y3_cluster {
  struct LC_LT_t_des {

    static Interp2D const lambda0_interp;
    static Interp2D const lambda1_interp;
    static Interp2D const lambda2_interp;
    static Interp2D const lambda3_interp;

    explicit LC_LT_t_des(const cosmosis::DataBlock&) {}
    LC_LT_t_des() {}

    double
    operator()(double, double lt, double zt) const
    {
       return lambda3_interp(lt, zt);
    }
  };
}

#endif
