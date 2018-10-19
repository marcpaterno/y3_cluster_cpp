#ifndef Y3_CLUSTER_MATTEOS_LC_LT_HH
#define Y3_CLUSTER_MATTEOS_LC_LT_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "read_vector.hh"

#include <iostream>
#include <sstream>
#include <string>

namespace y3_cluster {
  namespace {
    auto lts = read_vector("int_P_lob_ltr/l_in_table_v27_0.dat");
    auto zts = read_vector("int_P_lob_ltr/z_in_table_v27_0.dat");

    Interp2D
    make_Interp2D(std::size_t i)
    {
      std::ostringstream filename;
      filename << "int_P_lob_ltr/int_P_lob_ltr_ztr_Deltal_" << i
               << std::string("_v27_0.dat");
      const auto lc_lt_integrated = read_vector(std::string(filename.str()).c_str());
      if (lts.size() * zts.size() != lc_lt_integrated.size()) {
        lts.pop_back();
        zts.pop_back();
      }
      return Interp2D(lts, zts, lc_lt_integrated);
    }
  }

  struct MATTEOS_LC_LT_t {
    static std::array<Interp2D, 5> const lc_lt_tables;

    explicit MATTEOS_LC_LT_t(const cosmosis::DataBlock&) {}
    MATTEOS_LC_LT_t() {}

    double
    operator()(std::size_t loi, double lt, double zt) const
    {
      if (loi >= 5)
        throw std::runtime_error("Richness bin out of bounds");
      return lc_lt_tables[loi].eval(lt, zt);
    }
  };
}

const std::array<y3_cluster::Interp2D, 5>
y3_cluster::MATTEOS_LC_LT_t::lc_lt_tables = {{make_Interp2D(1),
                                            make_Interp2D(2),
                                            make_Interp2D(3),
                                            make_Interp2D(4),
                                            make_Interp2D(5)}};

#endif
