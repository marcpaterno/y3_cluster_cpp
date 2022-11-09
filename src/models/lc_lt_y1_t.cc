#include "lc_lt_y1_t.hh"
#include "utils/read_vector.hh"

#include <string>

// apparently reading in a table from a file

namespace {
  inline std::string
  augment_filename(char const* part)
  {
    return std::string("lc_lt_fits/DESY1A_v1.4/") + part;
  }

  inline std::vector<double>
  read_vector_from_fileroot(char const* part)
  {
    return read_vector(augment_filename(part));
  }

  auto const lamtrue_in = read_vector(augment_filename("lamtrue.txt"));
  auto const z_in = read_vector(augment_filename("z.txt"));
}

y3_cluster::LC_LT_Y1_t::LC_LT_Y1_t()
  : tau_interp(lamtrue_in, z_in, read_vector_from_fileroot("tau.txt"))
  , mu_interp(lamtrue_in, z_in, read_vector_from_fileroot("mu.txt"))
  , sigma_interp(lamtrue_in, z_in, read_vector_from_fileroot("sigma.txt"))
  , fmsk_interp(lamtrue_in, z_in, read_vector_from_fileroot("fmsk.txt"))
  , fprj_interp(lamtrue_in, z_in, read_vector_from_fileroot("fprj.txt"))
{}
