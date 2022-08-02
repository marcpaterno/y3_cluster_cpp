#include "int_lc_lt_des_t.cuh"
#include "models/interpolation_tables.cuh"

using namespace quad;

namespace {
  // make_short_vec creates an std::vector<double> from an std::array,
  // using all but the last element of the std::array.
  template <size_t N>
  std::vector<double>
  make_short_vec(std::array<double, N> const& a)
  {
    static_assert(N != 0, "make_short_vec requires a nonzero-length array");
    return {a.begin(), a.end() - 1};
  }

  // make_vec creates an std::vector<double> from an std::array, using all the
  // values in the std::array.
  template <size_t N>
  std::vector<double>
  make_vec(std::array<double, N> const& a)
  {
    return {a.begin(), a.end()};
  }

  // Create an Interp2D from an x-axis, y-axis, and z "matrix", with the matrix
  // unrolled into a one-dimenstional array.
  template <size_t M, std::size_t N>
  Interp2D
  make_Interp2D_aux(std::array<double, M> const& xs,
                    std::array<double, N> const& ys,
                    std::array<double, (N) * (M)> const& zs)
  {
    return {xs, ys, zs};
  }

  auto make_Interp2D = [](auto zs) {
    return make_Interp2D_aux(lt_bins, zt_bins, zs);
  };
}
