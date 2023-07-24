#ifndef Y3_CLUSTER_INTERP_1D_HH
#define Y3_CLUSTER_INTERP_1D_HH

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

// Interp1D is used for linear interpolation in 1 dimension.
// It uses the GSL library to do the actual interpolation.
// Interp1D objects do not allow extrapolation.

namespace y3_cluster {
  using gsl_interp_ptr = void*;

  class Interp1D {
  public:
    // Interpolator created from two arrays; compiler assures they are of the
    // same length.
    template <std::size_t N>
    Interp1D(std::array<double, N> const& xs, std::array<double, N> const& ys);

    // Interpolator created from two vectors: throws std::logic_error if the
    // vectors are not of the same length.
    Interp1D(std::vector<double>&& xs, std::vector<double>&& ys);

    // As above, but deep-copy the vectors instead of just moving them. */
    Interp1D(std::vector<double> const& xs, std::vector<double> const& ys)
      : Interp1D{std::vector<double>(xs), std::vector<double>(ys)}
    {}

    Interp1D(Interp1D const& other);
    Interp1D& operator=(Interp1D const& rhs);
    Interp1D(Interp1D&& other) noexcept;
    Interp1D& operator=(Interp1D&& rhs) noexcept;

    // Destructor must clean up allocated GSL resources.
    ~Interp1D() noexcept;

    void swap(Interp1D& other) noexcept;

    double operator()(double x) const;

    double
    eval(double x) const
    {
      return this->operator()(x);
    };

    double
    clamp(double x) const
    {
      double const xx = std::clamp(x, xs_.front(), xs_.back());
      return (*this)(xx);
    }

  private:
    std::vector<double> xs_;
    std::vector<double> ys_;
    gsl_interp_ptr interp_;

    void init_gsl_bits_();
  };
} // namespace y3_cluster

template <std::size_t N>
inline y3_cluster::Interp1D::Interp1D(std::array<double, N> const& xs,
                                      std::array<double, N> const& ys)
  : Interp1D({begin(xs), end(xs)}, {begin(ys), end(ys)})
{}

inline void
y3_cluster::Interp1D::swap(Interp1D& other) noexcept
{
  using std::swap;
  swap(xs_, other.xs_);
  swap(ys_, other.ys_);
  swap(interp_, other.interp_);
}

#endif
