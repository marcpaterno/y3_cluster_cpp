#include "interp_1d.hh"

#include "gsl/gsl_interp.h"
#include <cassert>
#include <stdexcept>

namespace {
  inline
  void free_if_nonnull(void* p) {
    if (p == nullptr) return;
    gsl_interp* interp = static_cast<gsl_interp*>(p);
    gsl_interp_free(interp);
  }
}

y3_cluster::Interp1D::Interp1D(std::vector<double>&& xs,
                               std::vector<double>&& ys)
  : xs_(std::move(xs)), ys_(std::move(ys)), interp_(nullptr)
{
  if (xs_.size() != ys_.size()) {
    throw std::logic_error("vector length mismatch in construction of Interp1D");
  }
  init_gsl_bits_();
}

y3_cluster::Interp1D::Interp1D(Interp1D const& other)
  : xs_(other.xs_), ys_(other.ys_), interp_(nullptr)
{
  init_gsl_bits_();
}

y3_cluster::Interp1D&
y3_cluster::Interp1D::operator=(Interp1D const& rhs) {
  Interp1D tmp(rhs);
  swap(tmp);
  return *this;
}

y3_cluster::Interp1D::Interp1D(Interp1D&& other) noexcept
  : xs_(std::move(other.xs_)),
    ys_(std::move(other.ys_)),
    interp_(other.interp_)
{
  other.interp_ = nullptr;
}

y3_cluster::Interp1D&
y3_cluster::Interp1D::operator=(Interp1D&& rhs) noexcept {
  if (this == &rhs) return *this;
  free_if_nonnull(interp_);
  xs_ = std::move(rhs.xs_);
  ys_ = std::move(rhs.ys_);
  interp_ = rhs.interp_;
  rhs.interp_ = nullptr;
  return *this;
}


y3_cluster::Interp1D::~Interp1D() noexcept
{
  free_if_nonnull(interp_);
}

double
y3_cluster::Interp1D::operator()(double x) const
{
  assert(interp_); // make sure we haven't been move'd from.
  double res;
  int rc = gsl_interp_eval_e(static_cast<gsl_interp*>(interp_),
      xs_.data(),
      ys_.data(),
      x,
      nullptr,
      &res);
  if (rc == 0)
    return res;
  throw std::domain_error("argument out of range in Interp1D");
}

void
y3_cluster::Interp1D::init_gsl_bits_() {
  gsl_interp* temp =  gsl_interp_alloc(gsl_interp_linear, xs_.size());
  gsl_interp_init(temp, xs_.data(), ys_.data(), xs_.size());
  interp_ = temp;
  assert(interp_);
}

