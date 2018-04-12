#ifndef Y3_CLUSTER_INTEGRATION_RANGE_HH
#define Y3_CLUSTER_INTEGRATION_RANGE_HH

#include <stdexcept>

namespace y3_cluster {
  class IntegrationRange {
  public:
    IntegrationRange(double a, double b) : _a(a), _range(b - a) {
      if (_range == 0.0) throw std::logic_error("zero-length IntegrationRange");
    }

    double
    jacobian() const
    {
      return _range;
    }

    double
    transform(double x) const
    {
      return _range * x + _a;
    }

  private:
    double _a;
    double _range;
  };
}

#endif