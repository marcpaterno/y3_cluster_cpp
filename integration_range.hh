#ifndef Y3_CLUSTER_INTEGRATION_RANGE_HH
#define Y3_CLUSTER_INTEGRATION_RANGE_HH

namespace y3_cluster {
  class IntegrationRange {
  public:
    IntegrationRange(double a, double b) : _a(a), _range(b - a) {}

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