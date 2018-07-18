#ifndef Y3_CLUSTER_INTEGRATION_RANGE_HH
#define Y3_CLUSTER_INTEGRATION_RANGE_HH

#include <iostream>
#include <stdexcept>
#include <string>
#include <datablock_reader.hh>

namespace y3_cluster {
  class IntegrationRange {
  public:
    IntegrationRange(double a, double b) : _a(a), _range(b - a)
    {
      if (_range == 0.0)
        throw std::logic_error("zero-length IntegrationRange");
    }

    IntegrationRange(cosmosis::DataBlock& sample, std::string var) {
      double b;
      std::string min = var + "_min";
      std::string max = var + "_max";
      _a = get_datablock<double>(sample, "cluster_abundance", min.c_str());
      b = get_datablock<double>(sample, "cluster_abundance", max.c_str());
      _range = b - _a;
      if (_range == 0.0)
        throw std::logic_error("zero-length IntegrationRange");
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

    friend std::ostream& operator<<(std::ostream& os, const IntegrationRange&);
  private:
    double _a;
    double _range;
  };

  std::ostream&
  operator<<(std::ostream& os, const y3_cluster::IntegrationRange& ir)
  {
      return os << "[" << ir._a << ", " << ir._a + ir._range << "]";
  }
}

#endif
