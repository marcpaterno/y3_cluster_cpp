#ifndef Y3_CLUSTER_ZO_ZT_DES_T_HH
#define Y3_CLUSTER_ZO_ZT_DES_T_HH

#include "utils/primitives.hh"
#include "models/sigma_photoz_des.hh"
// #include <cmath>

namespace y3_cluster {

    class ZO_ZT_DES_t {
    public:
        explicit ZO_ZT_DES_t() {}

        double
        operator()(double zo, double zt) const {
            double _sigma_zt = _sigma_photoz_des(zt);
            return y3_cluster::gaussian(zo, zt, _sigma_zt);
        }

    private:
        y3_cluster::SIGMA_PHOTOZ_DES_t _sigma_photoz_des;
    };
}

#endif
