#ifndef Y3_CLUSTER_SILLY_MODELS_HH
#define Y3_CLUSTER_SILLY_MODELS_HH

#include "default_models.hh"

// This is an example of how to define alternate models for Cython. The
// cython build script will see this file `silly_models.hh`, and because of the
// `_models` suffix, expect it to contain a type called `SillyModels`.
// It then generates code to build a python interface for
// `Gamma_T_Integrand<SillyModels>`.
namespace y3_cluster {
    struct SillyModels : DefaultModels {
    };
}

#endif
