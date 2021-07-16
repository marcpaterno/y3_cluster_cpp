#ifndef Y3_CLUSTER_CPP_CUDA_MODULE_MACROS_HH
#define Y3_CLUSTER_CPP_CUDA_MODULE_MACROS_HH
////////////////////////////////////////////////////////////////////////
//
// Defines the macros
// DEFINE_COSMOSIS_CUDA_INTEGRATION_MODULE(<module_classname>) to be used
// in XXX_module.cc to declare a CosmoSIS module that integrates a scalar-
// valued integrand using CUDA.
//
// Note: Libraries that include these symbol definitions cannot be
// linked into a main program as other libraries are.  This is because
// the "one definition" rule would be violated.
//
////////////////////////////////////////////////////////////////////////

#include "CosmoSISSICUDAModule.hh"

// Produce the injected functions

#define DEFINE_COSMOSIS_CUDA_INTEGRATION_MODULE(klass)                         \
  using module_type = y3_cluster::CosmoSISSICUDAModule<klass>;                 \
  extern "C" {                                                                 \
  void*                                                                        \
  setup(cosmosis::DataBlock* cfg)                                              \
  {                                                                            \
    return new module_type(*cfg);                                              \
  }                                                                            \
                                                                               \
  DATABLOCK_STATUS                                                             \
  execute(cosmosis::DataBlock* sample, void* module)                           \
  {                                                                            \
    auto mod = static_cast<module_type*>(module);                              \
    mod->execute(*sample);                                                     \
    return DBS_SUCCESS;                                                        \
  }                                                                            \
                                                                               \
  int                                                                          \
  cleanup(void* module)                                                        \
  {                                                                            \
    delete static_cast<module_type*>(module);                                  \
    return 0;                                                                  \
  }                                                                            \
  }

#endif
