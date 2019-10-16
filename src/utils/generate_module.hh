#ifndef __Y3_CLUSTER_GENERATE_COSMOSIS_MODULE_HH
// Nobody should need to touch this in a normal circumcstance.
// you shouldn't need to fiddle this if you are trying to set up another cluster
// cosmology implementation
#define __Y3_CLUSTER_GENERATE_COSMOSIS_MODULE_HH

#include "clusters_module.hh"
#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/entry.hh"
#include "cosmosis/datablock/section.hh"
#include "cosmosis/datablock/section_names.h"

#define GENERATE_Y3_COSMOSIS_MODULE(klass)                                     \
                                                                               \
  extern "C" {                                                                 \
                                                                               \
  void*                                                                        \
  setup(cosmosis::DataBlock* options)                                          \
  {                                                                            \
    return new klass(*options);                                                \
  }                                                                            \
                                                                               \
  DATABLOCK_STATUS                                                             \
  execute(cosmosis::DataBlock* block, void* config)                            \
  {                                                                            \
    auto mod = static_cast<klass*>(config);                                    \
    mod->execute(*block);                                                      \
                                                                               \
    DATABLOCK_STATUS status = DBS_SUCCESS;                                     \
                                                                               \
    return status;                                                             \
  }                                                                            \
                                                                               \
  int                                                                          \
  cleanup(void* config)                                                        \
  {                                                                            \
    delete static_cast<klass*>(config);                                        \
    return 0;                                                                  \
  }                                                                            \
  }

#endif
