#include "utils/CosmoSISVectorIntegrationModule.hh"
#include "ExampleVectorIntegrand.hh"

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cosmosis/datablock/entry.hh"
#include "/cosmosis/cosmosis/datablock/section.hh"
#include "/cosmosis/cosmosis/datablock/section_names.h"

using ExampleVectorIntegrationModule =
  y3_cluster::CosmoSISVectorIntegrationModule<ExampleVectorIntegrand>;

extern "C"
{
  void*
  setup(cosmosis::DataBlock* config)
  {
    return new ExampleVectorIntegrationModule(*config);
  }

  DATABLOCK_STATUS
  execute(cosmosis::DataBlock* sample, void* module)
  {
    auto mod = static_cast<ExampleVectorIntegrationModule*>(module);
    mod->execute(*sample);
    return DBS_SUCCESS;
  }

  int
  cleanup(void* module)
  {
    delete static_cast<ExampleVectorIntegrationModule*>(module);
    return 0;
  }
}
