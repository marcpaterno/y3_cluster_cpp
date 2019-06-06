#include "utils/CosmoSISIntegrationModule.hh"
#include "ExampleIntegrand.hh"

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cosmosis/datablock/entry.hh"
#include "/cosmosis/cosmosis/datablock/section.hh"
#include "/cosmosis/cosmosis/datablock/section_names.h"

using ExampleIntegrationModule =
  y3_cluster::CosmoSISIntegrationModule<ExampleIntegrand>;

extern "C"
{
  void*
  setup(cosmosis::DataBlock* config)
  {
    return new ExampleIntegrationModule(*config);
  }

  DATABLOCK_STATUS
  execute(cosmosis::DataBlock* sample, void* module)
  {
    auto mod = static_cast<ExampleIntegrationModule*>(module);
    mod->execute(*sample);
    return DBS_SUCCESS;
  }

  int
  cleanup(void* module)
  {
    delete static_cast<ExampleIntegrationModule*>(module);
    return 0;
  }
}
