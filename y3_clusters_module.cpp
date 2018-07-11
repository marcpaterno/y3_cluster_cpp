#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/section_names.h"
#include "test/clusters_module.hh"
#include "test/default_models.hh"

using Models = DefaultModels;

extern "C" {

void * setup(DataBlock * options)
{
	return new ClustersModule<Models>(cosmosis::DataBlock& config);
}

DATABLOCK_STATUS execute(DataBlock * block, void * config)
{
    // Config is whatever you returned from setup above
    // Block is the collection of parameters and calculations for
    // this set of cosmological parameters
	
	auto mod = static_cast<ClustersModule<Models> *>(config);
	mod->execute(*block);

    DATABLOCK_STATUS status = 0;

    return status;
}

int cleanup(void * config)
{
	delete static_cast<ClusterModule<Models> *>(config);
}

}
