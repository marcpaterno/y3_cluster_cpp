#include "test/default_models.hh"
#include "generate_module.hh"

// If you want to use something other than default models:
//
// struct Y3Models : y3_cluster::DefaultModels {
//     using MOR = y3_cluster::really_great_MOR,
//     ...
//     }
//
// GENERATE_Y3_COSMOSIS_MODULE(Y3Models)

GENERATE_Y3_COSMOSIS_MODULE(y3_cluster::DefaultModels)

/*
#include "test/clusters_module.hh"

extern "C" { 

void * setup(cosmosis::DataBlock * options) 
{ 
    return new y3_cluster::ClustersModule<y3_cluster::DefaultModels>(*options); 
} 

DATABLOCK_STATUS execute(cosmosis::DataBlock * block, void * config) 
{ 
    auto mod = static_cast<y3_cluster::ClustersModule<y3_cluster::DefaultModels> *>(config); 
    mod->execute(*block); 

    DATABLOCK_STATUS status = DBS_SUCCESS; 

    return status; 
} 

int cleanup(void * config) 
{ 
    delete static_cast<y3_cluster::ClustersModule<y3_cluster::DefaultModels> *>(config); 
    return 0;
} 

}
*/
