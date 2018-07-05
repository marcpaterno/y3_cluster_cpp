#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cosmosis/datablock/section_names.h"
#include "test/clusters_module.hh"
#define GENERATE_Y3_COSMOSIS_MODULE(models) \
\
extern "C" { \
\
void * setup(cosmosis::DataBlock * options) \
{ \
    return new y3_cluster::ClustersModule<y3_cluster::DefaultModels>(*options); \
} \
\
DATABLOCK_STATUS execute(cosmosis::DataBlock * block, void * config) \
{ \
    auto mod = static_cast<y3_cluster::ClustersModule<y3_cluster::DefaultModels> *>(config); \
    mod->execute(*block); \
\
    DATABLOCK_STATUS status = DBS_SUCCESS; \
\
    return status; \
} \
\
int cleanup(void * config) \
{ \
    delete static_cast<y3_cluster::ClustersModule<y3_cluster::DefaultModels> *>(config); \
    return 0;\
} \
\
}
