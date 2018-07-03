#define GENERATE_Y3_COSMOSIS_MODULE(models) \
#include "cosmosis/datablock/datablock.hh" \
#include "cosmosis/datablock/section_names.h" \
#include "my_calculation_code.h" \
\
extern "C" { \
\
void * setup(DataBlock * options) \
{ \
    return new y3_cluster::ClustersModule<(models)>(cosmosis::DataBlock& options); \
} \
\
DATABLOCK_STATUS execute(DataBlock * block, void * config) \
{ \
    auto mod = static_cast<y3_cluster::ClustersModule<(models)>> *>(config); \
    mod->execute(*block); \
\
    DATABLOCK_STATUS status = 0; \
\
    return status; \
} \
\
int cleanup(void * config) \
{ \
    delete static_cast<y3_cluster::ClustersModule<(models)> *>(config); \
} \
\
}
