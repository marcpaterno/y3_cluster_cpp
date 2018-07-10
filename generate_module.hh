#ifndef __Y3_CLUSTER_GENERATE_COSMOSIS_MODULE_HH
#define __Y3_CLUSTER_GENERATE_COSMOSIS_MODULE_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cosmosis/datablock/entry.hh"
#include "/cosmosis/cosmosis/datablock/section_names.h"
#include "/cosmosis/cosmosis/datablock/section.hh"
#include "test/clusters_module.hh"
#define GENERATE_Y3_COSMOSIS_MODULE(models, nradii) \
\
extern "C" { \
\
void * setup(cosmosis::DataBlock * options) \
{ \
    try { \
        return new y3_cluster::ClustersModule<y3_cluster::DefaultModels, nradii>(*options); \
    } catch (cosmosis::DataBlock::BadDataBlockAccess e) { \
        std::cerr << "cosmosis::DataBlock::BadDataBlock exception on setup of y3_cluster_cpp:\n"; \
        throw e; \
    } catch (cosmosis::Section::BadSectionAccess e) { \
        std::cerr << "cosmosis::Section::BadSectionAccess exception on setup of y3_cluster_cpp:\n"; \
        throw e; \
    } catch (cosmosis::Entry::BadEntry e) { \
        std::cerr << "cosmosis::Entry::BadEntry exception on setup of y3_cluster_cpp:\n"; \
        throw e; \
    } \
} \
\
DATABLOCK_STATUS execute(cosmosis::DataBlock * block, void * config) \
{ \
    auto mod = static_cast<y3_cluster::ClustersModule<y3_cluster::DefaultModels, nradii> *>(config); \
    mod->execute(*block); \
\
    DATABLOCK_STATUS status = DBS_SUCCESS; \
\
    return status; \
} \
\
int cleanup(void * config) \
{ \
    delete static_cast<y3_cluster::ClustersModule<y3_cluster::DefaultModels, nradii> *>(config); \
    return 0;\
} \
\
}

#endif
