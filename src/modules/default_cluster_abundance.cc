#include "models/default_models.hh"
#include "utils/generate_module.hh"

// If you want to use something other than default models:
// #include "really_great_MOR.hh"
// struct Y3Models : y3_cluster::DefaultModels {
//     using MOR = y3_cluster::really_great_MOR;
//     }
// using mymodule=y3_cluster::ClustersModule<Y3Models>;
// GENERATE_Y3_COSMOSIS_MODULE(mymodule)

using module = y3_cluster::ClustersModule<y3_cluster::DefaultModels>;
GENERATE_Y3_COSMOSIS_MODULE(module)
