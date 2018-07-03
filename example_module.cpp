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
