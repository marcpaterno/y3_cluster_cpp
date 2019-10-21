include_directories($ENV{COSMOSIS_SRC_DIR} ${CUBA_INCLUDE_DIR}) 
link_directories($ENV{COSMOSIS_SRC_DIR}/cosmosis/datablock)

add_library(example_vector_integrand SHARED ExampleVectorIntegrand.cc)
target_link_libraries(example_vector_integrand PRIVATE cosmosis)
target_include_directories(example_vector_integrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_library(example_scalar_integrand SHARED ExampleScalarIntegrand.cc)
target_link_libraries(example_scalar_integrand PRIVATE cosmosis)
target_include_directories(example_scalar_integrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

# Note
# The following two "modules" are abnormal, in that they are declared as a
# shared libraries rather than as modules. This is so they can be used in testing the
# module templates themselves. This should not be done for actual modules, which should
# be declared as modules.
add_library(ExampleVectorIntegrationModule SHARED ExampleVectorIntegration_module.cc)
target_link_libraries(ExampleVectorIntegrationModule
                      PRIVATE
                      cosmosis
                      example_vector_integrand
                      ${CUBA_LIBRARIES})
target_include_directories(ExampleVectorIntegrationModule
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_library(ExampleScalarIntegrationModule SHARED ExampleScalarIntegration_module.cc)
target_link_libraries(ExampleScalarIntegrationModule
                      PRIVATE
                      cosmosis
                      example_scalar_integrand
                      ${CUBA_LIBRARIES})
target_include_directories(ExampleScalarIntegrationModule
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_executable(test_ExampleVectorIntegrationModule test_ExampleVectorIntegrationModule.cc)
target_link_libraries(test_ExampleVectorIntegrationModule
                      PRIVATE
                      cosmosis
                      example_vector_integrand
                      ExampleVectorIntegrationModule
                      example_vector_integrand)
target_include_directories(test_ExampleVectorIntegrationModule
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src)
add_test(test_ExampleVectorIntegrationModule test_ExampleVectorIntegrationModule)

add_executable(test_ExampleScalarIntegrationModule test_ExampleScalarIntegrationModule.cc)
target_link_libraries(test_ExampleScalarIntegrationModule
                      PRIVATE
                      cosmosis
                      example_scalar_integrand
                      ExampleScalarIntegrationModule
                      example_scalar_integrand)
target_include_directories(test_ExampleScalarIntegrationModule
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src)
add_test(test_ExampleScalarIntegrationModule test_ExampleScalarIntegrationModule)




add_executable(ExampleVectorIntegrand_test ExampleVectorIntegrand.test.cc)
target_link_libraries(ExampleVectorIntegrand_test
                      PRIVATE
                      example_vector_integrand
                      cosmosis
                      ${CUBA_LIBRARIES}
                      testmain)
target_include_directories(ExampleVectorIntegrand_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR}
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(ExampleVectorIntegrand_test ExampleVectorIntegrand_test)

add_executable(ExampleScalarIntegrand_test ExampleScalarIntegrand.test.cc)
target_link_libraries(ExampleScalarIntegrand_test
                      PRIVATE
                      example_scalar_integrand
                      cosmosis
                      ${CUBA_LIBRARIES}
                      testmain)
target_include_directories(ExampleScalarIntegrand_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR}
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(ExampleScalarIntegrand_test ExampleScalarIntegrand_test)


add_library(snapshotsim_scalarintegrand MODULE Snapshotsim_ScalarIntegrand.cc)
target_link_libraries(snapshotsim_scalarintegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(snapshotsim_scalarintegrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})

add_library(sigma_nl MODULE sigma_nl.cc)
target_link_libraries(sigma_nl PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(sigma_nl
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})


add_library(model_sigmahm MODULE model_sigmahm.cc)
target_link_libraries(model_sigmahm PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(model_sigmahm
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})

add_library(snapshotsim_scalarintegrand_NC MODULE Snapshotsim_ScalarIntegrand_NC.cc)
target_link_libraries(snapshotsim_scalarintegrand_NC PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(snapshotsim_scalarintegrand_NC
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})


add_library(y1_analysis MODULE y1_analysis.cc)
target_link_libraries(y1_analysis PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(y1_analysis
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})

add_library(sdss_analysis MODULE sdss_analysis.cc)
target_link_libraries(sdss_analysis PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(sdss_analysis
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})

#add_library(default_cluster_abundance MODULE default_cluster_abundance.cc)
#target_link_libraries(default_cluster_abundance
#                      PRIVATE
#                      utils
#                      profiler
#                      cfitsio
#                      ${CUBA_LIBRARIES}
#                      ${GSL_LIBRARIES}
#                      cosmosis)
#target_include_directories(default_cluster_abundance
#                           PRIVATE
#                           ${CMAKE_SOURCE_DIR}
#                           ${CMAKE_SOURCE_DIR}/externals
#                           ${CMAKE_SOURCE_DIR}/src
#                           ${CUBA_INCLUDE_DIR}
#                           ${CUBACPP_DIR}
#                           ${EIGEN3_INCLUDE_DIR})

add_library(cluster_abundance_covariance MODULE cluster_abundance_covariance.cc)
target_link_libraries(cluster_abundance_covariance
                      PRIVATE
                      utils
                      models
                      ${CUBA_LIBRARIES}
                      ${GSL_LIBRARIES}
                      cosmosis)
target_include_directories(cluster_abundance_covariance
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CMAKE_SOURCE_DIR}/src
                           ${GSL_INCLUDE_DIRS}
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR})

add_library(ExampleVectorIntegration MODULE ExampleVectorIntegration_module.cc)
target_link_libraries(ExampleVectorIntegration
                     PRIVATE
                     cosmosis
                     example_vector_integrand
                     ${CUBA_LIBRARIES})                     

target_include_directories(ExampleVectorIntegration
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_library(sigma_miscent MODULE sigma_miscent.cc)
target_link_libraries(sigma_miscent PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(sigma_miscent
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})
 
add_library(sigma_cent MODULE sigma_cent.cc)
target_link_libraries(sigma_cent PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(sigma_cent
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})

add_library(buzzard_sigma_halos MODULE buzzard_sigma_halos.cc)
target_link_libraries(buzzard_sigma_halos PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(buzzard_sigma_halos
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})


add_library(sigma_cent_y1 MODULE sigma_cent_y1.cc)
target_link_libraries(sigma_cent_y1 PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(sigma_cent_y1
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})


add_library(sigma_miscent_y1 MODULE sigma_miscent_y1.cc)
target_link_libraries(sigma_miscent_y1 PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(sigma_miscent_y1
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})

add_library(kappa_cent_y1 MODULE kappa_cent_y1.cc)
target_link_libraries(kappa_cent_y1 PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(kappa_cent_y1
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})

add_library(kappa_miscent_y1 MODULE kappa_miscent_y1.cc)
target_link_libraries(kappa_miscent_y1 PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(kappa_miscent_y1
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})

add_library(compton_y_sims MODULE compton_y_sims.cc)
target_link_libraries(compton_y_sims PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(compton_y_sims
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})

add_library(sptxdes_abundance_integral sptxdes/abundance/abundance_integral_module.cc)
target_link_libraries(sptxdes_abundance_integral PRIVATE cosmosis utils models
                      ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(sptxdes_abundance_integral PRIVATE ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
