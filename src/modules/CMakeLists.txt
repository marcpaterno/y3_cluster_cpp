
add_library(ExampleVectorIntegrationModule SHARED ExampleVectorIntegration_module.cc)
target_link_libraries(ExampleVectorIntegrationModule
                      PRIVATE
                      /cosmosis/cosmosis/datablock/libcosmosis.so
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
                      /cosmosis/cosmosis/datablock/libcosmosis.so
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
                      /cosmosis/cosmosis/datablock/libcosmosis.so
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
                      /cosmosis/cosmosis/datablock/libcosmosis.so
                      example_scalar_integrand
                      ExampleScalarIntegrationModule
                      example_scalar_integrand)
target_include_directories(test_ExampleScalarIntegrationModule
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src)
add_test(test_ExampleScalarIntegrationModule test_ExampleScalarIntegrationModule)


add_library(example_vector_integrand SHARED ExampleVectorIntegrand.cc)
target_link_libraries(example_vector_integrand PRIVATE /cosmosis/cosmosis/datablock/libcosmosis.so)
target_include_directories(example_vector_integrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_library(example_scalar_integrand SHARED ExampleScalarIntegrand.cc)
target_link_libraries(example_scalar_integrand PRIVATE /cosmosis/cosmosis/datablock/libcosmosis.so)
target_include_directories(example_scalar_integrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_executable(ExampleVectorIntegrand_test ExampleVectorIntegrand.test.cc)
target_link_libraries(ExampleVectorIntegrand_test
                      PRIVATE
                      example_vector_integrand
                      /cosmosis/cosmosis/datablock/libcosmosis.so
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
                      /cosmosis/cosmosis/datablock/libcosmosis.so
                      ${CUBA_LIBRARIES}
                      testmain)
target_include_directories(ExampleScalarIntegrand_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR}
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(ExampleScalarIntegrand_test ExampleScalarIntegrand_test)


add_library(snapshotsim_scalarintegrand SHARED Snapshotsim_ScalarIntegrand.cc)
target_link_libraries(snapshotsim_scalarintegrand PRIVATE /cosmosis/cosmosis/datablock/libcosmosis.so utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(snapshotsim_scalarintegrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_library(snapshotsim_scalarintegrand_NC SHARED Snapshotsim_ScalarIntegrand_NC.cc)
target_link_libraries(snapshotsim_scalarintegrand_NC PRIVATE /cosmosis/cosmosis/datablock/libcosmosis.so utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(snapshotsim_scalarintegrand_NC
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})


add_library(y1_analysis SHARED y1_analysis.cc)
target_link_libraries(y1_analysis PRIVATE /cosmosis/cosmosis/datablock/libcosmosis.so utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(y1_analysis
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_library(sdss_analysis SHARED sdss_analysis.cc)
target_link_libraries(sdss_analysis PRIVATE /cosmosis/cosmosis/datablock/libcosmosis.so utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(sdss_analysis
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

#add_library(default_cluster_abundance SHARED default_cluster_abundance.cc)
#target_link_libraries(default_cluster_abundance
#                      PRIVATE
#                      utils
#                      profiler
#                      cfitsio
#                      ${CUBA_LIBRARIES}
#                      ${GSL_LIBRARIES}
#                      /cosmosis/cosmosis/datablock/libcosmosis.so)
#target_include_directories(default_cluster_abundance
#                           PRIVATE
#                           ${CMAKE_SOURCE_DIR}
#                           ${CMAKE_SOURCE_DIR}/externals
#                           ${CMAKE_SOURCE_DIR}/src
#                           ${CUBA_INCLUDE_DIR}
#                           ${CUBACPP_DIR}
#                           ${EIGEN3_INCLUDE_DIR})

add_library(cluster_abundance_covariance SHARED cluster_abundance_covariance.cc)
target_link_libraries(cluster_abundance_covariance
                      PRIVATE
                      utils
                      models
                      ${CUBA_LIBRARIES}
                      ${GSL_LIBRARIES}
                      /cosmosis/cosmosis/datablock/libcosmosis.so)
target_include_directories(cluster_abundance_covariance
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR})

add_library(ExampleVectorIntegration SHARED ExampleVectorIntegration_module.cc)
target_link_libraries(ExampleVectorIntegration
                     PRIVATE
                     /cosmosis/cosmosis/datablock/libcosmosis.so
                     example_vector_integrand
                     ${CUBA_LIBRARIES})                     

target_include_directories(ExampleVectorIntegration
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_library(sigma_miscent SHARED sigma_miscent.cc)
target_link_libraries(sigma_miscent PRIVATE /cosmosis/cosmosis/datablock/libcosmosis.so utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(sigma_miscent
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})
 
add_library(sigma_cent SHARED sigma_cent.cc)
target_link_libraries(sigma_cent PRIVATE /cosmosis/cosmosis/datablock/libcosmosis.so utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(sigma_cent
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_library(buzzard_sigma_halos SHARED buzzard_sigma_halos.cc)
target_link_libraries(buzzard_sigma_halos PRIVATE /cosmosis/cosmosis/datablock/libcosmosis.so utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(buzzard_sigma_halos
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})


add_library(sigma_cent_y1 SHARED sigma_cent_y1.cc)
target_link_libraries(sigma_cent_y1 PRIVATE /cosmosis/cosmosis/datablock/libcosmosis.so utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(sigma_cent_y1
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})


add_library(sigma_miscent_y1 SHARED sigma_miscent_y1.cc)
target_link_libraries(sigma_miscent_y1 PRIVATE /cosmosis/cosmosis/datablock/libcosmosis.so utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(sigma_miscent_y1
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_library(kappa_cent_y1 SHARED kappa_cent_y1.cc)
target_link_libraries(kappa_cent_y1 PRIVATE /cosmosis/cosmosis/datablock/libcosmosis.so utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(kappa_cent_y1
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_library(kappa_miscent_y1 SHARED kappa_miscent_y1.cc)
target_link_libraries(kappa_miscent_y1 PRIVATE /cosmosis/cosmosis/datablock/libcosmosis.so utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(kappa_miscent_y1
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})
