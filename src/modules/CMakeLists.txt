add_library(example_integrand SHARED ExampleIntegrand.cc)
target_link_libraries(example_integrand PRIVATE /cosmosis/cosmosis/datablock/libcosmosis.so)
target_include_directories(example_integrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_executable(ExampleIntegrand_test ExampleIntegrand.test.cc)
target_link_libraries(ExampleIntegrand_test
                      PRIVATE
                      example_integrand
                      /cosmosis/cosmosis/datablock/libcosmosis.so
                      ${CUBA_LIBRARIES}
                      testmain)
target_include_directories(ExampleIntegrand_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR}
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(ExampleIntegrand_test ExampleIntegrand_test)

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

add_library(ExampleIntegration SHARED ExampleIntegration_module.cc)
target_link_libraries(ExampleIntegration
                     PRIVATE
                     /cosmosis/cosmosis/datablock/libcosmosis.so
                     example_integrand
                     ${CUBA_LIBRARIES})                     

target_include_directories(ExampleIntegration
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

