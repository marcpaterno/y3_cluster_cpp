add_library(interp SHARED interp_1d.cc interp_2d.cc)
target_link_libraries(interp PRIVATE ${GSL_LIBRARIES})
    
add_library(lc_lt SHARED lc_lt_t.cc)
target_link_libraries(lc_lt PRIVATE ${GSL_LIBRARIES})
target_include_directories(lc_lt PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_library(example_integrand SHARED ExampleIntegrand.cc)
target_link_libraries(example_integrand PRIVATE /cosmosis/cosmosis/datablock/libcosmosis.so)
target_include_directories(example_integrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/../../cubacpp
                           ${EIGEN3_INCLUDE_DIR})

add_library(smp_var
            SHARED
            sample_variance.cc
            bessel_polynomial_integral.cc
            sin_cos_polynomial_integral.cc)
target_link_libraries(smp_var
                      PRIVATE
                      interp
                      ${GSL_LIBRARIES}
                  ${CUBA_LIBRARIES})
target_include_directories(smp_var
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CMAKE_SOURCE_DIR}/../../cubacpp
                           ${EIGEN3_INCLUDE_DIR})

add_library(default_cluster_abundance SHARED default_cluster_abundance.cc)
target_link_libraries(default_cluster_abundance
                      PRIVATE
                      interp
                      lc_lt
                      profiler
                      cfitsio
                      ${CUBA_LIBRARIES}
                      ${GSL_LIBRARIES}
                      /cosmosis/cosmosis/datablock/libcosmosis.so)
target_include_directories(default_cluster_abundance
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBA_INCLUDE_DIR}
                           ${CMAKE_SOURCE_DIR}/../../cubacpp
                           ${EIGEN3_INCLUDE_DIR})

add_library(cluster_abundance_covariance SHARED cluster_abundance_covariance.cc)
target_link_libraries(cluster_abundance_covariance
                      PRIVATE
                      interp
                      lc_lt
                      smp_var
                      ${CUBA_LIBRARIES}
                      ${GSL_LIBRARIES}
                      /cosmosis/cosmosis/datablock/libcosmosis.so)
target_include_directories(cluster_abundance_covariance
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBA_INCLUDE_DIR}
                           ${CMAKE_SOURCE_DIR}/../../cubacpp)

add_library(ExampleIntegration SHARED ExampleIntegration_module.cc)
target_link_libraries(ExampleIntegration
                     PRIVATE
                     example_integrand)
target_include_directories(ExampleIntegration
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CUBA_INCLUDE_DIR}
                           ${CMAKE_SOURCE_DIR}/../../cubacpp
                           ${EIGEN3_INCLUDE_DIR})