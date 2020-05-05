include_directories($ENV{COSMOSIS_SRC_DIR} ${CUBA_INCLUDE_DIR}) 
link_directories($ENV{COSMOSIS_SRC_DIR}/cosmosis/datablock)

add_library(example_vector_integrand SHARED ExampleVectorIntegrand.cc)
target_link_libraries(example_vector_integrand PRIVATE cosmosis utils)
target_include_directories(example_vector_integrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR}
                           ${GSL_INCLUDE_DIRS})

add_library(example_scalar_integrand SHARED ExampleScalarIntegrand.cc)
target_link_libraries(example_scalar_integrand PRIVATE cosmosis utils)
target_include_directories(example_scalar_integrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_library(example_oned_integrand SHARED ExampleOneDIntegrand.cc)
target_link_libraries(example_oned_integrand PRIVATE cosmosis utils)
target_include_directories(example_oned_integrand PRIVATE ${CMAKE_SOURCE_DIR}/src) 

# Note
# The following three "modules" are abnormal, in that they are declared as a
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
                           ${GSL_INCLUDE_DIRS}
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
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})

add_library(ExampleOneDIntegrationModule SHARED ExampleOneDIntegration_module.cc)
target_link_libraries(ExampleOneDIntegrationModule
                      PRIVATE
                      cosmosis
                      example_oned_integrand
                      ${CUBA_LIBRARIES}
                      ${GSL_LIBRARIES})
target_include_directories(ExampleOneDIntegrationModule
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR}
                           ${GSL_INCLUDE_DIRS})

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

add_executable(test_ExampleOneDIntegrationModule test_ExampleOneDIntegrationModule.cc)
target_link_libraries(test_ExampleOneDIntegrationModule
                      PRIVATE
                      cosmosis
                      example_oned_integrand
                      ExampleOneDIntegrationModule
                      example_oned_integrand)
target_include_directories(test_ExampleOneDIntegrationModule
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src)
add_test(test_ExampleOneDIntegrationModule test_ExampleOneDIntegrationModule)

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

add_executable(ExampleOneDIntegrand_test ExampleOneDIntegrand.test.cc)
target_link_libraries(ExampleOneDIntegrand_test
                      PRIVATE
                      example_oned_integrand
                      cosmosis
                      ${CUBA_LIBRARIES}
                      ExampleOneDIntegrationModule
                      testmain)
target_include_directories(ExampleOneDIntegrand_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR}
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(ExampleOneDIntegrand_test ExampleOneDIntegrand_test)

#### The above are example modules

add_subdirectory(mass_y1)
add_subdirectory(mass_mort_y1)
add_subdirectory(sigma_mort_y1)
add_subdirectory(y1_analysis)
add_subdirectory(sdss_analysis)
add_subdirectory(snapshotsim)
add_subdirectory(cluster_abundance_covariance)
add_subdirectory(DESxSPTModule)
add_subdirectory(sigma_kappa_y1)
add_subdirectory(buzzard_test)
add_subdirectory(compton_y_sims)
add_subdirectory(model_sigmahm)

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
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})



