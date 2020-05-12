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



