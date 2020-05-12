include_directories($ENV{COSMOSIS_SRC_DIR} ${CUBA_INCLUDE_DIR}) 
link_directories($ENV{COSMOSIS_SRC_DIR}/cosmosis/datablock)

# Note
# The following three "modules" are abnormal, in that they are declared as a
# shared libraries rather than as modules. This is so they can be used in testing the
# module templates themselves. This should not be done for actual modules, which should
# be declared as modules.

add_library(example_oned_integrand SHARED ExampleOneDIntegrand.cc)
target_link_libraries(example_oned_integrand PRIVATE cosmosis utils)
target_include_directories(example_oned_integrand PRIVATE ${CMAKE_SOURCE_DIR}/src)

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





