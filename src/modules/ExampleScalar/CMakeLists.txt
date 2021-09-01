include_directories($ENV{COSMOSIS_SRC_DIR} ${CUBA_INCLUDE_DIR}) 
link_directories($ENV{COSMOSIS_SRC_DIR}/cosmosis/datablock)


add_library(example_scalar_integrand SHARED ExampleScalarIntegrand.cc)
target_link_libraries(example_scalar_integrand PRIVATE cosmosis utils)
target_include_directories(example_scalar_integrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})


# Note
# The following three "modules" are abnormal, in that they are declared as a
# shared libraries rather than as modules. This is so they can be used in testing the
# module templates themselves. This should not be done for actual modules, which should
# be declared as modules.

add_library(ExampleScalarIntegrationModule SHARED ExampleScalarIntegration_module.cc)
target_link_libraries(ExampleScalarIntegrationModule
                      PRIVATE
                      cosmosis
                      example_scalar_integrand
                      utils
                      ${CUBA_LIBRARIES}
		      ${GSL_LIBRARIES})
target_include_directories(ExampleScalarIntegrationModule
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})



add_executable(test_ExampleScalarIntegrationModule test_ExampleScalarIntegrationModule.cc)
target_link_libraries(test_ExampleScalarIntegrationModule
                      PRIVATE
                      cosmosis
                      example_scalar_integrand
                      ExampleScalarIntegrationModule
                      example_scalar_integrand
                      ${CUBA_LIBRARIES}
		      ${GSL_LIBRARIES})
target_include_directories(test_ExampleScalarIntegrationModule
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src)
add_test(test_ExampleScalarIntegrationModule test_ExampleScalarIntegrationModule)



add_executable(ExampleScalarIntegrand_test ExampleScalarIntegrand.test.cc)
target_link_libraries(ExampleScalarIntegrand_test
                      PRIVATE
                      example_scalar_integrand
                      cosmosis
                      ${CUBA_LIBRARIES}
                      testmain
		      ${GSL_LIBRARIES})
target_include_directories(ExampleScalarIntegrand_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR}
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(ExampleScalarIntegrand_test ExampleScalarIntegrand_test)





