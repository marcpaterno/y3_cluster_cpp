add_executable(integrand_test main.cc integrand.test.cc)
target_link_libraries(integrand_test PRIVATE ${CUBA_LIBRARIES})
target_include_directories(integrand_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CMAKE_SOURCE_DIR}/../../cubacpp)
add_test(integrand_test integrand_test)

add_executable(integration_test main.cc integration.test.cc)
target_link_libraries(integration_test PRIVATE ${CUBA_LIBRARIES})
target_include_directories(integration_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CMAKE_SOURCE_DIR}/../../cubacpp)
add_test(integration_test integration_test)

add_executable(trivial_gamma_t trivial_gamma_t.cc interp_1d.cc)
target_link_libraries(trivial_gamma_t PRIVATE ${CUBA_LIBRARIES} ${GSL_LIBRARIES})
target_include_directories(trivial_gamma_t
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/examples
                           ${CUBA_INCLUDE_DIR}
                           ${CMAKE_SOURCE_DIR}/../../cubacpp)

add_executable(transform_test main.cc transform.test.cc)
target_include_directories(transform_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(transform_test transform_test)

add_executable(interp_1d_test main.cc interp_1d.test.cc interp_1d.cc)
target_link_libraries(interp_1d_test PRIVATE ${GSL_LIBRARIES})
target_include_directories(interp_1d_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/externals
                           ${GSL_INCLUDE_DIR})
add_test(interp_1d_test interp_1d_test)

add_executable(mz_power_law_test main.cc mz_power_law.test.cc)
target_include_directories(mz_power_law_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(mz_power_law_test mz_power_law_test)