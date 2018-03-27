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

add_executable(trivial_gamma_t trivial_gamma_t.cc)
target_link_libraries(trivial_gamma_t PRIVATE ${CUBA_LIBRARIES})
target_include_directories(trivial_gamma_t
	                   PRIVATE
			   ${CMAKE_SOURCE_DIR}
			   ${CMAKE_SOURCE_DIR}/examples
			   ${CUBA_INCLUDE_DIR}
			   ${CMAKE_SOURCE_DIR}/../../cubacpp)

