add_library(testmain SHARED main.cc)
target_include_directories(testmain
	                   PUBLIC ${CMAKE_SOURCE_DIR}
				${CMAKE_SOURCE_DIR}/../src
				${CMAKE_SOURCE_DIR}/../externals
				${CMAKE_SOURCE_DIR}/src)

add_executable(xi_nfw_test xi_nfw.test.cc)
target_link_libraries(xi_nfw_test PRIVATE testmain)
target_include_directories(xi_nfw_test
                           PRIVATE
			   ${CMAKE_SOURCE_DIR}
		           ${CMAKE_SOURCE_DIR}/../src
		           ${CMAKE_SOURCE_DIR}/../externals
			   ${CMAKE_SOURCE_DIR}/src)
add_test(xi_nfw_test xi_nfw_test)


