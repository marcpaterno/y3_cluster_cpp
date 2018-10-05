add_library(interp STATIC interp_1d.cc interp_2d.cc)
target_link_libraries(interp PRIVATE ${GSL_LIBRARIES})
set_property(TARGET interp PROPERTY POSITION_INDEPENDENT_CODE ON)

add_library(lc_lt STATIC lc_lt_t.cc)
target_link_libraries(lc_lt PRIVATE ${GSL_LIBRARIES})
target_include_directories(lc_lt PRIVATE ${CMAKE_SOURCE_DIR}/src)
set_property(TARGET lc_lt PROPERTY POSITION_INDEPENDENT_CODE ON)

add_library(smp_var
            STATIC
            sample_variance.cc
            bessel_polynomial_integral.cc
            sin_cos_polynomial_integral.cc)
target_link_libraries(smp_var PRIVATE ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(smp_var
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CMAKE_SOURCE_DIR}/../../cubacpp)
set_property(TARGET smp_var PROPERTY POSITION_INDEPENDENT_CODE ON)

add_library(example_module SHARED example_module.cc datablock_reader.cc interp_1d.cc interp_2d.cc)
target_link_libraries(example_module
		      PRIVATE
		      interp
		      lc_lt
		      ${CUBA_LIBRARIES}
		      ${GSL_LIBRARIES}
		      /cosmosis/cosmosis/datablock/libcosmosis.so)
target_include_directories(example_module
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBA_INCLUDE_DIR}
                           ${CMAKE_SOURCE_DIR}/../../cubacpp)
