add_library(interp STATIC interp_1d.cc interp_2d.cc)
target_link_libraries(interp PRIVATE ${GSL_LIBRARIES})
set_property(TARGET interp PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(interp
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/../externals
                           ${CMAKE_SOURCE_DIR}/../src
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/../../../cubacpp)

add_library(xihm SHARED xihm.cc interp_1d.cc interp_2d.cc)
target_link_libraries(xihm
                      PRIVATE
		      interp
                      /cosmosis/cosmosis/datablock/libcosmosis.so
		      ${CUBA_LIBRARIES})
target_include_directories(xihm
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/../externals
                           ${CMAKE_SOURCE_DIR}/../src
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/../../../cubacpp)
