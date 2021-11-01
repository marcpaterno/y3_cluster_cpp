include_directories($ENV{COSMOSIS_SRC_DIR} ${CUBA_INCLUDE_DIR})
link_directories($ENV{COSMOSIS_SRC_DIR}/cosmosis/datablock)

add_library(sdss_analysis MODULE sdss_analysis.cc)
target_link_libraries(
  sdss_analysis PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                        ${CUBA_LIBRARIES}
  )
target_include_directories(
  sdss_analysis PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                        ${GSL_INCLUDE_DIRS}
  )
