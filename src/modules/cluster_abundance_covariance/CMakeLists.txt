include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

add_library(cluster_abundance_covariance MODULE cluster_abundance_covariance.cc)
target_link_libraries(
  cluster_abundance_covariance PRIVATE utils models ${CUBA_LIBRARIES}
                                       ${GSL_LIBRARIES} cosmosis
  )
target_include_directories(
  cluster_abundance_covariance
  PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/externals
          ${CMAKE_SOURCE_DIR}/src ${GSL_INCLUDE_DIRS} ${CUBA_INCLUDE_DIR}
          ${CUBACPP_DIR}
  )
