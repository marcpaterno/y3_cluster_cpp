include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

add_library(sigma_nl MODULE sigma_nl.cc)
target_link_libraries(
  sigma_nl PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES}
  )
target_include_directories(
  sigma_nl PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS}
  )

add_library(model_sigmahm MODULE model_sigmahm.cc)
target_link_libraries(
  model_sigmahm PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                        ${CUBA_LIBRARIES}
  )
target_include_directories(
  model_sigmahm PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                        ${GSL_INCLUDE_DIRS}
  )
