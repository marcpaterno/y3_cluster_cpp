include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

add_library(compton_y_sims MODULE compton_y_sims.cc)
target_link_libraries(
  compton_y_sims PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                         ${CUBA_LIBRARIES}
  )
target_include_directories(
  compton_y_sims PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                         ${GSL_INCLUDE_DIRS}
  )
