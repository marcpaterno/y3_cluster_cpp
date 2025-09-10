include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

add_library(y1_analysis MODULE y1_analysis.cc)
target_link_libraries(
  y1_analysis PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES}
  )
target_include_directories(
  y1_analysis PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                      ${GSL_INCLUDE_DIRS}
  )

add_library(y1_analysis_mor MODULE y1_analysis_mor.cc)
target_link_libraries(
  y1_analysis_mor PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                          ${CUBA_LIBRARIES}
  )
target_include_directories(
  y1_analysis_mor PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                          ${GSL_INCLUDE_DIRS}
  )

add_library(y1_analysis_mor_2022 MODULE y1_analysis_mor_2022.cc)
target_link_libraries(
  y1_analysis_mor_2022 PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                          ${CUBA_LIBRARIES}
  )
target_include_directories(
  y1_analysis_mor_2022 PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                          ${GSL_INCLUDE_DIRS}
  )

add_library(y1_analysis_mor_logm_2022 MODULE y1_analysis_mor_logm_2022.cc)
target_link_libraries(
  y1_analysis_mor_logm_2022 PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                          ${CUBA_LIBRARIES}
  )
target_include_directories(
  y1_analysis_mor_logm_2022 PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                          ${GSL_INCLUDE_DIRS}
  )

