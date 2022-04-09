include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

add_library(NCMiscentY1ScalarIntegrand MODULE NC_miscent_y1_scalarintegrand.cc)
target_link_libraries(
  NCMiscentY1ScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                                     ${CUBA_LIBRARIES} ${GSL_LIBS}
  )
target_include_directories(
  NCMiscentY1ScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                     ${GSL_INCLUDE_DIRS}
  )

add_library(NCCentY1ScalarIntegrand MODULE NC_cent_y1_scalarintegrand.cc)
target_link_libraries(
  NCCentY1ScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                                  ${CUBA_LIBRARIES} ${GSL_LIBS}
  )
target_include_directories(
  NCCentY1ScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                  ${GSL_INCLUDE_DIRS}
  )

add_library(
  NCRadMiscentY1ScalarIntegrand MODULE NCRad_miscent_y1_scalarintegrand.cc
  )
target_link_libraries(
  NCRadMiscentY1ScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                                        ${CUBA_LIBRARIES} ${GSL_LIBS}
  )
target_include_directories(
  NCRadMiscentY1ScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                        ${GSL_INCLUDE_DIRS}
  )

add_library(
  MassMiscentY1ScalarIntegrand MODULE Mass_miscent_y1_scalarintegrand.cc
  )
target_link_libraries(
  MassMiscentY1ScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                                       ${CUBA_LIBRARIES} ${GSL_LIBS}
  )
target_include_directories(
  MassMiscentY1ScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                       ${GSL_INCLUDE_DIRS}
  )

add_library(MassCentY1ScalarIntegrand MODULE Mass_cent_y1_scalarintegrand.cc)
target_link_libraries(
  MassCentY1ScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                                    ${CUBA_LIBRARIES} ${GSL_LIBS}
  )
target_include_directories(
  MassCentY1ScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                    ${GSL_INCLUDE_DIRS}
  )
