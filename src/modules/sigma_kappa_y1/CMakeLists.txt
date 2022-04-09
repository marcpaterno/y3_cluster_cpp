include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

add_library(
  SigmaMiscentY1ScalarIntegrand MODULE sigma_miscent_y1_scalarintegrand.cc
  )
target_link_libraries(
  SigmaMiscentY1ScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                                        ${CUBA_LIBRARIES}
  )
target_include_directories(
  SigmaMiscentY1ScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                        ${GSL_INCLUDE_DIRS}
  )

add_library(SigmaCentY1ScalarIntegrand MODULE sigma_cent_y1_scalarintegrand.cc)
target_link_libraries(
  SigmaCentY1ScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                                     ${CUBA_LIBRARIES}
  )
target_include_directories(
  SigmaCentY1ScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                     ${GSL_INCLUDE_DIRS}
  )

add_library(sigma_miscent MODULE sigma_miscent.cc)
target_link_libraries(
  sigma_miscent PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                        ${CUBA_LIBRARIES}
  )
target_include_directories(
  sigma_miscent PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                        ${GSL_INCLUDE_DIRS}
  )

add_library(sigma_cent MODULE sigma_cent.cc)
target_link_libraries(
  sigma_cent PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES}
  )
target_include_directories(
  sigma_cent PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS}
  )

add_library(sigma_cent_y1 MODULE sigma_cent_y1.cc)
target_link_libraries(
  sigma_cent_y1 PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                        ${CUBA_LIBRARIES}
  )
target_include_directories(
  sigma_cent_y1 PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                        ${GSL_INCLUDE_DIRS}
  )

add_library(sigma_miscent_y1 MODULE sigma_miscent_y1.cc)
target_link_libraries(
  sigma_miscent_y1 PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                           ${CUBA_LIBRARIES}
  )
target_include_directories(
  sigma_miscent_y1 PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
  )

add_library(kappa_cent_y1 MODULE kappa_cent_y1.cc)
target_link_libraries(
  kappa_cent_y1 PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                        ${CUBA_LIBRARIES}
  )
target_include_directories(
  kappa_cent_y1 PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                        ${GSL_INCLUDE_DIRS}
  )

add_library(kappa_miscent_y1 MODULE kappa_miscent_y1.cc)
target_link_libraries(
  kappa_miscent_y1 PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                           ${CUBA_LIBRARIES}
  )
target_include_directories(
  kappa_miscent_y1 PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
  )
