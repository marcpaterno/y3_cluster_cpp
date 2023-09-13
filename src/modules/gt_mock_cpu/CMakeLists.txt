include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

add_library(avgCentFullScalarIntegrand MODULE avgCentFull.cc)
target_link_libraries(avgCentFullScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(avgCentFullScalarIntegrand PRIVATE ${PAGANI_DIR} ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})

add_library(avgMiscFullScalarIntegrand MODULE avgMiscFull.cc)
target_link_libraries(avgMiscFullScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(avgMiscFullScalarIntegrand PRIVATE ${PAGANI_DIR} ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})

add_library(avgMiscApprox1ScalarIntegrand MODULE avgMiscApprox1.cc)
target_link_libraries(avgMiscApprox1ScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(avgMiscApprox1ScalarIntegrand PRIVATE ${PAGANI_DIR} ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})

add_library(avgMiscApprox2ScalarIntegrand MODULE avgMiscApprox1.cc)
target_link_libraries(avgMiscApprox2ScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(avgMiscApprox2ScalarIntegrand PRIVATE ${PAGANI_DIR} ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})