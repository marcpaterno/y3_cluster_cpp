include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

add_library(numberCountsCardCPU MODULE numberCounts.cc)
target_link_libraries(numberCountsCardCPU PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES} fmt::fmt)
target_include_directories(numberCountsCardCPU PRIVATE ${PAGANI_DIR} ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})

add_library(avgCentCardCPU MODULE avgCentFull.cc)
target_link_libraries(avgCentCardCPU PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(avgCentCardCPU PRIVATE ${PAGANI_DIR} ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})

add_library(avgMiscFullCardCPU MODULE avgMiscFull.cc)
target_link_libraries(avgMiscFullCardCPU PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(avgMiscFullCardCPU PRIVATE ${PAGANI_DIR} ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})

add_library(avgMiscApprox1CardCPU MODULE avgMiscApprox1.cc)
target_link_libraries(avgMiscApprox1CardCPU PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(avgMiscApprox1CardCPU PRIVATE ${PAGANI_DIR} ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})

add_library(avgMiscApprox2CardCPU MODULE avgMiscApprox2.cc)
target_link_libraries(avgMiscApprox2CardCPU PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES} fmt::fmt)
target_include_directories(avgMiscApprox2CardCPU PRIVATE ${PAGANI_DIR} ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})