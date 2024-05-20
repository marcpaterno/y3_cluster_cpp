include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

add_library(avgCentSigmaCumParkCPU MODULE avgCentFull_cum_sigma.cc)
target_link_libraries(avgCentSigmaCumParkCPU PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES} fmt::fmt)
target_include_directories(avgCentSigmaCumParkCPU PRIVATE ${PAGANI_DIR} ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})