include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

add_library(NumCountsFullScalarIntegrand MODULE NumCountsFullScalarIntegrand.cc)
target_link_libraries(NumCountsFullScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(NumCountsFullScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                                                ${GSL_INCLUDE_DIRS})
