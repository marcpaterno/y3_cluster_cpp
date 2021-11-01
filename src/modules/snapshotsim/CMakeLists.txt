include_directories($ENV{COSMOSIS_SRC_DIR} ${CUBA_INCLUDE_DIR}) 
link_directories($ENV{COSMOSIS_SRC_DIR}/cosmosis/datablock)

add_library(snapshotsim_scalarintegrand_NC MODULE Snapshotsim_ScalarIntegrand_NC.cc)
target_link_libraries(snapshotsim_scalarintegrand_NC PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(snapshotsim_scalarintegrand_NC
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS})

add_library(snapshotsim_scalarintegrand MODULE Snapshotsim_ScalarIntegrand.cc)
target_link_libraries(snapshotsim_scalarintegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(snapshotsim_scalarintegrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS})


add_library(snapshotsim_scalarintegrand_sigma MODULE Snapshotsim_ScalarIntegrand_Sigma.cc)
target_link_libraries(snapshotsim_scalarintegrand_sigma PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(snapshotsim_scalarintegrand_sigma
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS})

add_library(snapshotsim_scalarintegrand_max MODULE Snapshotsim_ScalarIntegrand_max.cc)
target_link_libraries(snapshotsim_scalarintegrand_max PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(snapshotsim_scalarintegrand_max
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS})

add_library(snapshotsim_scalarintegrand_max_sigma MODULE Snapshotsim_ScalarIntegrand_maxxi_sigma.cc)
target_link_libraries(snapshotsim_scalarintegrand_max_sigma PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(snapshotsim_scalarintegrand_max_sigma
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS})


