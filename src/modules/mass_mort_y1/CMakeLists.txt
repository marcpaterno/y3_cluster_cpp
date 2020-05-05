include_directories($ENV{COSMOSIS_SRC_DIR} ${CUBA_INCLUDE_DIR}) 
link_directories($ENV{COSMOSIS_SRC_DIR}/cosmosis/datablock)


add_library(NCMiscentY1MortScalarIntegrand MODULE NC_miscent_y1_mort_scalarintegrand.cc)
target_link_libraries(NCMiscentY1MortScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(NCMiscentY1MortScalarIntegrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})

add_library(NCCentY1MortScalarIntegrand MODULE NC_cent_y1_mort_scalarintegrand.cc)
target_link_libraries(NCCentY1MortScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(NCCentY1MortScalarIntegrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})

add_library(NCRadMiscentY1MortScalarIntegrand MODULE NCRad_miscent_y1_mort_scalarintegrand.cc)
target_link_libraries(NCRadMiscentY1MortScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(NCRadMiscentY1MortScalarIntegrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})

add_library(MassMiscentY1MortScalarIntegrand MODULE Mass_miscent_y1_mort_scalarintegrand.cc)
target_link_libraries(MassMiscentY1MortScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(MassMiscentY1MortScalarIntegrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})

add_library(MassCentY1MortScalarIntegrand MODULE Mass_cent_y1_mort_scalarintegrand.cc)
target_link_libraries(MassCentY1MortScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(MassCentY1MortScalarIntegrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})


