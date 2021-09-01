include_directories($ENV{COSMOSIS_SRC_DIR} ${CUBA_INCLUDE_DIR}) 
link_directories($ENV{COSMOSIS_SRC_DIR}/cosmosis/datablock)


add_library(SigmaCentY1MortScalarIntegrand MODULE sigma_cent_y1_mort_scalarintegrand.cc)
target_link_libraries(SigmaCentY1MortScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(SigmaCentY1MortScalarIntegrand
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})

add_library(SigmaMiscentY1MortScalarIntegrand MODULE  sigma_miscent_y1_mort_scalarintegrand.cc)
target_link_libraries(SigmaMiscentY1MortScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(SigmaMiscentY1MortScalarIntegrand
                           PRIVATE
			   ${PAGANI_DIR}
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${EIGEN3_INCLUDE_DIR})



