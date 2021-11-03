include_directories($ENV{COSMOSIS_SRC_DIR} ${CUBA_INCLUDE_DIR})
link_directories($ENV{COSMOSIS_SRC_DIR}/cosmosis/datablock)

add_library(NCCentY1MortScalarIntegrand MODULE NC_cent_y1_mort_scalarintegrand.cc)
target_link_libraries(NCCentY1MortScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(NCCentY1MortScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                                               ${GSL_INCLUDE_DIRS})

add_library(NCMiscentY1MortScalarIntegrand MODULE NC_miscent_y1_mort_scalarintegrand.cc)
target_link_libraries(NCMiscentY1MortScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(NCMiscentY1MortScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                                                  ${GSL_INCLUDE_DIRS})

if (USE_CUDA)
  find_library(NVTX_LIBRARY nvToolsExt PATHS ENV LD_LIBRARY_PATH)

  add_library(NCCentY1MortCUDAScalarIntegrand MODULE NC_cent_y1_mort_cuda_scalarintegrand.cu)
  target_link_libraries(NCCentY1MortCUDAScalarIntegrand PRIVATE cosmosis ${NVTX_LIBRARY} utils)
  target_include_directories(NCCentY1MortCUDAScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                                                     ${PAGANI_DIR})
  set_target_properties(NCCentY1MortCUDAScalarIntegrand PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES
                                                                                                ${Y3GCC_TARGET_ARCH})
  target_compile_options(NCCentY1MortCUDAScalarIntegrand PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  add_library(NCMiscentY1MortCUDAScalarIntegrand MODULE NC_miscent_y1_mort_cuda_scalarintegrand.cu)
  target_link_libraries(NCMiscentY1MortCUDAScalarIntegrand PRIVATE cosmosis ${NVTX_LIBRARY} utils)
  target_include_directories(NCMiscentY1MortCUDAScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                                                        ${PAGANI_DIR})
  set_target_properties(NCMiscentY1MortCUDAScalarIntegrand PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES
                                                                                                   ${Y3GCC_TARGET_ARCH})
  target_compile_options(NCMiscentY1MortCUDAScalarIntegrand PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")
endif ()

add_library(NCRadMiscentY1MortScalarIntegrand MODULE NCRad_miscent_y1_mort_scalarintegrand.cc)
target_link_libraries(NCRadMiscentY1MortScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                                                                ${CUBA_LIBRARIES})
target_include_directories(NCRadMiscentY1MortScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                                                     ${GSL_INCLUDE_DIRS})

add_library(MassRadMiscentY1MortScalarIntegrand MODULE MassRad_miscent_y1_mort_scalarintegrand.cc)
target_link_libraries(MassRadMiscentY1MortScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                                                                  ${CUBA_LIBRARIES})
target_include_directories(MassRadMiscentY1MortScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                                                       ${GSL_INCLUDE_DIRS})

add_library(MassMiscentY1MortScalarIntegrand MODULE Mass_miscent_y1_mort_scalarintegrand.cc)
target_link_libraries(MassMiscentY1MortScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(MassMiscentY1MortScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                                                    ${GSL_INCLUDE_DIRS})

add_library(MassCentY1MortScalarIntegrand MODULE Mass_cent_y1_mort_scalarintegrand.cc)
target_link_libraries(MassCentY1MortScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(MassCentY1MortScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                                                 ${GSL_INCLUDE_DIRS})
