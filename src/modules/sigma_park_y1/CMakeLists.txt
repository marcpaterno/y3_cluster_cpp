include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

if (USE_CUDA)
  find_library(NVTX_LIBRARY nvToolsExt PATHS ENV LD_LIBRARY_PATH)

  # summedNumbersCentY1.cu
  add_library(summedNumbersCentY1   MODULE summedNumbersCentY1.cu)
  target_link_libraries(summedNumbersCentY1  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models)
  target_include_directories(summedNumbersCentY1  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(summedNumbersCentY1  PROPERTIES POSITION_INDEPENDENT_CODE ON 
                                                                             CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(summedNumbersCentY1  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # summedNumbersMisCentY1.cu
  add_library(summedNumbersMiscentY1   MODULE summedNumbersMiscentY1.cu)
  target_link_libraries(summedNumbersMiscentY1  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models)
  target_include_directories(summedNumbersMiscentY1  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(summedNumbersMiscentY1  PROPERTIES POSITION_INDEPENDENT_CODE ON 
                                                                                CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(summedNumbersMiscentY1  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # avgSigmaCentY1.cu
  add_library(avgSigmaCentY1   MODULE avgSigmaCentY1.cu)
  target_link_libraries(avgSigmaCentY1  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models)
  set_target_properties(avgSigmaCentY1  PROPERTIES POSITION_INDEPENDENT_CODE ON 
                                                                        CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_include_directories(avgSigmaCentY1  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  target_compile_options(avgSigmaCentY1  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # avgSigmaMiscentY1.cu
  add_library(avgSigmaMiscentY1   MODULE avgSigmaMiscentY1.cu)
  target_link_libraries(avgSigmaMiscentY1  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models)
  set_target_properties(avgSigmaMiscentY1  PROPERTIES POSITION_INDEPENDENT_CODE ON
                                                                           CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_include_directories(avgSigmaMiscentY1  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  target_compile_options(avgSigmaMiscentY1  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # avgWpCent.cu
  add_library(avgWpCentY1   MODULE avgWpCentY1.cu)
  target_link_libraries(avgWpCentY1  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models)
  set_target_properties(avgWpCentY1  PROPERTIES POSITION_INDEPENDENT_CODE ON 
                                                                   CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_include_directories(avgWpCentY1  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  target_compile_options(avgWpCentY1  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # avgWpMiscentY1.cu
  add_library(avgWpMiscentY1   MODULE avgWpMiscentY1.cu)
  target_link_libraries(avgWpMiscentY1  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models)
  set_target_properties(avgWpMiscentY1  PROPERTIES POSITION_INDEPENDENT_CODE ON
                                                                      CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_include_directories(avgWpMiscentY1  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  target_compile_options(avgWpMiscentY1  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")


endif ()
