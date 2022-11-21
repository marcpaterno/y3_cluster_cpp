include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

if (USE_CUDA)
  find_library(NVTX_LIBRARY nvToolsExt PATHS ENV LD_LIBRARY_PATH)

  # summedNumbersCentBu.cu
  add_library(summedNumbersCentBu   MODULE summedNumbersCentBu.cu)
  target_link_libraries(summedNumbersCentBu  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models)
  target_include_directories(summedNumbersCentBu  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(summedNumbersCentBu  PROPERTIES POSITION_INDEPENDENT_CODE ON 
                                                                             CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(summedNumbersCentBu  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # avgKappaCentBu.cu
  add_library(avgKappaCentBu   MODULE avgKappaCentBu.cu)
  target_link_libraries(avgKappaCentBu  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models)
  target_include_directories(avgKappaCentBu  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgKappaCentBu  PROPERTIES POSITION_INDEPENDENT_CODE ON 
                                                                             CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgKappaCentBu  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # avgSigmaCentBu.cu
  add_library(avgSigmaCentBu   MODULE avgSigmaCentBu.cu)
  target_link_libraries(avgSigmaCentBu  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models)
  set_target_properties(avgSigmaCentBu  PROPERTIES POSITION_INDEPENDENT_CODE ON 
                                                                        CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_include_directories(avgSigmaCentBu  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  target_compile_options(avgSigmaCentBu  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # avgPkDampCentBU
  add_library(avgPkDampCentBU   MODULE avgPkDampCentBU.cu)
  target_link_libraries(avgPkDampCentBU  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models)
  set_target_properties(avgPkDampCentBU  PROPERTIES POSITION_INDEPENDENT_CODE ON 
                                                                        CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_include_directories(avgPkDampCentBU  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  target_compile_options(avgPkDampCentBU  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # avgWpCentBu.cu
  add_library(avgWpCentBu   MODULE avgWpCentBu.cu)
  target_link_libraries(avgWpCentBu  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models)
  set_target_properties(avgWpCentBu  PROPERTIES POSITION_INDEPENDENT_CODE ON 
                                                                   CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_include_directories(avgWpCentBu  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  target_compile_options(avgWpCentBu  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")



endif ()
