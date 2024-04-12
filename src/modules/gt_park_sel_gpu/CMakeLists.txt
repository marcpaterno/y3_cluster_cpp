include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

if (USE_CUDA)
  find_library(NVTX_LIBRARY nvToolsExt PATHS ENV LD_LIBRARY_PATH)

  # avgCentFull_sigma.cu
  add_library(avgCentSigmaPark   MODULE avgCentFull_sigma.cu)
  target_link_libraries(avgCentSigmaPark  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgCentSigmaPark  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgCentSigmaPark  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgCentSigmaPark  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # avgMiscApprox2_sigma.cu
  add_library(avgMiscApprox2SigmaPark   MODULE avgMiscApprox2_sigma.cu)
  target_link_libraries(avgMiscApprox2SigmaPark  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgMiscApprox2SigmaPark  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgMiscApprox2SigmaPark  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgMiscApprox2SigmaPark  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # # avgCentFull_cum_sigma.cu
  add_library(avgCentSigmaCumPark   MODULE avgCentFull_cum_sigma.cu)
  target_link_libraries(avgCentSigmaCumPark  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgCentSigmaCumPark  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgCentSigmaCumPark  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgCentSigmaCumPark  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # # avgMiscApprox2.cu
  add_library(avgMiscApprox2SigmaCumPark   MODULE avgMiscApprox2_cum_sigma.cu)
  target_link_libraries(avgMiscApprox2SigmaCumPark  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgMiscApprox2SigmaCumPark  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgMiscApprox2SigmaCumPark  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgMiscApprox2SigmaCumPark  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")


endif ()
