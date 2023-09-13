include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

if (USE_CUDA)
  find_library(NVTX_LIBRARY nvToolsExt PATHS ENV LD_LIBRARY_PATH)

  # avgCentFull.cu
  add_library(avgCentFull   MODULE avgCentFull.cu)
  target_link_libraries(avgCentFull  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgCentFull  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgCentFull  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgCentFull  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # avgMiscFull.cu
  add_library(avgMiscFull   MODULE avgMiscFull.cu)
  target_link_libraries(avgMiscFull  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgMiscFull  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgMiscFull  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgMiscFull  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # # avgMiscApprox1.cu
  add_library(avgMiscApprox1   MODULE avgMiscApprox1.cu)
  target_link_libraries(avgMiscApprox1  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgMiscApprox1  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgMiscApprox1  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgMiscApprox1  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # # avgMiscApprox2.cu
  add_library(avgMiscApprox2   MODULE avgMiscApprox2.cu)
  target_link_libraries(avgMiscApprox2  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgMiscApprox2  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgMiscApprox2  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgMiscApprox2  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")


endif ()
