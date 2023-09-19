include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

if (USE_CUDA)
  find_library(NVTX_LIBRARY nvToolsExt PATHS ENV LD_LIBRARY_PATH)

  add_library(numberCountsCard   MODULE numberCounts.cu)
  target_link_libraries(numberCountsCard  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils)
  target_include_directories(numberCountsCard  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(numberCountsCard  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(numberCountsCard  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # avgCentFull.cu
  add_library(avgCentCard   MODULE avgCentFull.cu)
  target_link_libraries(avgCentCard  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgCentCard  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgCentCard  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgCentCard  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # avgMiscFull.cu
  add_library(avgMiscFullCard   MODULE avgMiscFull.cu)
  target_link_libraries(avgMiscFullCard  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgMiscFullCard  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgMiscFullCard  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgMiscFullCard  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # # avgMiscApprox1.cu
  add_library(avgMiscApprox1Card   MODULE avgMiscApprox1.cu)
  target_link_libraries(avgMiscApprox1Card  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgMiscApprox1Card  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgMiscApprox1Card  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgMiscApprox1Card  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # # avgMiscApprox2.cu
  add_library(avgMiscApprox2Card   MODULE avgMiscApprox2.cu)
  target_link_libraries(avgMiscApprox2Card  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgMiscApprox2Card  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgMiscApprox2Card  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgMiscApprox2Card  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")


endif ()
