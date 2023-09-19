include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

if (USE_CUDA)
  find_library(NVTX_LIBRARY nvToolsExt PATHS ENV LD_LIBRARY_PATH)

  add_library(numberCountsMock   MODULE numberCounts.cu)
  target_link_libraries(numberCountsMock  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils)
  target_include_directories(numberCountsMock  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(numberCountsMock  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(numberCountsMock  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # avgCentFull.cu
  add_library(avgCentMock   MODULE avgCentFull.cu)
  target_link_libraries(avgCentMock  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgCentMock  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgCentMock  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgCentMock  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # avgMiscFull.cu
  add_library(avgMiscFullMock   MODULE avgMiscFull.cu)
  target_link_libraries(avgMiscFullMock  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgMiscFullMock  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgMiscFullMock  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgMiscFullMock  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # # avgMiscApprox1.cu
  add_library(avgMiscApprox1Mock   MODULE avgMiscApprox1.cu)
  target_link_libraries(avgMiscApprox1Mock  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgMiscApprox1Mock  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgMiscApprox1Mock  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgMiscApprox1Mock  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # # avgMiscApprox2.cu
  add_library(avgMiscApprox2Mock   MODULE avgMiscApprox2.cu)
  target_link_libraries(avgMiscApprox2Mock  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgMiscApprox2Mock  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgMiscApprox2Mock  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgMiscApprox2Mock  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")


endif ()
