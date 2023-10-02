include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

if (USE_CUDA)
  find_library(NVTX_LIBRARY nvToolsExt PATHS ENV LD_LIBRARY_PATH)

  add_library(numberCountsCardTriax   MODULE numberCounts.cu)
  target_link_libraries(numberCountsCardTriax  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils)
  target_include_directories(numberCountsCardTriax  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(numberCountsCardTriax  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(numberCountsCardTriax  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  
  # avgCentFull.cu
  #add_library(avgCentCardTriax   MODULE avgCentFull.cu)
  #target_link_libraries(avgCentCardTriax  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  #target_include_directories(avgCentCardTriax  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  #set_target_properties(avgCentCardTriax  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  #target_compile_options(avgCentCardTriax  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # avgMiscFull.cu
  #add_library(avgMiscFullCardTriax   MODULE avgMiscFull.cu)
  #target_link_libraries(avgMiscFullCardTriax  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  #target_include_directories(avgMiscFullCardTriax  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  #set_target_properties(avgMiscFullCardTriax  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  #target_compile_options(avgMiscFullCardTriax  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # # avgMiscApprox1.cu
  #add_library(avgMiscApprox1CardTriax   MODULE avgMiscApprox1.cu)
  #target_link_libraries(avgMiscApprox1CardTriax  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  # target_include_directories(avgMiscApprox1CardTriax  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  # set_target_properties(avgMiscApprox1CardTriax  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  # target_compile_options(avgMiscApprox1CardTriax  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # # # avgMiscApprox2.cu
  # add_library(avgMiscApprox2CardTriax   MODULE avgMiscApprox2.cu)
  # target_link_libraries(avgMiscApprox2CardTriax  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  # target_include_directories(avgMiscApprox2CardTriax  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  # set_target_properties(avgMiscApprox2CardTriax  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  # target_compile_options(avgMiscApprox2CardTriax  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")


endif ()
