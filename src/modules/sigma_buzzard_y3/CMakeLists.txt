include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

if (USE_CUDA)
  find_library(NVTX_LIBRARY nvToolsExt PATHS ENV LD_LIBRARY_PATH)


  # dummy.cu
  add_library(dummy   MODULE dummy.cu)
  target_link_libraries(dummy  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models)
  target_include_directories(dummy  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(dummy  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(dummy  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")


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

  # avgGammaFullBu.cu
  add_library(avgGammaFullBu   MODULE avgGammaFullBu.cu)
  target_link_libraries(avgGammaFullBu  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgGammaFullBu  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgGammaFullBu  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgGammaFullBu  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # # avgGammaCentBu.cu
  add_library(avgGammaCentBu   MODULE avgGammaCentBu.cu)
  target_link_libraries(avgGammaCentBu  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgGammaCentBu  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgGammaCentBu  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgGammaCentBu  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # # avgGammaMiscBu.cu
  add_library(avgGammaMiscBu   MODULE avgGammaMiscBu.cu)
  target_link_libraries(avgGammaMiscBu  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgGammaMiscBu  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgGammaMiscBu  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgGammaMiscBu  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # avgGammaProjBu.cu
  add_library(avgGammaProjBu   MODULE avgGammaProjBu.cu)
  target_link_libraries(avgGammaProjBu  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models cuda_utils fmt::fmt)
  target_include_directories(avgGammaProjBu  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(avgGammaProjBu  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(avgGammaProjBu  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  # avgSigmaCentBu.cu
  add_library(avgSigmaCentBu   MODULE avgSigmaCentBu.cu)
  target_link_libraries(avgSigmaCentBu  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models)
  set_target_properties(avgSigmaCentBu  PROPERTIES POSITION_INDEPENDENT_CODE ON 
                                                                        CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_include_directories(avgSigmaCentBu  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  target_compile_options(avgSigmaCentBu  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")
  
  # avgSigmaMiscBu.cu
  add_library(avgSigmaMiscBu   MODULE avgSigmaMiscBu.cu)
  target_link_libraries(avgSigmaMiscBu  PRIVATE cosmosis ${NVTX_LIBRARY} utils cuda_models)
  set_target_properties(avgSigmaMiscBu  PROPERTIES POSITION_INDEPENDENT_CODE ON 
                                                                        CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_include_directories(avgSigmaMiscBu  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  target_compile_options(avgSigmaMiscBu  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")
  
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
