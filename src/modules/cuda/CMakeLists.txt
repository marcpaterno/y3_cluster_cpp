find_library(NVTX_LIBRARY nvToolsExt PATHS ENV LD_LIBRARY_PATH)

add_library(SigmaMiscentY1MortCUDAScalarIntegrand MODULE sigma_miscent_y1_mor_scalarintegrand.cu)
target_link_libraries(SigmaMiscentY1MortCUDAScalarIntegrand PRIVATE cosmosis ${NVTX_LIBRARY} utils)
set_target_properties(SigmaMiscentY1MortCUDAScalarIntegrand PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
target_include_directories(SigmaMiscentY1MortCUDAScalarIntegrand PRIVATE
  ${CMAKE_SOURCE_DIR}/src
 ${CUBACPP_DIR}
 ${PAGANI_DIR})
target_compile_options(SigmaMiscentY1MortCUDAScalarIntegrand PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

add_executable(constexpr_test constexpr_test.cu)
set_target_properties(constexpr_test PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
#target_compile_options(constexpr_test PRIVATE "--expt-relaxed-constexpr")
