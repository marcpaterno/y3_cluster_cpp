include_directories($ENV{COSMOSIS_SRC_DIR} ${CUBA_INCLUDE_DIR})
link_directories($ENV{COSMOSIS_SRC_DIR}/cosmosis/datablock)

add_library(
  SigmaCentY1MortScalarIntegrand MODULE sigma_cent_y1_mort_scalarintegrand.cc
  )
target_link_libraries(
  SigmaCentY1MortScalarIntegrand PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                                         ${CUBA_LIBRARIES}
  )
target_include_directories(
  SigmaCentY1MortScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                         ${GSL_INCLUDE_DIRS}
  )

add_library(
  SigmaMiscentY1MortScalarIntegrand MODULE
  sigma_miscent_y1_mort_scalarintegrand.cc
  )
target_link_libraries(
  SigmaMiscentY1MortScalarIntegrand PRIVATE cosmosis utils models
                                            ${GSL_LIBRARIES} ${CUBA_LIBRARIES}
  )
target_include_directories(
  SigmaMiscentY1MortScalarIntegrand
  PRIVATE ${PAGANI_DIR} ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
          ${GSL_INCLUDE_DIRS}
  )

if (USE_CUDA)
  find_library(NVTX_LIBRARY nvToolsExt PATHS ENV LD_LIBRARY_PATH)
  add_library(
    SigmaMiscentY1MortCUDAScalarIntegrand MODULE
    sigma_miscent_y1_mort_cuda_scalarintegrand.cu
    )
  target_link_libraries(
    SigmaMiscentY1MortCUDAScalarIntegrand PRIVATE cosmosis ${NVTX_LIBRARY}
                                                  utils
    )
  set_target_properties(
    SigmaMiscentY1MortCUDAScalarIntegrand
    PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES
                                            ${Y3GCC_TARGET_ARCH}
    )
  target_include_directories(
    SigmaMiscentY1MortCUDAScalarIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src
                                                  ${CUBACPP_DIR} ${PAGANI_DIR}
    )
  target_compile_options(
    SigmaMiscentY1MortCUDAScalarIntegrand PRIVATE "--expt-relaxed-constexpr"
                                                  "-lineinfo" "-DCURAND"
    )
endif ()
