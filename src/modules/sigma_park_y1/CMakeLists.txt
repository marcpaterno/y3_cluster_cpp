include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

#add_library(summedNumbersCentY1   MODULE summedNumbersCentY1.cc)
#target_link_libraries(summedNumbersCentY1 PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
#target_include_directories(summedNumbersCentY1 PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
#                                                               ${GSL_INCLUDE_DIRS})

#add_library(summedNumbersMiscentY1  MODULE summedNumbersMiscentY1.cc)
#target_link_libraries(summedNumbersMiscentY1  PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
#target_include_directories(summedNumbersMiscentY1  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
#                                                                  ${GSL_INCLUDE_DIRS})

#add_library(avgSigmaCentY1  MODULE avgSigmaCentY1.cc)
#target_link_libraries(avgSigmaCentY1  PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
#target_include_directories(avgSigmaCentY1  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
#                                                                  ${GSL_INCLUDE_DIRS})

#add_library(avgSigmaMiscentY1 MODULE avgSigmaMiscentY1.cc)
#target_link_libraries(avgSigmaMiscentY1  PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
#target_include_directories(avgSigmaMiscentY1  PRIVATE ${PAGANI_DIR} ${CMAKE_SOURCE_DIR}/src
#                                                                     ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})

if (USE_CUDA)
  find_library(NVTX_LIBRARY nvToolsExt PATHS ENV LD_LIBRARY_PATH)

  add_library(summedNumbersCentY1GPU   MODULE summedNumbersCentY1GPU.cu)
  target_link_libraries(summedNumbersCentY1GPU  PRIVATE cosmosis ${NVTX_LIBRARY} utils)
  target_include_directories(summedNumbersCentY1GPU  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                                                     ${PAGANI_DIR})
  set_target_properties(summedNumbersCentY1GPU  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES
                                                                                                ${Y3GCC_TARGET_ARCH})
  target_compile_options(summedNumbersCentY1GPU  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  add_library(summedNumbersMiscentY1GPU   MODULE summedNumbersMiscentY1GPU.cu)
  target_link_libraries(summedNumbersMiscentY1GPU  PRIVATE cosmosis ${NVTX_LIBRARY} utils)
  target_include_directories(summedNumbersMiscentY1GPU  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                                                        ${PAGANI_DIR})
  set_target_properties(summedNumbersMiscentY1GPU  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES
                                                                                                   ${Y3GCC_TARGET_ARCH})
  target_compile_options(summedNumbersMiscentY1GPU  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  add_library(avgSigmaCentY1GPU   MODULE avgSigmaCentY1GPU.cu)
  target_link_libraries(avgSigmaCentY1GPU  PRIVATE cosmosis ${NVTX_LIBRARY} utils)
  set_target_properties(avgSigmaCentY1GPU  PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES
                                                                                                   ${Y3GCC_TARGET_ARCH})
  target_include_directories(avgSigmaCentY1GPU  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                                                        ${PAGANI_DIR})
  target_compile_options(avgSigmaCentY1GPU  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  add_library(avgSigmaMiscentY1GPU   MODULE avgSigmaMiscentY1GPU.cu)
  target_link_libraries(avgSigmaMiscentY1GPU  PRIVATE cosmosis ${NVTX_LIBRARY} utils)
  set_target_properties(avgSigmaMiscentY1GPU  PROPERTIES POSITION_INDEPENDENT_CODE ON
                                                                         CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_include_directories(avgSigmaMiscentY1GPU  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                                                           ${PAGANI_DIR})
  target_compile_options(avgSigmaMiscentY1GPU  PRIVATE "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")

  find_library(NVTX_LIBRARY nvToolsExt PATHS ENV LD_LIBRARY_PATH)

endif ()
