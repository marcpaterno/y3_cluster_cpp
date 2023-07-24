link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

add_library(
  utils SHARED
  bessel_polynomial_integral.cc
  datablock_reader.cc
  debugging.cc
  # fits_loader.cc
  interp_1d.cc
  interp_2d.cc
  make_ifstream.cc
  sin_cos_polynomial_integral.cc
  mpi_support.cc
  )

# If we ever want fits_loader.cc to be linked into the library, we have to get
# the cfitsio library into target_link_libraries. For this to work, we have to
# get a working Find<...>.cmake for this package. So far, no code has ever used
# the functions in this file.

target_link_libraries(utils PRIVATE ${GSL_LIBRARIES} PUBLIC cosmosis MPI::MPI_CXX)
message(STATUS "CUBACPP_DIR is: ${CUBACPP_DIR}")
target_include_directories(
  utils PRIVATE ${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR} ${GSL_INCLUDE_DIRS}
                ${CUBACPP_DIR} ${VEGAS_DIR}
  )
 

if (USE_CUDA)
  add_library(cuda_utils SHARED cuda_interp_1d.cu cuda_interp_2d.cu)
  target_link_libraries(cuda_utils PRIVATE utils cosmosis)
  target_include_directories(cuda_utils PRIVATE ${CMAKE_SRC_DIR}/src ${PAGANI_DIR} ${COSMOSIS_BASE_DIR})
  set_target_properties(cuda_utils PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURE ${Y3GCC_TARGET_ARCH})
  target_compile_options(cuda_utils PRIVATE "--expt-relaxed-constexpr" "-lineinfo")
endif()

