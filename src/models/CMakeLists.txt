link_directories($ENV{COSMOSIS_SRC_DIR}/cosmosis/datablock)

add_library(models SHARED
            sample_variance.cc
            lc_lt_t.cc
            lc_lt_y1_t.cc
            int_lc_lt_des_t.cc)


target_link_libraries(models
                      PRIVATE
                      utils
                      cosmosis
                      ${GSL_LIBRARIES}
                      ${CUBA_LIBRARIES})

target_include_directories(models
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           $ENV{COSMOSIS_SRC_DIR}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR}
                           ${GSL_INCLUDE_DIR})

if (USE_CUDA)
add_library(cuda_models SHARED
            int_lc_lt_des_t.cu)

target_link_libraries(cuda_models
            PRIVATE
            cosmosis
            utils)
set_target_properties(cuda_models PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  CUDA_ARCHITECTURES 70-real)
target_include_directories(cuda_models PRIVATE
  ${CMAKE_SOURCE_DIR}/src
  $ENV{COSMOSIS_SRC_DIR}
  ${PAGANI_DIR})
target_compile_options(cuda_models PRIVATE "--expt-relaxed-constexpr")
endif()