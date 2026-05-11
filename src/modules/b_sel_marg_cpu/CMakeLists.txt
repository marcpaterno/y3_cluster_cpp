include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

# ---------------------------------------------------------------------------
# PAGANI variants (reference, slow -- kept for benchmarking vs the fixed-GL
# evaluator below).
# ---------------------------------------------------------------------------
add_library(BSelMargP1PaganiIntegrand MODULE P1PaganiIntegrand.cc)
target_link_libraries(BSelMargP1PaganiIntegrand PRIVATE cosmosis utils models
                      ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(BSelMargP1PaganiIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})

add_library(BSelMargI1PaganiIntegrand MODULE I1PaganiIntegrand.cc)
target_link_libraries(BSelMargI1PaganiIntegrand PRIVATE cosmosis utils models
                      ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(BSelMargI1PaganiIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})

add_library(BSelMargI2PaganiIntegrand MODULE I2PaganiIntegrand.cc)
target_link_libraries(BSelMargI2PaganiIntegrand PRIVATE cosmosis utils models
                      ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
target_include_directories(BSelMargI2PaganiIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})

# ---------------------------------------------------------------------------
# Fixed-GL evaluator -- single .so that co-computes P1, I1, I2 and writes
# them to sections b_sel_marg_P1 / _I1 / _I2 via
# CosmoSISScalarEvaluatorModule.
# ---------------------------------------------------------------------------
add_library(BSelMargIntegrand MODULE BSelMargIntegrand.cc)
target_link_libraries(BSelMargIntegrand PRIVATE cosmosis utils models
                      ${GSL_LIBRARIES} ${CUBA_LIBRARIES} m)
target_include_directories(BSelMargIntegrand PRIVATE ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})
# Vectorize the MOR Gaussian / exp loops.  -ffast-math (with libmvec
# linked) lets gcc turn scalar std::exp into a vector call; -fopenmp-simd
# enables #pragma omp simd without activating OpenMP threading.
# Scope is local to this module so we don't relax FP semantics globally.
target_compile_options(BSelMargIntegrand PRIVATE
                       -ffast-math -fopenmp-simd -march=native)
