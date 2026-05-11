include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

# Merged fixed-GL evaluator: Sigma_prj + DSigma_prj + gamma_t^prj in
# one pass, outer theta with per-R breakpoint log-GL, inner (z, lnM)
# via ring+fg/bg GL.  Publishes {vals, rnd, cl} subfields under
# sigma_prj, dsigma_prj, shear_prj.
add_library(ShearPrjEvaluator MODULE ShearPrjEvaluator.cc)
target_link_libraries(ShearPrjEvaluator PRIVATE cosmosis utils models
                      ${GSL_LIBRARIES} ${CUBA_LIBRARIES} fmt::fmt)
target_include_directories(ShearPrjEvaluator PRIVATE ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})
# -ffast-math lets gcc vectorise the scalar std::sqrt/std::exp in the
# hot (zt, theta) loops; -fopenmp-simd enables #pragma omp simd
# without turning on OpenMP threading.  Scope is local so global FP
# semantics stay strict.
target_compile_options(ShearPrjEvaluator PRIVATE
                       -ffast-math -fopenmp-simd -march=native)

# Adaptive reference path: same theta outer (log-GL with breakpoints)
# but inner (z, lnM) driven by cubacpp::Cuhre / cubacpp::Vegas.  Used
# for head-to-head convergence diffs against ShearPrjEvaluator and
# Python SigmaPrj.  Writes sigma_prj_cuhre / dsigma_prj_cuhre /
# shear_prj_cuhre.
add_library(ShearPrjCuhre MODULE ShearPrjCuhre.cc)
target_link_libraries(ShearPrjCuhre PRIVATE cosmosis utils models
                      ${GSL_LIBRARIES} ${CUBA_LIBRARIES} fmt::fmt)
target_include_directories(ShearPrjCuhre PRIVATE ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})

# GSL-driven variant: QAGP over z (explicit breakpoint at z=zob),
# fixed GL on (lnM, theta).  Legacy -- still on the old theta recipe
# (2-segment split at theta_R).  Kept as a z-adaptive diagnostic.
add_library(ShearPrjGsl MODULE ShearPrjGsl.cc)
target_link_libraries(ShearPrjGsl PRIVATE cosmosis utils models
                      ${GSL_LIBRARIES} ${CUBA_LIBRARIES} fmt::fmt)
target_include_directories(ShearPrjGsl PRIVATE ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})
target_compile_options(ShearPrjGsl PRIVATE
                       -ffast-math -fopenmp-simd -march=native)
