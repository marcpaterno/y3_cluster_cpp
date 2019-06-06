add_library(utils SHARED
            interp_1d.cc
            interp_2d.cc
            bessel_polynomial_integral.cc
            sin_cos_polynomial_integral.cc)
target_link_libraries(utils PRIVATE ${GSL_LIBRARIES})
target_include_directories(utils PRIVATE ${CUBACPP_DIR} ${EIGEN3_INCLUDE_DIR})
