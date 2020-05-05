link_directories($ENV{COSMOSIS_SRC_DIR}/cosmosis/datablock)

add_library(utils SHARED
            bessel_polynomial_integral.cc
            datablock_reader.cc
            interp_1d.cc
            interp_2d.cc
            make_ifstream.cc
            sin_cos_polynomial_integral.cc)
target_link_libraries(utils PRIVATE ${GSL_LIBRARIES} cosmosis)
message(STATUS "CUBACPP_DIR is: ${CUBACPP_DIR}")
target_include_directories(utils PRIVATE
                           $ENV{COSMOSIS_SRC_DIR}
                           ${CUBA_INCLUDE_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})
