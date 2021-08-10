link_directories($ENV{COSMOSIS_SRC_DIR}/cosmosis/datablock)

add_library(utils SHARED
            bessel_polynomial_integral.cc
            datablock_reader.cc
            debugging.cc
            #fits_loader.cc
            interp_1d.cc
            interp_2d.cc
            make_ifstream.cc
            sin_cos_polynomial_integral.cc)

# If we ever want fits_loader.cc to be linked into the library, we have to get
# the cfitsio library into target_link_libraries. For this to work, we have to
# get a working Find<...>.cmake for this package. So far, no code has ever used
# the functions in this file.

target_link_libraries(utils PRIVATE ${GSL_LIBRARIES} cosmosis)
message(STATUS "CUBACPP_DIR is: ${CUBACPP_DIR}")
target_include_directories(utils PRIVATE
                           $ENV{COSMOSIS_SRC_DIR}
                           ${CUBA_INCLUDE_DIR}
                           ${GSL_INCLUDE_DIRS}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})
