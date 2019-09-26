include_directories($ENV{COSMOSIS_SRC_DIR} ${GSL_INCLUDE_DIRS})
link_directories($ENV{COSMOSIS_SRC_DIR}/cosmosis/datablock)
add_library(testmain SHARED main.cc)
target_include_directories(testmain
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/externals)

add_executable(make_grid_points_test make_grid_points.test.cc)
target_link_libraries(make_grid_points_test
                      PRIVATE
                      cosmosis
                      testmain)
target_include_directories(make_grid_points_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR}
                       ${CMAKE_SOURCE_DIR}/externals)
add_test(make_grid_points_test make_grid_points_test)

add_executable(make_integration_volume_test make_integration_volume.test.cc)
target_link_libraries(make_integration_volume_test
                      PRIVATE
                      cosmosis
                      testmain)
target_include_directories(make_integration_volume_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR}
                       ${CMAKE_SOURCE_DIR}/externals)
add_test(make_integration_volume_test make_integration_volume_test)


add_executable(integration_test integration.test.cc)
target_link_libraries(integration_test PRIVATE testmain ${CUBA_LIBRARIES})
target_include_directories(integration_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})
add_test(integration_test integration_test)

#add_executable(trivial_gamma_t trivial_gamma_t.cc)
#target_link_libraries(trivial_gamma_t PRIVATE ${CUBA_LIBRARIES} ${GSL_LIBRARIES} profiler utils cfitsio models)
#target_include_directories(trivial_gamma_t
#                           PRIVATE
#                           ${CMAKE_SOURCE_DIR}/src
#                           ${CUBA_INCLUDE_DIR}
#                           ${CUBACPP_DIR}
#                           ${EIGEN3_INCLUDE_DIR})

add_executable(transform_test transform.test.cc)
target_link_libraries(transform_test PRIVATE testmain)
target_include_directories(transform_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(transform_test transform_test)

add_executable(interp_1d_test interp_1d.test.cc)
target_link_libraries(interp_1d_test PRIVATE testmain ${GSL_LIBRARIES} utils)
target_include_directories(interp_1d_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${GSL_INCLUDE_DIR})
add_test(interp_1d_test interp_1d_test)

add_executable(interp_2d_test interp_2d.test.cc)
target_link_libraries(interp_2d_test PRIVATE testmain ${GSL_LIBRARIES} utils)
target_include_directories(interp_2d_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${GSL_INCLUDE_DIR})
add_test(interp_2d_test interp_2d_test)

add_executable(mz_power_law_test mz_power_law.test.cc)
target_link_libraries(mz_power_law_test PRIVATE testmain)
target_include_directories(mz_power_law_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(mz_power_law_test mz_power_law_test)

add_executable(point_3d_test point_3d.test.cc)
target_link_libraries(point_3d_test PRIVATE testmain)
target_include_directories(point_3d_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(point_3d_test point_3d_test)

add_executable(param_space_explorer_test param_space_explorer.test.cc)
target_link_libraries(param_space_explorer_test PRIVATE testmain)
target_include_directories(param_space_explorer_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(param_space_explorer_test param_space_explorer_test)

add_executable(fpsupport_test fpsupport.test.cc)
target_link_libraries(fpsupport_test PRIVATE testmain)
target_include_directories(fpsupport_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(fpsupport_test fpsupport_test)

add_executable(omega_z_sdss_test omega_z_sdss.test.cc)
target_link_libraries(omega_z_sdss_test PRIVATE testmain)
target_include_directories(omega_z_sdss_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(omega_z_sdss_test omega_z_sdss_test)

add_executable(sample_variance_test sample_variance.test.cc)
target_link_libraries(sample_variance_test PRIVATE testmain ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils models)
target_include_directories(sample_variance_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})
#add_test(sample_variance_test sample_variance_test)

#add_executable(hmf_t_test hmf_t.test.cc)
#target_link_libraries(hmf_t_test PRIVATE testmain ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils)
#target_include_directories(hmf_t_test
#                           PRIVATE
#                           ${CMAKE_SOURCE_DIR}/src
#                           ${CMAKE_SOURCE_DIR}/externals)
#add_test(hmf_t_test hmf_t_test)

add_executable(lo_lc_t_test lo_lc_t.test.cc)
target_link_libraries(lo_lc_t_test PRIVATE testmain)
target_include_directories(lo_lc_t_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(lo_lc_t_test lo_lc_t_test)

add_executable(ez_sq_test ez_sq.test.cc)
target_link_libraries(ez_sq_test PRIVATE testmain)
target_include_directories(ez_sq_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(ez_sq_test ez_sq_test)

add_executable(ez_test ez.test.cc)
target_link_libraries(ez_test PRIVATE testmain)
target_include_directories(ez_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(ez_test ez_test)

add_executable(roffset_t_test roffset_t.test.cc)
target_link_libraries(roffset_t_test PRIVATE testmain)
target_include_directories(roffset_t_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(roffset_t_test roffset_t_test)

add_executable(lc_lt_t_test lc_lt_t.test.cc)
target_link_libraries(lc_lt_t_test PRIVATE testmain ${GSL_LIBRARIES} utils models)
target_include_directories(lc_lt_t_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(lc_lt_t_test lc_lt_t_test)

add_executable(lc_lt_probability_test lc_lt_probability_unity.test.cc)
target_link_libraries(lc_lt_probability_test PRIVATE testmain ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils models)
target_include_directories(lc_lt_probability_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})
add_test(lc_lt_probability_test lc_lt_probability_test)

add_executable(dv_do_dz_t_test dv_do_dz_t.test.cc)
target_link_libraries(dv_do_dz_t_test PRIVATE testmain ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils)
target_include_directories(dv_do_dz_t_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBACPP_DIR})
add_test(dv_do_dz_t_test dv_do_dz_t_test)

add_executable(polynomial_test polynomial.test.cc)
target_link_libraries(polynomial_test PRIVATE testmain ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils)
target_include_directories(polynomial_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBACPP_DIR})
add_test(polynomial_test polynomial_test)

add_executable(integer_pow_test integer_pow.test.cc)
target_link_libraries(integer_pow_test PRIVATE testmain)
target_include_directories(integer_pow_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(integer_pow_test integer_pow_test)

add_executable(lo_lc_probability_unity_test lo_lc_probability_unity.test.cc)
target_link_libraries(lo_lc_probability_unity_test PRIVATE testmain ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils)
target_include_directories(lo_lc_probability_unity_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})
add_test(lo_lc_probability_unity_test lo_lc_probability_unity_test)

add_executable(mor_t2_test mor_t2.test.cc)
target_link_libraries(mor_t2_test PRIVATE testmain ${GSL_LIBRARIES} utils)
target_include_directories(mor_t2_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(mor_t2_test mor_t2_test)

add_executable(mor_probability_unity_test mor_probability_unity.test.cc)
target_link_libraries(mor_probability_unity_test PRIVATE testmain ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils)
target_include_directories(mor_probability_unity_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})
add_test(mor_probability_unity_test mor_probability_unity_test)

add_executable(mor2_probability_test mor2_probability_unity.test.cc)
target_link_libraries(mor2_probability_test PRIVATE testmain ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils)
target_include_directories(mor2_probability_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})
add_test(mor2_probability_test mor2_probability_test)

add_executable(roffset_probability_test roffset_probability_unity.test.cc)
target_link_libraries(roffset_probability_test PRIVATE testmain ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils)
target_include_directories(roffset_probability_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})
add_test(roffset_probability_test roffset_probability_test)

add_executable(integrate_delta_sigma integrate_delta_sigma.cc)
target_link_libraries(integrate_delta_sigma PRIVATE ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils models)
target_include_directories(integrate_delta_sigma
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_executable(sigma_crit_inv_test sigma_crit_inv.test.cc)
target_link_libraries(sigma_crit_inv_test PRIVATE testmain utils)
target_include_directories(sigma_crit_inv_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(sigma_crit_inv_test sigma_crit_inv_test)

add_executable(sin_cos_analytic_integral_test sin_cos_analytic_integral.test.cc)
target_link_libraries(sin_cos_analytic_integral_test PRIVATE testmain ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils models)
target_include_directories(sin_cos_analytic_integral_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR})
add_test(sin_cos_analytic_integral_test sin_cos_analytic_integral_test)

add_executable(bessel_analytic_integral_test bessel_analytic_integral.test.cc)
target_link_libraries(bessel_analytic_integral_test PRIVATE testmain ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils models)
target_include_directories(bessel_analytic_integral_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR})
add_test(bessel_analytic_integral_test bessel_analytic_integral_test)

add_executable(bessel_integral_equivalence_test bessel_integral_equivalence.test.cc)
target_link_libraries(bessel_integral_equivalence_test PRIVATE testmain ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils models)
target_include_directories(bessel_integral_equivalence_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR})
add_test(bessel_integral_equivalence_test bessel_integral_equivalence_test)

#add_executable(simultaneous_bin_comparison_test simultaneous_bin_comparison.test.cc)
#target_link_libraries(simultaneous_bin_comparison_test PRIVATE testmain ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils models cfitsio)
#target_include_directories(simultaneous_bin_comparison_test
#                           PRIVATE
#                           ${CMAKE_SOURCE_DIR}/src
#                           ${CMAKE_SOURCE_DIR}/externals
#                           ${CUBA_INCLUDE_DIR}
#                           ${CUBACPP_DIR}
#                           ${EIGEN3_INCLUDE_DIR})

add_executable(integrate_lc_lt integrate_lc_lt.cc)
target_link_libraries(integrate_lc_lt PRIVATE ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils models)
target_include_directories(integrate_lc_lt
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_executable(integrate_lo_lc integrate_lo_lc.cc)
target_link_libraries(integrate_lo_lc PRIVATE ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils)
target_include_directories(integrate_lo_lc
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_executable(integrate_mor integrate_mor.cc)
target_link_libraries(integrate_mor PRIVATE ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils)
target_include_directories(integrate_mor
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_executable(integrate_roffset integrate_roffset.cc)
target_link_libraries(integrate_roffset PRIVATE ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils)
target_include_directories(integrate_roffset
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

#add_executable(simultaneous_bin_comparison simultaneous_bin_comparison.cc)
#target_link_libraries(simultaneous_bin_comparison PRIVATE ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils models cfitsio)
#target_include_directories(simultaneous_bin_comparison
#                           PRIVATE
#                           ${CMAKE_SOURCE_DIR}/src
#                           ${CMAKE_SOURCE_DIR}/externals
#                           ${CUBA_INCLUDE_DIR}
#                           ${CUBACPP_DIR}
#                           ${EIGEN3_INCLUDE_DIR})

add_executable(sample_variance sample_variance.cc)
target_link_libraries(sample_variance PRIVATE ${CUBA_LIBRARIES} ${GSL_LIBRARIES} utils models profiler)
target_include_directories(sample_variance
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals
                           ${CUBA_INCLUDE_DIR}
                           ${CUBACPP_DIR}
                           ${EIGEN3_INCLUDE_DIR})

add_executable(pxizeta_t_test pxizeta_t.test.cc)
target_link_libraries(pxizeta_t_test PRIVATE testmain)
target_include_directories(pxizeta_t_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(pxizeta_t_test pxizeta_t_test)

add_executable(omega_z_y3xspt_test omega_z_y3xspt.test.cc)
target_link_libraries(omega_z_y3xspt_test PRIVATE testmain)
target_include_directories(omega_z_y3xspt_test
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/externals)
add_test(omega_z_y3xspt_test omega_z_y3xspt_test)
