add_custom_target(cython_cluster_abundance
                  COMMAND python setup.py build_ext --inplace
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/cython)
