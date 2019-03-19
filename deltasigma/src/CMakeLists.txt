add_library(deltasigma SHARED deltasigma.cc datablock_reader.cc)
target_link_libraries(deltasigma
                      PRIVATE
                      /cosmosis/cosmosis/datablock/libcosmosis.so)
target_include_directories(deltasigma
                           PRIVATE
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/../externals
                           ${CMAKE_SOURCE_DIR}/../src
                           ${CMAKE_SOURCE_DIR}/src
                           ${CMAKE_SOURCE_DIR}/../../../cubacpp)
