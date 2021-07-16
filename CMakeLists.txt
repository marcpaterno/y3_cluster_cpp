cmake_minimum_required(VERSION 3.19 FATAL_ERROR)
project(y3_cluster_cpp VERSION 0.2.0 LANGUAGES CXX CUDA)
enable_testing()

# These flags are used by everything in the project. Put anything that affects
# the ABI here -- e.g., the standard level.
#add_compile_options(-W -Wall -Wextra -Werror)
#add_compile_options(-Wnon-virtual-dtor -Wcast-align -Woverloaded-virtual -Wpedantic -Wdouble-promotion)
#SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-math-errno -march=native -mtune=native" )
SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-math-errno")

if (CMAKE_CXX_COMPILER STREQUAL "gcc")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-signaling-nans -fexcess-precision=fast")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")
endif()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_SHARED_MODULE_SUFFIX ".so") # We want this on all operating systems
set(CMAKE_SHARED_MODULE_PREFIX "")    # We don't need the leading lib
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON) # for relocatable builds

find_package(CUBA REQUIRED)
message("CUBA_LIBRARIES = ${CUBA_LIBRARIES}")

find_package(GSL REQUIRED)
# Fix sometime-broken GSL discovery.
# Documentation says GSL_INCLUDE_DIRS (not the final 'S') should be set.
# cmake v10.1 contains a FindGSL.cmake that says this, but under some conditions
# does not seem to do so. However, It does set 'GSL_INCLUDE_DIR'.
message("GSL_FOUND = ${GSL_FOUND}")
message("GSL_INCLUDE_DIR = ${GSL_INCLUDE_DIR}")
message("GSL_INCLUDE_DIRS = '${GSL_INCLUDE_DIRS}'")
message("GSL_LIBRARIES = ${GSL_LIBRARIES}")
if ("${GSL_INCLUDE_DIRS}" STREQUAL "")
  message("GSL_INCLUDE_DIRS was not set by FindGSL.cmake")
  set(GSL_INCLUDE_DIRS "${GSL_INCLUDE_DIR}")
  message("We have set GSL_INCLUDE_DIRS to ${GSL_INCLUDE_DIRS}")
endif()

add_subdirectory(src)
add_subdirectory(test)
