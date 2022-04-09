cmake_minimum_required(VERSION 3.19 FATAL_ERROR)
project(y3_cluster_cpp VERSION 0.4.0 LANGUAGES CXX C)

if (USE_CUDA)
  enable_language(CUDA)
  message("We will build the CUDA modules")
elseif ()
  message("We will not build the CUDA modules")
endif ()


enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (USE_CUDA)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  if (NOT DEFINED Y3GCC_TARGET_ARCH)
    message(
      FATAL_ERROR
        "You must specify -DY3GCC_TARGET_ARCH=<arch> to generate code for the right GPU e.g. '=70-real' for V100")
  endif ()
endif ()

# These flags are used by everything in the project. Put anything that affects the ABI here -- e.g., the standard level.
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-math-errno")

if (CMAKE_CXX_COMPILER STREQUAL "gcc")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-signaling-nans -fexcess-precision=fast")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")
endif ()


if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(LINKER_OPTIONS "LINKER:-z,defs")
endif()

set(CMAKE_SHARED_MODULE_SUFFIX ".so") # We want this on all operating systems
set(CMAKE_SHARED_MODULE_PREFIX "") # We don't need the leading lib
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON) # for relocatable builds

include(MPITest)

find_package(CUBA REQUIRED)
message("CUBA_LIBRARIES = ${CUBA_LIBRARIES}")

find_package(GSL REQUIRED)
# Fix sometime-broken GSL discovery. Documentation says GSL_INCLUDE_DIRS (note the final 'S') should be set. cmake v10.1
# contains a FindGSL.cmake that says this, but under some conditions does not seem to do so. However, It does set
# 'GSL_INCLUDE_DIR'.
message("GSL_FOUND = ${GSL_FOUND}")
message("GSL_INCLUDE_DIR = ${GSL_INCLUDE_DIR}")
message("GSL_INCLUDE_DIRS = '${GSL_INCLUDE_DIRS}'")
message("GSL_LIBRARIES = ${GSL_LIBRARIES}")
if ("${GSL_INCLUDE_DIRS}" STREQUAL "")
  message("GSL_INCLUDE_DIRS was not set by FindGSL.cmake")
  set(GSL_INCLUDE_DIRS "${GSL_INCLUDE_DIR}")
  message("We have set GSL_INCLUDE_DIRS to ${GSL_INCLUDE_DIRS}")
endif ()

set(MPI_CXX_SKIP_MPICXX)
find_package(MPI REQUIRED COMPONENTS C CXX)

# Make sure we know where CosmoSIS is found
if ("{COSMOSIS_SRC_DIR}" STREQUAL "")
  message(FATAL_ERROR "You must set up CosmoSIS before configuring the build of y3_cluster_cpp")
endif()
set(COSMOSIS_BASE_DIR "$ENV{COSMOSIS_SRC_DIR}/..")

get_property(link_flags TARGET MPI::MPI_C PROPERTY INTERFACE_LINK_OPTIONS)
string(REPLACE "," "$<COMMA>" link_flags "${link_flags}")
set_property(TARGET MPI::MPI_C PROPERTY INTERFACE_LINK_OPTIONS "$<HOST_LINK:${link_flags}>")
message(STATUS      "CMAKE_CXX_IMPLICIT_LINK_LIBRARIES: ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES}")
set(CMAKE_CUDA_IMPLICIT_LINK_LIBRARIES "mpi;${CMAKE_CUDA_IMPLICIT_LINK_LIBRARIES}")
message(STATUS "CMAKE_CUDA_IMPLICIT_LINK_LIBRARIES: ${CMAKE_CUDA_IMPLICIT_LINK_LIBRARIES}")

add_subdirectory(src)
add_subdirectory(test)
