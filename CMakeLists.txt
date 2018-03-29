cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(y3_cluster_cpp VERSION 0.1.0 LANGUAGES CXX)
enable_testing()

# These flags are used by everything in the project. Put anything that affects
# the ABI here -- e.g., the standard level.
add_compile_options(-W -Wall) # -Werror)
SET ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math -march=native -mtune=native" )
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

find_package(cuba REQUIRED)

add_subdirectory(test)

