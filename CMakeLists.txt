cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(y3_cluster_cpp VERSION 0.2.0 LANGUAGES CXX)
enable_testing()

# These flags are used by everything in the project. Put anything that affects
# the ABI here -- e.g., the standard level.
add_compile_options(-W -Wall) # -Werror)
SET ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-math-errno  -fno-signaling-nans -fexcess-precision=fast -march=native -mtune=native" )
#SET ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-math-errno  -fno-signaling-nans -fexcess-precision=fast -march=native -mtune=native -fsanitize=address" )
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

find_package(cuba REQUIRED)
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)



add_subdirectory(src)
add_subdirectory(test)
