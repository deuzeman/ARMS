cmake_minimum_required (VERSION 2.6)
project (ARMS)

set (ARMS_VERSION_MAJOR 4)
set (ARMS_VERSION_MINOR 0)
set (CMAKE_CXX_COMPILER mpicxx)
set (CMAKE_CXX_COMPILE_FLAGS -O3 -march=native)

include_directories ("${PROJECT_SOURCE_DIR}")
include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("/usr/src/linux-headers-2.6.38-8/arch/ia64/include")

add_subdirectory (Random)
add_subdirectory (src)

add_executable(rms rms.cpp)
target_link_libraries(rms arms)

add_executable(rmfit rmfit.cpp)
target_link_libraries(rmfit arms)
