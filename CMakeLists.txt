cmake_minimum_required (VERSION 2.6)
project (ARMS)

set (ARMS_VERSION_MAJOR 2)
set (ARMS_VERSION_MINOR 1)

include_directories ("${PROJECT_SOURCE_DIR}/")
include_directories ("${PROJECT_SOURCE_DIR}/Data")
include_directories ("${PROJECT_SOURCE_DIR}/Eigen")
include_directories ("${PROJECT_SOURCE_DIR}/Minim")
include_directories ("${PROJECT_SOURCE_DIR}/Params")
include_directories ("${PROJECT_SOURCE_DIR}/RanMat")
include_directories ("/usr/src/linux-headers-2.6.38-8/arch/ia64/include")

add_subdirectory (Data)
add_subdirectory (Minim)
add_subdirectory (Params)
add_subdirectory (Random)
add_subdirectory (RanMat)

add_executable(arms arms.cpp)
target_link_libraries(arms random ranmat params)

add_executable(rmfit rmfit.cpp)
target_link_libraries(rmfit data minim)

