cmake_minimum_required (VERSION 2.4)
project (ARMS)
# The version number.
set (Tutorial_VERSION_MAJOR 0)
set (Tutorial_VERSION_MINOR 1)

include_directories ("${PROJECT_SOURCE_DIR}/")
include_directories ("${PROJECT_SOURCE_DIR}/Eigen")
include_directories ("${PROJECT_SOURCE_DIR}/RanMat")
include_directories ("${PROJECT_SOURCE_DIR}/Minim")
include_directories ("/usr/src/linux-headers-2.6.38-8/arch/ia64/include")
add_subdirectory (Random) 
add_subdirectory (RanMat)
add_subdirectory (Minim)

# add the executable
add_executable(arms arms.cpp parseInput.cpp)
target_link_libraries(arms random ranmat)
