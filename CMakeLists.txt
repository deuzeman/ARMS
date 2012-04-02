cmake_minimum_required (VERSION 2.6)
project (ARMS)

set (ARMS_VERSION_MAJOR 2)
set (ARMS_VERSION_MINOR 5)
# set (CMAKE_C_COMPILER icc)
# set (CMAKE_CXX_COMPILER icpc)
# set (CMAKE_CXX_COMPILE_FLAGS -O3 -w0 -march=native -ip ${MPI_COMPILE_FLAGS})
# set (CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})


include_directories ("${PROJECT_SOURCE_DIR}/")
include_directories ("${PROJECT_SOURCE_DIR}/Data")
include_directories ("${PROJECT_SOURCE_DIR}/Eigen")
include_directories ("${PROJECT_SOURCE_DIR}/Minim")
include_directories ("${PROJECT_SOURCE_DIR}/Params")
include_directories ("${PROJECT_SOURCE_DIR}/RanMat")
include_directories ("/usr/src/linux-headers-2.6.38-8/arch/ia64/include")
# include_directories(MPI_INCLUDE_PATH)

add_subdirectory (Data)
add_subdirectory (Minim)
add_subdirectory (Params)
add_subdirectory (Random)
add_subdirectory (RanMat)

add_executable(arms arms.cpp)
target_link_libraries(arms random ranmat params)

add_executable(rmzero rmzero.cpp)
target_link_libraries(rmzero random ranmat params)

add_executable(rmfit rmfit.cpp)
target_link_libraries(rmfit data minim params)

add_executable(rmgof rmgof.cpp)
target_link_libraries(rmgof data ranmat params)
