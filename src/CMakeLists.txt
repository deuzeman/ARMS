add_library(arms Comparator Data Minim Params Point RanMat Simplex)
target_link_libraries(arms random)

add_executable(rmfit rmfit.cpp)
target_link_libraries(rmfit data minim params)
