add_library(arms Comparator Data Minim Log Params Point RanMat Simplex)
target_link_libraries(arms random)

add_executable(rmfit rmfit.cpp)
target_link_libraries(rmfit arms)
