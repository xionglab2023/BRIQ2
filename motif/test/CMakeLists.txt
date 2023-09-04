link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(rnaGraph rnaGraph.cpp)
add_executable(mtfPDB writeSingleMotifPDB.cpp)
add_executable(testAlign TestAlignment.cpp)
add_executable(simpleAlign TestSimpleAlign.cpp)
add_executable(testAlignScore TestAlignScore.cpp)
add_executable(calAlignScore CalculateAllToAllScore.cpp)


target_link_libraries (calAlignScore briqx)
target_link_libraries (testAlignScore briqx)
target_link_libraries (testAlign briqx)
target_link_libraries (simpleAlign briqx)
target_link_libraries (mtfPDB briqx)
target_link_libraries (rnaGraph briqx)