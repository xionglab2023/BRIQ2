link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(rnaGraph rnaGraph.cpp)
add_executable(mtfPDB writeSingleMotifPDB.cpp)
add_executable(testAlign TestAlignment.cpp)
add_executable(simpleAlign TestSimpleAlign.cpp)
add_executable(testAlignScore TestAlignScore.cpp)
add_executable(calAlignScore CalculateAllToAllScore.cpp)


target_link_libraries (calAlignScore loop)
target_link_libraries (testAlignScore loop)
target_link_libraries (testAlign loop)
target_link_libraries (simpleAlign loop)
target_link_libraries (mtfPDB loop)
target_link_libraries (rnaGraph loop)