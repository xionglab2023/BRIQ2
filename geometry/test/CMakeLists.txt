link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(testCs TestCoordSystem.cpp)
add_executable(testTT TestTime.cpp)

target_link_libraries (testCs briqx)
target_link_libraries (testTT briqx)