link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(testKDTree test.cpp)

target_link_libraries (testKDTree briqx)