link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(testSASA TestSSAndSasa.cpp)
add_executable(testDM TestDM.cpp)
add_executable(testBaseLib TestBaseLib.cpp)
add_executable(testBBLib TestBackboneRotamer.cpp)
add_executable(testAtLib TestAtomLib.cpp)

target_link_libraries (testAtLib briqx)
target_link_libraries (testBBLib briqx)
target_link_libraries (testSASA briqx)
target_link_libraries (testDM briqx)
target_link_libraries (testBaseLib briqx)