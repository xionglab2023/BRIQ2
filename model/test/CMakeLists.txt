link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(testSASA TestSSAndSasa.cpp)
add_executable(testBaseLib TestBaseLib.cpp)
add_executable(testBasePairLib TestBasePairLib.cpp)
add_executable(testBBLib TestBackboneRotamer.cpp)
add_executable(testRotLib TestRotLib.cpp)
add_executable(testStacking TestStacking.cpp)
add_executable(testRNASS TestRNASS.cpp)

target_link_libraries (testRNASS briqx)
target_link_libraries (testStacking briqx)
target_link_libraries (testRotLib briqx)
target_link_libraries (testBBLib briqx)
target_link_libraries (testSASA briqx)
target_link_libraries (testBaseLib briqx)
target_link_libraries (testBasePairLib briqx)

install(TARGETS testStacking testBaseLib testBasePairLib testRNASS
    RUNTIME DESTINATION bin)