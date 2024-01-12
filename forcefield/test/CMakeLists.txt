link_directories(${PROJECT_BINARY_DIR}/lib)

# add_executable(testBBCG TestBackboneCG.cpp)
add_executable(testClashEnergy TestAtomicEnergy.cpp)
add_executable(rnaEne CalRNAEnergy.cpp)
add_executable(testBaseAtomName TestBaseAtomNameAndSep.cpp)
add_executable(basePairEnergy TestBasePairEnergy.cpp)
add_executable(testBP6DEneTabIO TestBP6DEneTabIO.cpp)
add_executable(testBP6DEneTabCGIO TestBP6DEneTabCGIO.cpp)

# target_link_libraries (testBBCG briqx)
target_link_libraries (basePairEnergy briqx)
target_link_libraries (testClashEnergy briqx)
target_link_libraries (testBaseAtomName briqx)
target_link_libraries (rnaEne briqx)
target_link_libraries (testBP6DEneTabIO briqx)
target_link_libraries (testBP6DEneTabCGIO briqx)

# install(TARGETS rnaEne testBBCG
install(TARGETS rnaEne testBP6DEneTabIO testBP6DEneTabCGIO
    RUNTIME DESTINATION bin)