link_directories(${PROJECT_BINARY_DIR}/lib)

 add_executable(g09ec gaussEnergyCalculator2.cpp)
 add_executable(xtbec xtbEnergyCalculator2.cpp)
 add_executable(xtbec1 xtbEnergyCalculator.cpp)
 add_executable(g09ec1 gaussEnergyCalculator.cpp)
 add_executable(cmSP500 writeMoveSetSP500.cpp)

 target_link_libraries (g09ec briqx)
 target_link_libraries (g09ec1 briqx)
 target_link_libraries (xtbec1 briqx)
 target_link_libraries (xtbec briqx)
 target_link_libraries (cmSP500 briqx)

 install(TARGETS g09ec xtbec cmSP500 xtbec1 g09ec1
     RUNTIME DESTINATION bin)