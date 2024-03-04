link_directories(${PROJECT_BINARY_DIR}/lib)

 add_executable(g09ec gaussEnergyCalculator2.cpp)
 add_executable(xtbec xtbEnergyCalculator2.cpp)
 add_executable(cmSP500 writeMoveSetSP500.cpp)

 target_link_libraries (g09ec briqx)
 target_link_libraries (xtbec briqx)
 target_link_libraries (cmSP500 briqx)

 install(TARGETS g09ec xtbec cmSP500
     RUNTIME DESTINATION bin)