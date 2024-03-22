link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(testCs TestCoordSystem.cpp)
add_executable(testTT TestTime.cpp)
add_executable(indexToPDB IndexToPDB.cpp)


target_link_libraries (indexToPDB briqx)
target_link_libraries (testCs briqx)
target_link_libraries (testTT briqx)

install(TARGETS indexToPDB
    RUNTIME DESTINATION bin)