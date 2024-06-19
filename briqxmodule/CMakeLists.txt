add_executable(bmalign Main.cpp)
add_executable(bmassign bmAssign.cpp)
add_executable(getMotifPDB getMotifPDBbyBM.cpp)
add_executable(findRedundantPDB findRedundantPDB.cpp)

target_link_libraries(bmalign briqx)
target_link_libraries(bmassign briqx)
target_link_libraries(getMotifPDB briqx)
target_link_libraries(findRedundantPDB briqx)

install(TARGETS bmalign bmassign getMotifPDB findRedundantPDB
    RUNTIME DESTINATION bin)
