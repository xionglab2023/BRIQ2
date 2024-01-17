add_executable(bmalign Main.cpp)
add_executable(bmassign bmAssign.cpp)

target_link_libraries(bmalign briqx)
target_link_libraries(bmassign briqx)

install(TARGETS bmalign bmassign
    RUNTIME DESTINATION bin)
