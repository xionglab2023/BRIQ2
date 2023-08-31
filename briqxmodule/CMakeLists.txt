add_executable(bmalign Main.cpp)
target_link_libraries(bmalign briqx)
install(TARGETS bmalign
    RUNTIME DESTINATION bin)
