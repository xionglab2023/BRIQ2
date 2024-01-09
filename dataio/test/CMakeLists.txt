add_executable(txt2binaryMove txt2binaryMove.cpp)

target_link_libraries(txt2binaryMove briqx)

install(TARGETS txt2binaryMove
    RUNTIME DESTINATION bin)