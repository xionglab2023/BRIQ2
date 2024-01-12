add_executable(txt2binaryMove txt2binaryMove.cpp)
add_executable(txt2binaryBP6DEneTab txt2binaryBP6DEneTab.cpp)

target_link_libraries(txt2binaryMove briqx)
target_link_libraries(txt2binaryBP6DEneTab briqx)

install(TARGETS txt2binaryMove txt2binaryBP6DEneTab
    RUNTIME DESTINATION bin)