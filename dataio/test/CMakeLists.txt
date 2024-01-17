add_executable(txt2binaryMove txt2binaryMove.cpp)
add_executable(txt2binaryBP6DEneTab txt2binaryBP6DEneTab.cpp)
add_executable(txt2binaryBP6DEneTabCG txt2binaryBP6DEneTabCG.cpp)

target_link_libraries(txt2binaryMove briqx)
target_link_libraries(txt2binaryBP6DEneTab briqx)
target_link_libraries(txt2binaryBP6DEneTabCG briqx)

install(TARGETS txt2binaryMove txt2binaryBP6DEneTab txt2binaryBP6DEneTabCG
    RUNTIME DESTINATION bin)