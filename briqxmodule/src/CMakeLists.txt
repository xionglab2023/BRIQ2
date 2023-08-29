add_library(NSPbm
BriqxModule.cpp
BMAlign.cpp
BMScore.cpp)

link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(bmalign Main.cpp)

target_link_libraries(bmalign onepiece)
