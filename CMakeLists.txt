cmake_minimum_required (VERSION 3.0)
project (BMAlign VERSION 0.0.1 LANGUAGES CXX)
configure_file(version.h.in "${PROJECT_SOURCE_DIR}/version.h")

set(CMAKE_CXX_STANDARD 20)

SET( EXECUTABLE_OUTPUT_PATH  "${PROJECT_BINARY_DIR}/bin" )
SET( LIBRARY_OUTPUT_PATH  "${PROJECT_BINARY_DIR}/lib" )

include_directories("${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/include/")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++20 -pthread -L ${LIBRARY_OUTPUT_PATH}" )

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

file(GLOB_RECURSE BR_SRC
    briqxmodule/src/*.cpp
    dataio/src/*.cpp
    forcefield/src/*.cpp
    geometry/src/*.cpp
    math/src/*.cpp
    model/src/*.cpp
    motif/src/*.cpp
    para/src/*.cpp
    predNA/src/*.cpp
    tools/src/*.cpp
    )
add_library(briqx ${BR_SRC})

add_executable(bmalign briqxmodule/Main.cpp)
target_link_libraries(bmalign briqx)