cmake_minimum_required (VERSION 3.0)
project (BRIQX VERSION 1.0 LANGUAGES CXX)
configure_file(version.h.in "${PROJECT_SOURCE_DIR}/version.h")

set(CMAKE_CXX_STANDARD 20)

SET( EXECUTABLE_OUTPUT_PATH  "${PROJECT_BINARY_DIR}/bin" )
SET( LIBRARY_OUTPUT_PATH  "${PROJECT_BINARY_DIR}/lib" )

include_directories("${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/include/")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++20 -pthread -L ${LIBRARY_OUTPUT_PATH} -Wl,-rpath,$ENV{LIBRARY_PATH}") # bypassing CMake rpath mechanism
# SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++20 -pthread -L ${LIBRARY_OUTPUT_PATH}")

# set(CMAKE_BUILD_RPATH "$ENV{LIBRARY_PATH}")
# set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")  # add internal lib to install rpath
endif("${isSystemDir}" STREQUAL "-1")


if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g -gdwarf-3")
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
    protein/src/*.cpp
    predNA/src/*.cpp
    tools/src/*.cpp
    )
add_library(briqx SHARED ${BR_SRC})
install(TARGETS briqx
        LIBRARY DESTINATION lib)

add_subdirectory (briqxmodule)
add_subdirectory (model/test)
add_subdirectory (forcefield/test)
add_subdirectory (protein/test)
add_subdirectory (predNA/test)
add_subdirectory (motif/test)
add_subdirectory (tools/test)
