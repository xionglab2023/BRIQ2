#!/bin/bash

if [ "$HOSTNAME"X == MBOOKX ]; then
    echo "compiling on $HOSTNAME"
    export CC=/usr/bin/gcc-12
    export CXX=/usr/bin/g++-12
    # cmake -B cmake-build-BRIQ -DCMAKE_BUILD_TYPE=Debug
    # cmake --build cmake-build-BRIQ
    cmake -B cmake-build -DCMAKE_BUILD_TYPE=Debug
    cmake --build cmake-build
fi
