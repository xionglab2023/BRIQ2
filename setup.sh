#!/bin/bash

if [ "$HOSTNAME"X == MBOOKX ] || [ "$HOSTNAME"X == ABOOKX ]; then
    export CC=/usr/bin/gcc-12
    export CXX=/usr/bin/g++-12
    # cmake -B cmake-build-BRIQ -DCMAKE_BUILD_TYPE=Debug
    # cmake --build cmake-build-BRIQ
elif [ "$HOSTNAME"X == adminX ]; then
    module load gcc/12.3.0
    module load binutils/2.35.2
    module load cmake/3.26.5
    export CC=/public/share/software/gcc/12.3.0/bin/gcc
    export CXX=/public/share/software/gcc/12.3.0/bin/g++
    export BRIQX_DATAPATH=/public/share/pengx_share/briqx/data
fi
echo "compiling on $HOSTNAME"
# cmake -B cmake-build-shared -DCMAKE_BUILD_TYPE=Debug
# cmake --build cmake-build-shared
cmake -B cmake-build-install --install-prefix $HOME/apps/BRIQX
cd cmake-build-install
make -j64 || exit 1
cd ~-
cmake --install cmake-build-install || exit 1
