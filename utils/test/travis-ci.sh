#!/usr/bin/env bash
set -ex

if [ "${WITH_MPI}" == "false" ]; then
   WITH_ARPACK=true
else
   WITH_ARPACK=false
fi

cmake_options=(
   "-DWITH_DFTD3=true"
   "-DWITH_TRANSPORT=true"
   "-DWITH_ARPACK=${WITH_ARPACK}"
   "-DWITH_MPI=${WITH_MPI}"
   "-DWITH_API=true"
   "-DFYPP_FLAGS='-DTRAVIS'"
)

if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
   BUILD_TYPE="Debug"
else
   BUILD_TYPE="Release"
fi

mkdir -p _build
pushd _build
cmake -DCMAKE_BUILD_TYPE=${BUILD_TYPE} "${cmake_options[@]}" ..
make -j 2
ctest -j 2
