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
   "-DFYPP_FLAGS='-DTRAVIS'"
   "-DCMAKE_TOOLCHAIN_FILE=../sys/gnu.cmake"
)

mkdir -p _build
pushd _build

if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
   cmake -DCMAKE_BUILD_TYPE=Debug "${cmake_options[@]}" ..
else
   cmake -DCMAKE_BUILD_TYPE=Release "${cmake_options[@]}" ..
fi

make -j 2
ctest -j 2
