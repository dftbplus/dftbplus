#!/usr/bin/env bash
set -ex

if [ "${WITH_MPI}" == "false" ]; then
   WITH_ARPACK=true
else
   WITH_ARPACK=false
fi

if [ "${WITH_ELSI}" == "false" ]; then
   ELSI_VERSION=0
else
   ELSI_VERSION="${elsi_VERSION}"
fi

if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
   BUILD_TYPE="Debug"
else
   BUILD_TYPE="Release"
fi

SOURCE_DIR="${PWD}"
BUILD_DIR="${PWD}/_build"
INSTALL_DIR="${PWD}/_install"

cmake_options=(
   "-DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}"
   "-DWITH_DFTD3=true"
   "-DWITH_MBD=true"
   "-DWITH_TRANSPORT=true"
   "-DWITH_ARPACK=${WITH_ARPACK}"
   "-DWITH_MPI=${WITH_MPI}"
   "-DELSI_VERSION=${ELSI_VERSION}"
   "-DSCALAPACK_LIBRARY='scalapack-openmpi'"
   "-DWITH_API=true"
   "-DFYPP_FLAGS='-DTRAVIS'"
   "-DCMAKE_BUILD_TYPE=${BUILD_TYPE}"
)

cmake -B ${BUILD_DIR}  "${cmake_options[@]}" .
cmake --build ${BUILD_DIR} -- -j
pushd ${BUILD_DIR}
ctest -j
popd
cmake --install ${BUILD_DIR}


CMAKE_PREFIX_PATH="${INSTALL_DIR}:${CMAKE_PREFIX_PATH}" \
    ./test/integration/cmake/runtest.sh _build_cmake \
    -DCMAKE_MODULE_PATH="${SOURCE_DIR}/cmake;${SOURCE_DIR}/external/scalapackfx/origin/cmake" \
    -DSCALAPACK_LIBRARY='scalapack-openmpi'

PKG_CONFIG_PATH="${INSTALL_DIR}/lib/pkgconfig:$PKG_CONFIG_PATH" \
    ./test/integration/pkgconfig/runtest.sh _build_pkgconfig

