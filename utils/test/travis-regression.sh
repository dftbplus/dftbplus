#!/usr/bin/env bash
#
# Expects following env-vars:
#
# FC, CC, WITH_MPI ("false"/"true"), SCALAPACK_LIBRARY (for MPI-builds)
# SOURCE_DIR, BUILD_DIR, INSTALL_DIR
#
set -ex

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

cmake_options=(
  "-DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}"
  "-DWITH_DFTD3=true"
  "-DWITH_MBD=true"
  "-DWITH_TRANSPORT=true"
  "-DWITH_ARPACK=${WITH_ARPACK}"
  "-DWITH_MPI=${WITH_MPI}"
  "-DELSI_VERSION=${ELSI_VERSION}"
  "-DSCALAPACK_LIBRARY='${SCALAPACK_LIBRARY}'"
  "-DWITH_API=true"
  "-DFYPP_FLAGS='-DTRAVIS'"
  "-DHYBRID_CONFIG_METHODS='Submodule'"
  "-DCMAKE_BUILD_TYPE=${BUILD_TYPE}"
)

cmake -B ${BUILD_DIR} "${cmake_options[@]}" ${SOURCE_DIR}
cmake --build ${BUILD_DIR} -- -j
pushd ${BUILD_DIR}
ctest -j --output-on-failure
popd
cmake --install ${BUILD_DIR}
