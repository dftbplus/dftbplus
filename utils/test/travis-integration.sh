#!/usr/bin/env bash
#
# Must be called, after DFTB+ had been installed
#
# Expects following env-vars:
#
# FC, CC, ("false"/"true"), SCALAPACK_LIBRARY (when DFTB+ build is MPI-enabled)
# SOURCE_DIR (DFTB+ source dir), BUILD_DIR,
# INSTALL_DIR (DFTB+ install dir with already installed DFTB+)
#
set -ex

# Integration test for CMake builds
CMAKE_PREFIX_PATH="${INSTALL_DIR}:${CMAKE_PREFIX_PATH}" \
    ${SOURCE_DIR}/test/integration/cmake/runtest.sh ${BUILD_DIR}/cmake \
    -DSCALAPACK_LIBRARY="${SCALAPACK_LIBRARY}"

# Integration test for PKG-CONFIG builds
PKG_CONFIG_PATH="${INSTALL_DIR}/lib/pkgconfig:$PKG_CONFIG_PATH" \
    ${SOURCE_DIR}/test/integration/pkgconfig/runtest.sh ${BUILD_DIR}/pkgconfig
