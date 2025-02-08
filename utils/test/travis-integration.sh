#!/usr/bin/env bash
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

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
    ${SOURCE_DIR}/test/src/dftbp/integration/cmake/runtest.sh ${BUILD_DIR}/cmake \
    -DSCALAPACK_LIBRARY="${SCALAPACK_LIBRARY}"

# Integration test for PKG-CONFIG builds
PKG_CONFIG_PATH="${INSTALL_DIR}/lib/pkgconfig:$PKG_CONFIG_PATH" \
    ${SOURCE_DIR}/test/src/dftbp/integration/pkgconfig/runtest.sh ${BUILD_DIR}/pkgconfig
