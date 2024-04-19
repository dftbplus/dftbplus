#!/usr/bin/env bash
#
# Build external components first and uses them as external dependencies during build
#
# Expects following env-vars:
#
# WITH_MPI ("false"/"true"), SCALAPACK_LIBRARY, SOURCE_DIR, BUILD_DIR, INSTALL_DIR
#
set -ex

if [ "${WITH_MPI}" == "true" ]; then

  cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
        -DBUILD_EXPORTED_TARGETS_ONLY=True -DCMAKE_BUILD_TYPE=Debug \
        -B ${BUILD_DIR}/mpifx ${SOURCE_DIR}/external/mpifx/origin
  cmake --build ${BUILD_DIR}/mpifx -- -j
  cmake --install ${BUILD_DIR}/mpifx

  cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
        -DSCALAPACK_LIBRARY="${SCALAPACK_LIBRARY}" -DBUILD_EXPORTED_TARGETS_ONLY=True\
        -DCMAKE_BUILD_TYPE=Debug \
        -B ${BUILD_DIR}/scalapackfx ${SOURCE_DIR}/external/scalapackfx/origin
  cmake --build ${BUILD_DIR}/scalapackfx -- -j
  cmake --install ${BUILD_DIR}/scalapackfx

  cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
        -DENABLE_SCALAPACK_MPI=True \
        -DSCALAPACK_LIBRARIES="${SCALAPACK_LIBRARY}" \
        -DCMAKE_BUILD_TYPE=Debug \
        -DBUILD_SHARED_LIBS=False \
        -DBUILD_TESTING=False \
        -B ${BUILD_DIR}/mbd ${SOURCE_DIR}/external/mbd/origin
  cmake --build ${BUILD_DIR}/mbd -- -j
  cmake --install ${BUILD_DIR}/mbd

  cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
        -DWITH_MPI=True \
        -DCMAKE_BUILD_TYPE=Debug \
        -DBUILD_TESTING=False \
        -B ${BUILD_DIR}/negf ${SOURCE_DIR}/external/libnegf/origin
  cmake --build ${BUILD_DIR}/negf -- -j
  cmake --install ${BUILD_DIR}/negf

  cmake -DWITH_MPI=True -DWITH_MBD=True -DWITH_TRANSPORT=True  \
        -DHYBRID_CONFIG_METHODS='Find' \
        -DSCALAPACK_LIBRARY="${SCALAPACK_LIBRARY}" \
        -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
        -DCMAKE_BUILD_TYPE=Debug \
        -B ${BUILD_DIR}/dftbplus ${SOURCE_DIR}
  cmake --build ${BUILD_DIR}/dftbplus -- VERBOSE=1 dftb+

else

  cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
        -DENABLE_SCALAPACK_MPI=False \
        -DCMAKE_BUILD_TYPE=Debug \
        -DBUILD_TESTING=False \
        -DBUILD_SHARED_LIBS=False \
        -B ${BUILD_DIR}/mbd ${SOURCE_DIR}/external/mbd/origin
  cmake --build ${BUILD_DIR}/mbd -- -j
  cmake --install ${BUILD_DIR}/mbd
  
  cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
        -DWITH_MPI=False \
        -DCMAKE_BUILD_TYPE=Debug \
        -DBUILD_TESTING=False \
        -B ${BUILD_DIR}/negf ${SOURCE_DIR}/external/libnegf/origin
  cmake --build ${BUILD_DIR}/negf -- -j
  cmake --install ${BUILD_DIR}/negf

  cmake -DWITH_MPI=False -DWITH_MBD=True -DWITH_TRANSPORT=True  \
        -DHYBRID_CONFIG_METHODS='Find' \
        -DBUILD_SHARED_LIBS=False \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
        -B ${BUILD_DIR}/dftbplus ${SOURCE_DIR}
  cmake --build ${BUILD_DIR}/dftbplus -- VERBOSE=1 dftb+

fi
