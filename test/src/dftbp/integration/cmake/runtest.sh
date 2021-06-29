#!/bin/bash
#
# Tests whether the installed DftbPlus library can be used within a CMake project.
#
# Arguments:
#
#   - building directory (will be created, should not exist)
#
# Requirements:
#
#   - Environment variables FC & CC contain the same compiler as used for DftbPlus
#
#   - Environment variable CMAKE_PREFIX_PATH contains the DftbPlus install root.
#
SCRIPTDIR=$(dirname $0)
SCRIPTNAME=$(basename $0)
BUILDDIR=$1
shift

if [ -d ${BUILDDIR} ]; then
  echo "${SCRIPTNAME}: Test build directory '${BUILDDIR}' already exists." >&2
  exit 1
fi

FC=$FC CC=$CC cmake -B ${BUILDDIR} ${SCRIPTDIR} -DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH} "$@" \
  || { echo "Configuration step failed" >&2; exit 1; }
cmake --build ${BUILDDIR} -- VERBOSE=1 || { echo "Build step failed" >&2; exit 1; }
echo "CMake build succeeded!"
