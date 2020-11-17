#!/bin/bash
#
# Tests whether the installed DFTB+ library can be used with pkg-config based builds.
#
# Arguments:
#
#   - building directory (will be created if it does not exist)
#
# Requirements:
#
#   - Environment variables FC and CC contain the same Fortran and C compilers as used for DFTB+
#
#   - Environment variable PKG_CONFIG_PATH contains the lib/pkgconfig folder within
#     the installed DFTB+ tree.
#
#   - You pass all linker options as arguments, which are needed to link an MPI-binary
#     with your compiler. Alternatively, you can specify the name of the MPI-wrappers
#     as your compilers in FC and CC.
#
SCRIPTDIR=$(dirname $0)
SCRIPTNAME=$(basename $0)
BUILDDIR=$1
shift
CUSTOMLIBS=$*

if [ ! -d ${BUILDDIR} ]; then
  mkdir ${BUILDDIR} || { echo "Could not create build dir '${BUILDDIR}'" >&2; exit 1; }
fi

# Make sure, scriptdir is absoulte
cd ${SCRIPTDIR}
SCRIPTDIR=${PWD}
cd -

cd ${BUILDDIR} || { echo "Could not change to build dir '${BUILDDIR}'" >&2; exit 1; }
pkg-config --exists dftbplus || { echo "No PKG-CONFIG found for DFTB+" >&2;  exit 1; }

cflags=$(pkg-config --cflags dftbplus)
libs=$(pkg-config --libs dftbplus)

cmd="${FC} ${cflags} ${SCRIPTDIR}/test_build.f90 ${libs} ${CUSTOMLIBS}"

echo "Build command: ${cmd}"
${cmd} || { echo "Build command failed" >&2;  exit 1; }
echo "PKG-CONFIG build succeeded."
