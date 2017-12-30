#!/bin/bash
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2017  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

#
# Creates a source archive from the current git-checkout.
#
# Usage:
#
#  make_archive.sh ARCHIVE_NAME TEMPDIR
#
#  where ARCHIVE_NAME specifies the name of the archive (without extension)
#  and TEMPDIR a temporary directory where packing/unpacking can be done. The script
#  must be invoked from the root of the git repository.
#

archive=$1
if [ -z "$archive" ]; then
  echo "Missing first argument (archive name)" >&2
  exit 1
fi

tmpdir=$2
if [ -z "$tmpdir" ]; then
  echo "Missing second argument (temporary directory)" >&2
  exit 1
fi

if [ ! -d "$tmpdir" ]; then
  echo "Could not create temporary directory $tmpdir"
  exit
fi

echo "Temporary directory: $tmpdir"
echo "Archiving repository with prefix: $archive/"
git archive --format=tar --prefix $archive/ HEAD | tar -C $tmpdir -x -f -
if [ -e .gitmodules ]; then
  paths=$(grep 'path =' .gitmodules | sed -s 's/path = //g')
  for path in $paths; do
    echo "Archiving sub-module with prefix: $archive/$path/"
    cd $path
    git archive --format=tar --prefix $archive/$path/ HEAD | tar -C $tmpdir -x -f -
    cd -
  done
fi
./utils/build/update_release $tmpdir/$archive/RELEASE

echo "Creating archive ${archive}.tar.xz"
tar -C $tmpdir -c -J -f $PWD/${archive}.tar.xz $archive/
