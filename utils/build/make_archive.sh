#!/bin/bash
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

#
# Creates a source archive from the current git-checkout.
#
# Usage:
#
#  make_archive.sh ARCHIVE_NAME
#
#  where ARCHIVE_NAME specifies the name of the archive. The script
#  must be invoked from the root of the git repository.
#

archive=$1
if [ -z "$archive" ]; then
  echo "Missing first argument (archive name)" >&2
  exit 1
fi

tmpdir="$(mktemp -d --tmpdir make_archive.XXXXXX)"
if [ -z "$tmpdir" ]; then
  echo "Temporary directory '$tmpdir' could not be created" >&2
  exit 1
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
echo "Deleting temporary directory $tmpdir"
rm -rf $tmpdir
