#!/usr/bin/env bash
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#


refdir=$1
curdir=$2
shift 2
tests="$*"

if [ refdir = '-h' -o -z "$refdir" -o -z "$curdir" -o -z "$tests" ]; then
  echo "$0 refdir curdir test [test [test ...]]"
  echo ""
  echo "Compares the text output files between tests in refdir and curdir."
  echo ""
  echo "refdir:  Root of the reference autotest test calculations"
  echo "curdir: Root of the current autotest test calculations"
  echo "test: Test dir to compare (relative to refdir and curdir)"
  exit
fi

files="output detailed.out band.out md.out hessian.out results.tag detailed.xml autotest.tag geo_end.xyz geo_end.gen"

for test in $tests; do
  echo "TEST: $test"
  for file in $files; do
    if [ -f $curdir/$test/$file ]; then
      echo "checking: $test/$file"
      mydiff=$(diff -w -B -b $curdir/$test/$file $refdir/$test/$file)
      if [ -n "$mydiff" ]; then
        echo "*** DIFF in $test/$file"
        echo "$mydiff"
        echo
      fi
    else
      echo "missing: $test/$file"
    fi
  done
  echo "----------------------------------------------------------------------"
done

