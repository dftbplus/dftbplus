#!/bin/bash

PYTHON_TESTS="grids repeatgen straingen xyz2gen gen2xyz gen2cif dp_bands dp_dos"

function abspath() {
  cd $1
  echo $PWD
}

testdir=$(dirname $0)
testdir=$(abspath $testdir)

packagedir="$testdir/../../../tools/dptools"
packagedir=$(abspath $packagedir)

workdir="./"
workdir=$(abspath $workdir)

if [ $# -gt 0 ]; then
  pythons=$*
else
  pythons="python3"
fi

pythonpath="$packagedir/src"
if [ -z "$PYTHONPATH" ]; then
  export PYTHONPATH="$pythonpath"
else
  export PYTHONPATH="$pythonpath:$PYTHONPATH"
fi

echo $PYTHONPATH
cd $workdir
failed="0"
failing_tests=""
for python in $pythons; do
  echo -e "\n* Testing with interpreter $python"
  for python_test in $PYTHON_TESTS; do
    echo -e "\n** Executing $python_test tests"
    $python $testdir/test_${python_test}.py
    exitcode=$?
    if [ $exitcode != 0 ]; then
      failed="$(($failed + 1))"
      if [ -z "$failing_tests" ]; then
	failing_tests="$python_test ($python)"
      else
	failing_tests="$failing_tests, $python_test ($python)"
      fi
    fi
  done
done

echo
echo "********************************************************************************"
if [ $failed -gt 0 ]; then
  echo "Failing test runs: $failed" >&2
  echo "Failing test(s): $failing_tests" >&2
  exit 1
else
  echo "All test runs finished successfully"
  exit 0
fi
