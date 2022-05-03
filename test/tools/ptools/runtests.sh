#!/usr/bin/env bash

PYTHON_TESTS="detailedout eigenvecout hsd2geo hsdinput resultstag"

function abspath() {
  cd $1
  echo $PWD
}

testdir=$(dirname $0)
testdir=$(abspath ${testdir})

pythonpath=$(abspath "$testdir/../../../tools/dptools/src")
pythonpath2=$(abspath "$testdir/../../../tools/python/ptools/src")
pythonpath3=$(abspath "$testdir/../")
pythonpath="${pythonpath}:${pythonpath2}:${pythonpath3}"


workdir="./"
workdir=$(abspath $workdir)

if [ $# -gt 0 ]; then
  pythons=$*
else
  pythons="python3"
fi

if [ -z "$PYTHONPATH" ]; then
  PYTHONPATH="$pythonpath"
  export PYTHONPATH
else
  PYTHONPATH="$pythonpath:$PYTHONPATH"
  export PYTHONPATH
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
