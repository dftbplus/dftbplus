#!/usr/bin/env bash
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

#
# Collect the output of failed tests into a tarball
#
# Usage: collect_failed_output.sh <build_dir> <archive_name>
#
# build_dir: The build directory where code was built and tested with CMake
# archive_name: The name of the tarball to create
#
# The program assumes that the test names in the failure log are of the form *:<test_name>
# <test_name> is the name of the test directory that is below <build_dir>/test.
#

BUILD_DIR="$1"
ARCHIVE_NAME="$2"

if [ "$#" -ne "2" ]; then
    echo "Usage: $0 <build_dir> <archive_name>"
    exit 1
fi

if [ ! -d ${BUILD_DIR} ]; then
    echo "No build directory found at ${BUILD_DIR}"
    exit 0
fi

TEST_FAILURE_LOG="${BUILD_DIR}/Testing/Temporary/LastTestsFailed.log"

if [ ! -f ${TEST_FAILURE_LOG} ]; then
    echo "No test failure log found at ${TEST_FAILURE_LOG}"
    exit 0
fi

TEST_ROOT_DIR="test"

echo "Checking ${TEST_FAILURE_LOG} for failed tests"
FAILED_TEST_DIRS=()
for test in $(cat ${TEST_FAILURE_LOG}); do
  dir="${TEST_ROOT_DIR}/${test#*:}"
  if [ -e ${BUILD_DIR}/${dir} ]; then
    FAILED_TEST_DIRS+=(${dir})
  else
    echo "Ignoring non-existent test directory ${dir}"
  fi
done

if [ ${#FAILED_TEST_DIRS[@]} -eq 0 ]; then
    echo "No failed test directories to archive, exiting..."
    exit 0
fi

tar -C ${BUILD_DIR} -c -f ${ARCHIVE_NAME} ${FAILED_TEST_DIRS[*]}

echo "Archive with failed test directories created at ${ARCHIVE_NAME}"
