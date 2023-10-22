#!/usr/bin/env bash
RUN_CMD="$@"
echo $RUN_CMD > stdout.log
exec $RUN_CMD ../../testers/test_mpiexpdmhsz | tee -a stdout.log
