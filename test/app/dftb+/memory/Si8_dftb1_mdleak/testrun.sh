#!/usr/bin/env bash
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2023  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

############################################################################
#
#  Script for starting up a long serial calculation and monitoring
#  memory, generating a synthetic regression tag file which has true
#  values if the calculation 1) stays at constant memory use (up to
#  the granularity reported) 2) responds correctly to a stop file to
#  terminate.
#
############################################################################

DFTBPLUS_CMD=$*
STOP_FILE="stop_driver"
MEMLOG_FILE="log.txt"
OUTPUT_FILE="output"
TAG_FILE="autotest.tag"

# Time in seconds to wait for the code to respond meaningfully
SLEEP_TIME=6
# cycles of sampling the test memory
TEST_CYCLES=10

rm -f $MEMLOG_FILE $STOP_FILE $OUTPUT_FILE $TAG_FILE

# run the actual calculation in the background, keeping the process ID
sh -c "echo $$ > subprocess.pid; exec $DFTBPLUS_CMD > $OUTPUT_FILE" &

for t in $(seq $TEST_CYCLES)
do
    # Wait a little
    sleep $SLEEP_TIME
    ps -p $(cat subprocess.pid) -o %mem | grep -v '%' >> $MEMLOG_FILE
done

# Bring down the calculation cleanly
touch $STOP_FILE
sleep $SLEEP_TIME

echo "output_file         :logical:0:" >> $TAG_FILE
if [ -e $OUTPUT_FILE ]
then
    echo " T" >> $TAG_FILE
else
    echo " F" >> $TAG_FILE
fi

echo "clean_code_halt     :logical:0:" >> $TAG_FILE
if [ $(grep -c 'Molecular dynamics completed' $OUTPUT_FILE) -eq 1 ]
then
    echo " T" >> $TAG_FILE
else
    echo " F" >> $TAG_FILE
fi

echo "memory_sizes        :integer:0:" >> $TAG_FILE
echo -e " " >> $TAG_FILE
echo $(uniq $MEMLOG_FILE | wc -l) >> $TAG_FILE
