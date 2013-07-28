#!/bin/bash
#
#
NUMBER_PER_LINE=`awk '{print NF}' $1 | uniq | wc -l  | cut -d " " -f 1`
REPEAT_POSITIONS=`cat $1 | cut -d " " -f 3 | uniq -d`
if [ "${NUMBER_PER_LINE}" == 1 ] && [ "${REPEAT_POSITIONS}" == ""]; then
    echo "Correct Data File"
else
    echo "wc per line" 
    echo ${NUMBER_PER_LINE}
    echo "repeat positions in haps"
    echo ${REPEAT_POSITIONS}
fi
