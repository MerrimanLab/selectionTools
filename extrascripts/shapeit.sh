#!/bin/bash

#
# Shape it Load leveler
#


#
# $1 ped
# $2 map
# $3 genetic map
# $4 output
# $5 threads
# $6 wall clock
# $7 memory_total
# 



help(){
cat << EOF
	create_ihh_jobs_nesi.sh 
        arguments are a follows
	    1 = ped
	    2 = map
	    3 = genetic map
	    4 = output
	    5 = threads
        6 = wall clock
        7 = memory total
EOF
}

EXPECTED_ARGS=7

if [ $# -ne $EXPECTED_ARGS ] ; then
    help
    exit 1
fi

PED=$1
MAP=$2
GENETIC_MAP=$3
OUTPUT_FILE=$4
THREADS=$5
WALL_CLOCK=$6
MEMORY_TOTAL=$7

#calculate ulimit
limit=`echo "${MEMORY_TOTAL} * 1024 * 1024" | bc`
#Calculate memory per thread
mem_per_thread=`echo "${MEMORY_TOTAL} / ${THREADS}" | bc`

echo "#@ shell = /bin/bash
#@ environment = COPY_ALL
#@ job_name = shapeit_run
#@ job_type = serial
#@ group = nesi
#@ class = default
#@ notification = never
#@ wall_clock_limit = ${WALL_CLOCK}
#@ resources = ConsumableMemory(${mem_per_thread}gb) ConsumableVirtualMemory(${mem_per_thread}gb)
#@ output = \$(jobid).out 
#@ error = \$(jobid).err
#@ parallel_threads = ${THREADS}
#@ notification = complete
#@ queue

ulimit -v ${limit} -m ${limit}

shapeit --input-ped ${PED} ${MAP} -M ${GENETIC_MAP} --output-max ${OUTPUT_FILE} --thread ${THREADS} " > shapeit.job

sync
llsubmit shapeit.job
















