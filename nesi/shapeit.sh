#!/bin/bash

#$1 input file
#$2 chromosome

echo "
#@ shell = /bin/bash
#@ job_name= ${1}
#@ job_type = serial
#@ group = nesi
#@ class = default
#@ wall_clock_limit = 167:00:00
#@ resources = ConsumableMemory(40gb) ConsumableVirtualMemory(40gb)
#@ notification = never
#@ output = ${1}.\$(jobid).out
#@ error = ${1}.\$(jobid).err
#@ parralel_threads = 10
#@queue

ulimit -v 41943040 -m 41943040
shapeit --input-ped ${1}.ped ${1}.map -M ${GENETIC_DIR}${2} --output-max ${1} --threads 10
" > ${1}.job
llsubmit ${1}.job

