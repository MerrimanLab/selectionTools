#!/bin/bash

#$1 input file
#$2 chromosome

echo '
#@ shell = /bin/bash
#@ job_name= ${1}
#@ job_type = serial
#@ group = nesi
#@ class = default
#@ wall_clock_limit = 167:00:00
#@ resources = ConsumableMemory(40gb) ConsumableVirtualMemory(40gb)
#@ notification = never
#@ output = ${1}.\${jobid}.out
#@ error = ${1}.\${jobid}.err
#@ parallel_threads = 10
#@queue

ulimit -v 41943040 -m 41943040
shapeit_dir="${HOME}/shapeit.v2.r644.linux.x86_64"

${shapeit_dir}/shapeit.v2.r644.linux.x86_64 --input-ped ${1}_chr${2}_missing_removed_hwe0-001_maf0-05.ped ${1}_chr${2}_missing_removed_hwe0-001_maf0-05.map -M "${shapeit_dir}/genetic_maps_b37/genetic_map_chr${2}_combined_b37.txt" --output-max ${1}_chr${2}_hwe0-001_maf0-05.phased --thread 10


${shapeit_dit}/shapeit --input-ped ${1}.ped ${1}.map -M ${GENETIC_DIR}${2} --output-max ${1} --threads 10
' > ${1}.job
llsubmit ${1}.job

