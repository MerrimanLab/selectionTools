#!/bin/bash

#
# $1 large_input_file
# $2 window
# $3 overlap
# $4 parralel cores
# $5 Chromosome
#
# The overlap is the overlap for each chunk in this case 10* the overlap for
# the job run
# 
# 

help(){
cat << EOF
	create_ihh_jobs_nesi.sh 
        arguments are a follows
	1 = input_haps_file
	2 = Window size
	3 = overlap between all the chunks
	4 = parrallel cores
	5 = Chromosome 
    6 = population
    7 = MAF filter
    8 = Wall Clock hh:mm:ss
    9 = Memory per thread (gb)
EOF
}

HAPS=$1
WINDOW=$2
OVERLAP=$3
PARRALEL_CORES=$4
CHROM=$5
POP=$6
MAF=$7
WALL_CLOCK=$8
MEM_PER_THREAD=$9



if [ "$1" == "" ] ; then
	help
	exit 1
fi
echo "$@" > command.out
bigWindow=`echo "(${WINDOW}-${OVERLAP}) * (${PARRALEL_CORES}) + ${OVERLAP}" | bc`
echo $bigWindow
max_line=`tail -1 $HAPS | awk '{ print $3 }'`
#let "limit = 30000 * 1024"
mem_in_gigs=`echo "${MEM_PER_THREAD} * ${PARRALEL_CORES}" | bc`
limit=`echo "${MEM_PER_THREAD} * ${PARRALEL_CORES} * 1024 * 1024" | bc`


noFolders=`echo "(${max_line}+${OVERLAP})/(${bigWindow}-${OVERLAP}) + 1" | bc`
echo $noFolders 
python prepare_files_aa.py $HAPS $bigWindow $OVERLAP $POP
for i in $(eval echo "{1..${noFolders}}") ; do
     let "offset = ${i} * 10"
     offset=`echo "${offset} - 9" | bc`
     echo "${i}" >> folderlist
     echo "Processing $i in $working_dir"
     echo "#@ shell = /bin/bash
     #@ environment = COPY_ALL
     #@ job_name = ihs_CEU_${i}
     #@ job_type = serial
     #@ group = nesi
     #@ class = default
     #@ notification = never
     #@ wall_clock_limit = 1:00:00
     #@ resources = ConsumableMemory(${mem_in_gigs}) ConsumableVirtualMemory(${mem_in_gigs})
     #@ output = ${i}/\$(jobid).out
     #@ error = ${i}/\$(jobid).err
     #@ parallel_threads =${PARRALEL_CORES} 
     #@ notification = complete
     #@ queue
 		 
     module load R/3.0.1
     # ulimit sets memory constraints for jobs running on single nodes (to prevent the job
     # from consuming too much memory).
     # The first argument is in KB and should equal ConsumableMemory.
     ulimit -v ${limit} -m ${limit}
     mkdir $i
     # Call R with the input file as a command line argument
     Rscript multicore_iHH.R ${POP} ${POP}${i}.phaps ${CHROM} ${WINDOW} ${OVERLAP} ${PARRALEL_CORES} $i $offset ${MAF}
	" > ${i}.job
     sync
     llsubmit ${i}.job
 
     #remove the temp file
     #rm ${i}.job	
done

NUMBER_FINISHED=0

while true; do 
    NUMBER_FINISHED=`ls */*ihh | wc -l`
    if [ "${NUMBER_FINISHED}" == "${noFolders}" ]; then
        Rscript nesi_recombine_ihh.R ${POP} ${CHROM} ${WINDOW} ${OVERLAP} ${OVERLAP} ${PARRALEL_CORES}
        break
    fi
    sleep 5m
done

