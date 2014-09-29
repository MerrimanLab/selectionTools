#!/bin/bash
#SBATCH -J OpenMP_JOB
#SBATCH -A          # Project Account
#SBATCH --time=23:00:00     # Walltime
#SBATCH --mem-per-cpu=8196  # memory/cpu (in MB)
#SBATCH --cpus-per-task=8   # 8 OpenMP Threads
#SBATCH -C sb               # sb=Sandybridge,wm=Westmere




module load Python/2.7.4
module load R/3.0.3-goolf-1.5.14

srun python ~/.local/bin/multipop_selection_pipeline -p CEU_ids.txt -p YRI_ids.txt -i CEU_YRI_lactase.vcf --config-file ~/SelectionPipelineTestData/pipeline_test.cfg --a "--phased-vcf" -c 2
