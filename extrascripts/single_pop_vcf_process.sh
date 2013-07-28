#!/bin/bash
#$1 = pop name
#$2 = vcf file of population
#$3 = chr

shapeit_dir="${HOME}/Murray/src/shapeit.v2.r644.linux.x86_64"


vcftools --gzvcf $2 --plink --out ${1}_chr${3}
plink --noweb --file ${1}_chr${3} --geno 0.99 --out ${1}_chr${3}_missing_removed --recode
shapeit.v2.r644.linux.x86_64 --input-ped ${1}_chr${3}_missing_removed.ped ${1}_chr${3}_missing_removed.map -M "${shapeit_dir}/genetic_maps_b37/genetic_map_chr${3}_combined_b37.txt" --output-max ${1}_chr${3}.phased --thread 10
