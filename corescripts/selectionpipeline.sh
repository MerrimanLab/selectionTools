#!/bin/bash

# Murray Cadzow
# July 2013
# University of Otago

# $1 filtered bgzipped vcf file
# $2 chr number
# $3 pop1 text file
# $4 pop2 text file
# $5 name of pop1
# $6 name of pop2


shapeit_dir="${HOME}/Murray/src/shapeit.v2.r644.linux.x86_64"
readCountFilter="${HOME}/Murray/Scripts/readCountFilter/readCountFilter.py"
rehh="./rehh.r"


#remove indels
vcftools --gzvcf $1 --remove-indels --recode
bgzip out.recode.vcf
tabix -p vcf out.recode.vcf.gz
echo point 1

#calculate hapmap  Fst between pops
vcftools --gzvcf out.recode.vcf.gz --hapmap-fst-pop $3 --hapmap-fst-pop $4
echo point 2
#calculate weir and cockerham Fst between pops
vcftools --gzvcf out.recode.vcf.gz --weir-fst-pop $3 --weir-fst-pop $4
echo point 3

#filter vcf to make genotypes called by single read be unknown
#only for our sequenced data
#python "$readCountFilter"  out.recode.vcf > temp.filtered.vcf

#subset out the populations
pop1=`cat $3 | tr '\n' ','`
pop2=`cat $4 | tr '\n' ','`
vcf-subset -c `echo $pop1` out.recode.vcf.gz | bgzip -c > pop1.vcf.gz && tabix -p vcf pop1.vcf.gz
vcf-subset -c `echo $pop2` out.recode.vcf.gz | bgzip -c > pop2.vcf.gz && tabix -p vcf pop2.vcf.gz 
echo point 4

#calculate Tajima's D for each population
vcftools --gzvcf pop1.vcf.gz --TajimaD 1000 --out pop1
vcftools --gzvcf pop2.vcf.gz --TajimaD 1000 --out pop2

#recode vcf to be ped/map
#vcftools --vcf temp.filtered.vcf --plink
vcftools --gzvcf pop1.vcf.gz --plink --out pop1
vcftools --gzvcf pop2.vcf.gz --plink --out pop2
echo point 5

#filter out snps with no genotypes
plink --noweb --file pop1 --geno 0.99 --out pop1_missing_removed --recode
plink --noweb --file pop2 --geno 0.99 --out pop2_missing_removed --recode
echo point 6


#merge with htslib vcfmerge

#phase using shapeit

# rename

#pop1

shapeit --input-ped pop1_missing_removed.ped pop2_missing_removed.map -M "${shapeit_dir}/genetic_maps_b37/genetic_map_chr${2}_combined_b37.txt" --output-max pop1.phased --thread 6
#pop2
shapeit --input-ped pop2_missing_removed.ped pop2_missing_removed.map -M "${shapeit_dir}/genetic_maps_b37/genetic_map_chr${2}_combined_b37.txt" --output-max pop2.phased --thread 6
echo point 7




#rehh
Rscript ./$rehh $5 ./pop1.phased.haps $6 ./pop2.phased.haps $2

echo point 8
echo done
