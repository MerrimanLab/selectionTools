#!/usr/bin/env python
import sys
import os
import re

from optparse import OptionParser
from pyfasta import Fasta
#
# $1 haps format to be converted to hapmap format
# $2 sample format file
# 


# 11 columns of precursor

# regex for determining we have a valid SNP #


def aa_fasta(options):
    f = Fasta(options.ancestral_fasta)
    chroms = {}
    for key in f.keys():
        if(options.ref_fasta != None):
            chroms[(key.split(' '))[0]] = key
        else:
            chroms[(key.split(':'))[2]] = key
    aaSeq = f[chroms[options.chromosome]]
    return(aaSeq)

def main():
    header = "rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode"
    parser = OptionParser()
    parser.add_option('-i',dest="haps_file",help="Haps Input File")
    parser.add_option('-s',dest="sample_file",help="Sample Input File")
    parser.add_option('-c',dest="chromosome",help="Chromosome")
    parser.add_option('-o',dest="output_file_name",help="Output File name")
    parser.add_option('-a',dest="ancestral_fasta",help="Outgroup fasta file")
    parser.add_option('--id',dest="ancestral_indivdual_id",help="Name of the ancestral Individual")
    parser.add_option('--ref-fasta',dest="ref_fasta",help="Reference fasta")
    (options,args) = parser.parse_args() 
    options.chromosome = str(options.chromosome)
    print(options)
    # Set default ancestral ID#
    if (options.ancestral_indivdual_id is None):
        options.ancestral_indivdual_id = 'ANCESTOR' 
    sample_ids = []
    output = open(options.output_file_name,'w')
    failed_snps = open('failed_snps.txt','w')
    aaSeq=aa_fasta(options)
    with open(options.sample_file,'r') as f:
        for i, line in enumerate(f):
            if(i > 1):
                line = line.split()
                sample_ids.append(line[1])
    # Construct the header line.
    sample_ids.append(options.ancestral_indivdual_id)
    header = header + ' ' + ' '.join(sample_ids) + '\n'
    output.write(header)
    with open(options.haps_file,'r')  as f: 
        for line in f:
            output_line = ''
            line=line.split()
            rsid = line[1]
            pos = line[2]
            ancestral_allele = aaSeq[int(pos)]
            if not (re.match('[ACTGactg]',ancestral_allele)):
               failed_snps.write(rsid + ' ' + pos + '\n')
            else: 
                a1 = line[3]
                a2 = line[4]
                ancestral_genotypes = ancestral_allele + ancestral_allele
                change_alleles=map(lambda x: a1 if(int(x) == 0) else a2,line[5:])
                zipa=change_alleles[0::2]
                zipb=change_alleles[1::2]
                change_alleles = zip(zipa,zipb)
                change_alleles = [''.join(row) for row in change_alleles]
                output_line = rsid + ' ' + a1 + '/' + a2 + ' '+ options.chromosome + ' ' + pos
                output_line = output_line + ' + -9 -9 -9 -9 -9 -9 ' + ' '.join(change_alleles) + ' ' + ancestral_genotypes
                output.write(output_line + '\n')
    output.close() 
    failed_snps.close()

if __name__=="__main__":main()
