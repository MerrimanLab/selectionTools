import sys
import os

from optparse import OptionParser

#
# $1 haps format to be converted to hapmap format
# $2 sample format file
# 


# 11 columns of precursor




def main():
    header = "rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode"
    parser = OptionParser()
    parser.add_option('-i',dest="haps_file",help="Haps Input File")
    parser.add_option('-s',dest="sample_file",help="Sample Input File")
    parser.add_option('-c',dest="chromosome",help="Chromosome")
    parser.add_option('-o',dest="output_file_name",help="Output File name")
    (options,args) = parser.parse_args() 
    options.chromosome = str(options.chromosome)
    sample_ids = []
    output = open(options.output_file_name,'w') 
    with open(options.sample_file,'r') as f:
        for i, line in enumerate(f):
            if(i > 1):
                line = line.split()
                sample_ids.append(line[1])
    # Construct the header line.
    header = header + ' ' + ' '.join(sample_ids) + '\n'
    output.write(header)
    with open(options.haps_file,'r')  as f: 
        for line in f:
            output_line = ''
            line=line.split()
            rsid = line[1]
            pos = line[2]
            a1 = line[3]
            a2 = line[4]
            change_alleles=map(lambda x: a1 if(x == 0) else a2,line[5:])
            zipa=change_alleles[0::2]
            zipb=change_alleles[1::2]
            change_alleles = zip(zipa,zipb)
            change_alleles = [''.join(row) for row in change_alleles]
            output_line = rsid + ' ' + a1 + '/' + a2 + ' '+ options.chromosome + ' ' + pos
            output_line = output_line + ' + -9 -9 -9 -9 -9 -9 ' + ' '.join(change_alleles)
            output.write(output_line + '\n')
    output.close() 

if __name__=="__main__":main()
