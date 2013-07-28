import sys,re

#Import vcf Library

import vcf

from optparse import OptionParser
from pyfasta import Fasta

#
# Command Line Arguments
#
# --haps phased haps
# --aa ancestral allele annotation
# --chr Chromosome
# --output Output file
# --format |High|Low|
#

#
# 10000 genomesfasta Ancestral allele file can be downloaded from
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.tar.bz2
# 


#
# Takes a FASTA File Containing the Ancestral Alleles.
# Using the same format as the 1000 genomes ancestral alleles file.
# ATCG high confidence call
# actg low confidence call
# N : failure
# - : the exant species contains a insertion at this position.
# . : no coverage at this alignment.
#  
#

# annotate haps file.

#Annotate a phased_vcf file that is a file from 1000 genomes with the ancestral alleles 
# Could potentially take a population list of ids also but for now just creates it from the 
# vcf file it is given

# TODO Shrink the classes down so that there is far less repeated code.
#


def annotate_vcf(options):
    f = Fasta(options.ancestralfasta)
    keyz = (f.keys())
    chroms = {}
    if(options.output != None):
        output = open(options.output, 'w')
    for key in keyz:
        chroms[(key.split(':'))[2]] = key
    
    aaSeq = f[chroms[options.chromosome]]
    vcf_reader = vcf.Reader(open(options.vcf_file), 'r')
    for record in vcf_reader:
        if(record.ID == None):
            record.ID = str(record.POS)
        line = record.ID + ' ' + record.ID + ' ' + str(record.POS) + ' ' + str(record.REF) + ' ' + str(record.ALT[0])
        for samples in record.samples:
            gt = samples['GT']
            # Need to skip any snps that have any missing phase data to increase certainty of our results.
            # If every snp will indeed be phased
            if( '|'  in gt):
                gtSplit = gt.split('|')
                line = line  + ' ' + gtSplit[0] + ' ' + gtSplit[1]
            else:
                line = None
                break
        if(line != None):
            output_line = aa_check(aaSeq[record.POS],record.REF,record.ALT,options.format,line)
            if(output_line != None):
                if (options.output != None):
                        output.write(output_line)
                else:
                        print(output_line)



def aa_check(realAA,ref,alt,format,line):
    if(re.match('[ACTGactg]',realAA)):
        if(realAA.islower() and format == "high"):
            return None
        else:
            if(realAA == ref):
                return line.strip() 
            elif(realAA == alt):
                newLine = line.split()
#                print newLine
                newLine[3] = alt
                newLine[4] = ref
#                print(newLine)
                for i in range(5,len(newLine)):
                    if((newLine[i]) == "1"):
                        newLine[i] = '0'
                    else:
                        newLine[i] = '1'
            else:
                newLine = line.split()
                newLine[3] = realAA
                newLine[4] = ref
                #print(newLine)
                for i in range(5,len(newLine)):
                        newLine[i] = '1'
            return ' '.join(newLine)
            #return line
    else:
        return None


def annotate_haps(options):
    f = Fasta(options.ancestralfasta)
    keyz = (f.keys())
    chroms = {}
    for key in keyz:
        chroms[(key.split(':'))[2]] = key
    
    aaSeq = f[chroms[options.chromosome]]
    output = None
    if(options.output != None):
        output = open(options.output, 'w')
    with open(options.haps, 'r') as haps:
        for line in haps:
            lineSplit = line.split()
            pos = int(lineSplit[2])
            ref = lineSplit[3]
            alt = lineSplit[4]
            tempSeq=aaSeq[pos]
            #print(tempSeq)
            outputLine = aa_check(tempSeq,ref,alt,options.format,line) 
            if(outputLine != None):
                if(options.output != None):
                    output.write(outputLine + "\n")
                else:
                    print(outputLine)


def main():
    parser = OptionParser()
    parser.add_option('-i','--haps',dest='haps',help="Haplotype File (.haps)")
    parser.add_option('-a','--aa', dest='ancestralfasta',help="Ancestral Allele FastA file")
    parser.add_option('-c','--chr',dest="chromosome",help="Chromosome")
    parser.add_option('-o','--output',dest="output",help="Output File (optional)")
    parser.add_option('-f','--format',dest="format",help="Format us High or use Low & High")
    parser.add_option('-v','--phased-vcf',dest="vcf_file",help="Phased VCF file (.vcf)")
    parser.add_option('-V','--phased-vcf-gz',dest="vcf_gz", help="Gzipped VCF not implemented")
    (options,args) = parser.parse_args()
    #print(options)
    
    # Will annotate the haps file with exactly what is required
    # More options could be added later covering a wider range of file types 
    # andy maybe different input ancestral alleles.
    if(options.haps != None):
        annotate_haps(options)
    elif(options.vcf_file != None):
        annotate_vcf(options)


    
    

if __name__=="__main__":main()
