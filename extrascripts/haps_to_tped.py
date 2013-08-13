#!/usr/bin/env python
#
# Script for converting from haps to tped
# 
# haps to tped.
#
import os
import logging
import shutil
from optparse import OptionParser

logging.basicConfig(format="%(asctime)s %(message)s")

log = logging.getLogger(__name__)
log.setLevel(level=logging.INFO)

def make_allelic_lambda(ref,alt):
    return lambda x: ref if 0 == int(x) else alt 

def hap_to_tped(basename,output_basename,chromosome):
    log.debug("Started converting from HAPS to TPED")
    tped_out = open(output_basename + '.tped','w')
    # Use shutil to copy sample to tfam
    #default centimorgan values.
    i = 0
    tfam_out = open(output_basename + '.tfam','w')
    with open(basename + '.sample','r') as f:
        for line in f:
            if i > 1:
                newline=line.split()
                newline=line[0:5]
                tfam_out.write(' '.join(newline) +'\n')
            i = i + 1
    tfam_out.close()
    centimorgans = '0'
    with open(basename + '.haps','r') as f:
        log.debug("Read in TPED file")
        for line in f:
            split_line = line.split()
            rsid = split_line[0]
            pos = split_line[2]
            ref = split_line[3]
            alt = split_line[4]
            allele_lambda = make_allelic_lambda(ref,alt)
            alleles = map(allele_lambda,split_line[5:])
            tped_out.write(chromosome + ' ' + rsid + ' ' + centimorgans + ' ' + pos + ' ' +' '.join(alleles) + '\n') 
    tped_out.close()
    log.debug("Finished converting from HAPS to TPED")

def main():
    parser= OptionParser()
    parser.add_option("-i",'--file',dest="haps_files",help="Base names of haps / sample pair")
    parser.add_option('-v','--verbose',action="store_true",dest="verbose",help="Verbosity of the logging outputs")
    parser.add_option('-o','--output',dest="output",help="Output file basename (w+) used so better clean up directory after running program")
    parser.add_option('-c','--chromosome',dest="chromosome",help="Chromosome haps file originates from")
    (options,args) = parser.parse_args()
    if(options.verbose == None):
        log.setLevel(logging.ERROR)
    else:
        log.setLevel(logging.DEBUG)
    assert (options.haps_files) is not None, "Atleast two haps file base names are required"
    assert options.output is not None, "Output base name needs to be specified"
    hap_to_tped(options.haps_files,options.output,options.chromosome)

if __name__=="__main__":main() 
