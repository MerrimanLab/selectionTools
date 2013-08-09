#!/usr/bin/env python
#
# Script for concatenating haps files. 
#
# @author James Boocock
# @date 09/08/2013
#

from optparse import OptionParser
import logging

logging.basicConfig(format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)
log.setLevel(level=logging.INFO)

#If there are any repeat ids between the two columns use the first one#
#
# Haps have no header line.
# file looks like.
# RSID RSID POS REF ALT HAPLOTYPES....

def merge_haps(one_haps,two_haps,output_file):
    merged_haps = open(output_file + '.haps','w')
    merged_sample = open(output_file + '.sample','w')
    with open(one_haps + '.sample') as one_f:
        with open(two_haps + '.sample') as two_f:
            for line in one_f:
                merged_sample.write(line)
            i = 0
            for line in two_f:
                if ( i >= 2):
                    merged_sample.write(line)
                i = i + 1
    log.debug("Finished creating Merged sample file for " + one_haps +" " + two_haps)
    with open(one_haps+'.haps','r') as one_f:
        with open(two_haps+'.haps','r') as two_f:
            one_first_line=one_f.readline()
            two_first_line=two_f.readline()
            while(one_first_line != '' and two_first_line != ''):
                one_split = one_first_line.strip().split()
                two_split = two_first_line.strip().split()
                pos_one = one_split[2]
                pos_two = two_split[2]
                if(int(pos_two) == int(pos_two)):
                    merged_haps.write(' '.join(one_split) + ' ' + ' '.join(two_split[5:]) + '\n')
                    one_first_line=one_f.readline()
                    two_first_line=two_f.readline()
                elif (int(pos_one) > int(pos_two)):
                    one_first_line=one_f.readline()
                else: 
                    two_first_line=two_f.readline()
    log.debug("Finished merging haps files for " + one_haps+  " " + two_haps)
    merged_haps.close()
    return output_file



def main():
    parser = OptionParser()
    parser.add_option("-i",'--file',dest="haps_files",action="append",help="Base names of haps / sample pairs atleast two files are required")
    parser.add_option('-v','--verbose',action="store_true",dest="verbose",help="Verbosity of the logging outputs")
    parser.add_option('-d','--delete-missing-snps',action="store_true",dest="remove",help="Remove the snps that do not occur in both files")
    parser.add_option('-o','--output',dest="output",help="Output file basename (w+) used so better clean up directory after running program")
    (options,args) = parser.parse_args()
    if(options.verbose == None):
        log.setLevel(logging.ERROR)
    else:
        log.setLevel(logging.DEBUG)
    assert len(options.haps_files) > 1, "Atleast two haps file base names are required"
    assert options.output is not None, "Output base name needs to be specified"
    haps_files = options.haps_files
    first=haps_files[0]
    log.debug(haps_files)
    for i in range(1,len(haps_files)):
        second=haps_files[i]
        first = merge_haps(first,second,options.output)


if __name__=="__main__":main()

