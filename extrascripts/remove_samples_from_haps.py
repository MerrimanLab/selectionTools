import os
import sys

from optparse import OptionParser

def index_ids(sample_file,id,basename):
    ids = []
    comma_ids = id.split(',')
    samples = []
    samples_line = []
    new_samples = open(basename+'.sample','w')
    with open(sample_file,'r') as f:
        i = 0
        for line in f:
            if (i >= 2): 
                samples.append(line.split()[0])
                samples_line.append(line)
            else:
                new_samples.write(line)
            i = i + 1
    for id in comma_ids:
        ids.append(samples.index(id))
    ids.sort()
    for id in ids:
        new_samples.write(samples_line[id])
    new_samples.close()
    return ids

def haps_keep_samples(haps_file,indexed_ids,basename):
    new_haps = open(basename+'.haps','w')
    with (open(haps_file,'r')) as f:
        for line in f:
            lineSplit = line.split()
            new_haps.write(' '.join((lineSplit[:5]))+' ')
            i = 0
            i_person = 0
            person = indexed_ids[i_person] 
            for item in lineSplit[5:]:
                if((person == i/2) and (i % 2 == 0)):
                    new_haps.write(item+' ')
                elif((person <= i/2)):
                    new_haps.write(item+' ')
                    i_person = i_person + 1
                    if(i_person == len(indexed_ids)):
                        break
                    person = indexed_ids[i_person]
                i = i + 1 
            new_haps.write('\n')
    new_haps.close() 
              
def main():
    parser=OptionParser()
    parser.add_option('-i','--haps',dest='haps',help="Haplotype File (.haps)")
    parser.add_option('-s','--sample',dest="sample",help="Sample File (.sample)")
    parser.add_option('-c','--keep',dest="keep" ,help="Comma seperated list of IDs to keep") 
    parser.add_option('-o','--basename',dest="basename",help="Output haps file")
    (options,args) = parser.parse_args()
    indexed_ids = index_ids(options.sample,options.keep,options.basename)
    haps_keep_samples(options.haps,indexed_ids,options.basename)
    
if __name__=="__main__":main()
