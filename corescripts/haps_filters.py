# Implement's haps filters for HWE MISSINGNESS
# and MAF

import argparse
import scipy.stats.fisher_exact

def filter_haps_file(args):
    with open(args.haps,'r') as input_haps:
        with open(args.output,'r') as output_haps:
            for line in input_haps:
                line = line.split()[5:]
                total = len(line.split)
                if(sum(line)/total < args.maf):
                    continue
                if((1- (line.count('?')/total)) <  args.missing):
                    continue
                zipa = line[0::2]
                zipb = line[0::1]
                #count the occurence of AA AB BB alleles
                countAA = 0
                countAB = 0
                countBB = 0
                for a ,b in zip(zipa , zipb):
                    if ( a == '0' and b == '0'):
                        countAA += 1
                    elif(a == '1' and b == '1'):
                        countBB += 1
                    elif(a == '1' and b == '0'):
                        countAA += 1
                    elif(a == '0' and b == '1'):
                        countAB += 1
                countAA = float(countAA)
                countAB = float(countAB)
                countBB = float(countBB)
                p = (2 * countAA + countAB) / (2*(countAA + countAB + countBB))
                q = 1 - p


def main():
    parser = argparse.ArgumentParser(description='Preform filtering on a haps file')
    parser.add_argument('--hwe',dest='hwe')
    parser.add_argument('--haps',dest='haps')
    parser.add_argument('--output',dest='output')
    parser.add_argument('--maf',dest='hwe')
    parser.add_argument('--missing',dest='missing')    
    args = parser.parse_args()
    if(args.hwe is None):
        args.hwe = 0.0
    else:
        args.hwe = float(args.hwe)
    if(args.maf is None):
        args.maf = 0.0
    else:
        args.maf = float(args.maf)
    if(args.missing is None):
        args.missing = 1
    else 
        args.missing = float(args.missing)
    assert args.haps is not None, \
        "Haps file is required to run haps filters"    assert args.output is not None, \
        "Output file name is required to haps filter" 
    filter_haps_file(args) 
    

if __name__ == "__main__": 
    main()
