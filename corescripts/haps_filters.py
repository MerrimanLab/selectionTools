# Implement's haps filters for HWE MISSINGNESS
# and MAF

import argparse



def main():
    parser = argparse.ArgumentParser(description='Preform filtering on a haps file')
    parser.add_argument('--hwe',dest='hwe')
    parser.add_argument('--haps',dest='haps')
    parser.add_argument('--output',dest='output')
    parser.add_argument('--maf',dest='hwe')
    parser.add_argument('--missing',dest='missing')    args = parser.parse_args()
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
    


if __name__ == "__main__": 
    main()
