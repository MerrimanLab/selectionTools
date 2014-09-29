# Implement's haps filters for HWE MISSINGNESS
# and MAF

import argparse
try:
    from  scipy import stats
except ImportError as inst:
    print ("Could not import scipy install to use haps filters")
    raise inst


def hardy_weinberg_asymptotic(obs_het, obs_a , obs_b):
    obs_het = float(obs_het)
    obs_a = float(obs_a)
    obs_b = float(obs_b)
    sample_size = obs_het + obs_a + obs_b
    p = (((2 * obs_a) + obs_het) / ( 2 * (sample_size)))
    q = 1 - p 
    exp_a = p * p * sample_size 
    exp_b = q * q * sample_size
    exp_ab = 2 * p * q * sample_size
    
    # get chiSquare values
    if(exp_a == 0):
        chi_a = 0
    else:
        chi_a = ((obs_a - exp_a) * 2.0) / exp_a
    if(exp_b == 0):
        chi_b = 0
    else:
        chi_b = ((obs_b - exp_b) * 2.0) / exp_b
    if(exp_ab == 0):
        chi_ab = 0
    else:
        chi_ab = ((obs_het - exp_ab) * 2.0 ) / exp_ab
    chi_sq_total = chi_a + chi_b + chi_ab
    return stats.chisqprob(chi_sq_total, 1)    
    
def hardy_weinberg_exact(obs_het, obs_a, obs_b):
   raise Exception("Not Implemented")

def filter_haps_file(args):
    with open(args.haps,'r') as input_haps:
        with open(args.output,'w') as output_haps:
            for snp_data in input_haps:
                line = snp_data.split()[5:]
                total = float(len(line))
                if(((line.count('?')/total)) > args.missing):
                    continue
                question_marks = line.count('?')
                # TriAllellic Message
                if( (question_marks + line.count('0') + line.count('1'))!= int(total)):
                    continue
                p = line.count('0')
                q = line.count('1')
                major = p
                minor = q
                if(q > major):
                    major = q
                    minor = p
                total = float(p + q)
                if(minor/total < args.maf):
                    continue
                zipa = line[0::2]
                zipb = line[1::2]
                countAA = 0
                countAB = 0
                countBB = 0
                for a ,b in zip(zipa , zipb):
                    if ( a == '0' and b == '0'):
                        countAA += 1
                    elif(a == '1' and b == '1'):
                        countBB += 1
                    elif(a == '1' and b == '0'):
                        countAB += 1
                    elif(a == '0' and b == '1'):
                        countAB += 1
                countAA = float(countAA)
                countAB = float(countAB)
                countBB = float(countBB)
                if(args.chi_square):
                    hwe_pvalue = \
                         hardy_weinberg_asymptotic(countAB, countAA, countBB)
                else:
                    hwe_pvalue = \
                        hardy_weinberg_exact(countAB,countAA,countBB)
                if(hwe_pvalue <= args.hwe):
                    continue
                output_haps.write(snp_data)
                
                
    
                


def main():
    parser = argparse.ArgumentParser(description='Preform filtering on a haps file')
    parser.add_argument('--hwe',dest='hwe')
    parser.add_argument('--haps',dest='haps')
    parser.add_argument('--output',dest='output')
    parser.add_argument('--maf',dest='maf')
    parser.add_argument('--missing',dest='missing')    
    parser.add_argument('--chi-sq',action='store_true',
                        dest='chi_square',default=False,
                        help="Use a chi-square test instead of an exact test")
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
    else: 
        args.missing = float(args.missing)
    assert args.haps is not None, \
        "Haps file is required to run haps filters"    
    assert args.output is not None, \
        "Output file name is required to haps filter" 
    filter_haps_file(args) 
    

