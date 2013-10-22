#
#
# Multipopulation script calls the selectio
# pipeline for each population that we need
# to do then zips up and runs a script to p# each of the cross population statistics once
# it has all finished.
# institution: University of Otago
# author: James Boocock
#
#
# requires that the selection pipeline is 
# installed.
#

import sys
import os
from optparse import OptionParser

def get_populations(populations):
    pops = {}
    for pop in populations:
        with open(pop, 'r') as f:
            for i, line in enumerate(f):
                if ( i == 0 ):
                    pops[pop]=[]
                    pop_name = pop
                else:
                    pops[pop_name].append(pop)
    return pops

def main():
    parser=OptionParser()
    parser.add_option('-p','--population',action='append',dest="populations",help='population_files')
    (options,args) = parser.parse_args()
    pop_names=get_populations(options.populations)
     




if __name__=="__main__":main()
