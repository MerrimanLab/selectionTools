# Python Script to perform for running the single process for our pipeline
#
# Murray Cadzow 
# July 2013
# University Of Otago
#
# James Boocock
# July 2013
# University Of Otago
# 

import os
import sys
import re

#For matching the file names
import fnmatch

from optparse import OptionParser
import configparser

## Subprocess import clause required for running commands on the shell##
import subprocess
import logging

#Import standard run 
from .load_leveler_run import LoadLevelerRun 
from .standard_run import StandardRun

logging.basicConfig(format='%(asctime)s %(message)s')
logger = logging.getLogger(__name__)
logger.setLevel(level=logging.INFO)

SUBPROCESS_FAILED_EXIT=10


def parse_config(options):
    config = configparser.ConfigParser()
    config.read(options.config_file)
    config_parsed = {}
    logger.debug(config.sections())
    for section in config.sections():
        logger.debug(section)
        opts = config.options(section)
        config_parsed[section] = {}
        for op in opts:
            logger.debug(op)
            try:
                config_parsed[section][op] = config.get(section,op)
            except:
                logger.info("exception on {0}".format(op))
                config_parsed[section][op] = None
    return config_parsed

def parse_arguments():
    parser = OptionParser()
    parser.add_option('-v','--verbose',action="store_true",dest='verbose',help="Print debug messages")
    parser.add_option('-q','--silent',action="store_false",dest='verbose', help="Run Silently")
    parser.add_option('-i','--vcf',dest='vcf_input',help="VCF input file")
    parser.add_option('-o','--out',dest='output_prefix', help="Output file prefix")
    parser.add_option('-c','--chromosome',dest='chromosome',help="Chromosome")
    parser.add_option('-l','--log-fire',dest='log_file',help="Log file for the pipeline process")
    parser.add_option('--maf',dest='maf',help='Minor allele-frequency filter') 
    parser.add_option('--hwe',dest='hwe',help="Hardy-Weinberg Equillibrium filter proportion")
    parser.add_option('--daf',dest='daf',help="Derived Allele Frequency filter proportion")
    parser.add_option('--remove-missing',dest="remove_missing",help="Remove missing genotypes") 
    parser.add_option('--config-file',dest="config_file", help="Config file")
    parser.add_option('--phased-vcf',action="store_true",dest="phased_vcf",help="Phased vcf file")
    parser.add_option('--population',dest="population",help="Population Code ")
    parser.add_option('--imputation',action="store_true", dest="imputation",help="Imputation")
    parser.add_option('--full-process',action="store_true",dest="full_process",help="Run Entire Process")
    parser.add_option('--gzvcf',action="store_true",dest="vcf_gz",help="VCF input is in GZ file (optional)")
    parser.add_option('--TajimaD',dest='tajimas_d',help="Output Tajima's D statistic in bins of size <int>")
      
    (options, args) = parser.parse_args()
    if(options.verbose != None):
        if(options.verbose):
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.ERROR)
    # Obligatory arguments
    assert options.vcf_input is not None, "No VCF file has been specified as input"
    assert options.chromosome is not None, "No chromosome has been specified to the script"
    assert options.output_prefix is not None, "Output file prefix has not been specified."
    assert options.population is not None, "Population code has not been specified."
    #Optional arguments using sane defaults  

    if(options.tajimas_d is None):
        # default tajimas D 1000 bin size
        options.tajimas_d = str(1000)
    else:
        options.tajimas_d = str(options.tajimas_d)
    if(options.imputation is None):
        options.imputation = False
    if(options.hwe is None):
        options.hwe = 0.001
    if(options.maf is None):
        options.maf = 0.01
    if(options.daf is None):
        options.daf = 0.00
    if(options.remove_missing is None):
        options.remove_missing = 0.99
    if (options.config_file == None):
        options.config_file = "defaults.cfg"
    if (options.phased_vcf == None):
        options.phased_vcf = False
    if (options.full_process == None):
        options.full_process= False
    if (options.vcf_gz == None):
        options.vcf_gz = False
    return options 
     
# Calls a subprocess to run vcf tools


def main():
    options = parse_arguments()
    config = parse_config(options)
    if(config['system']['nesi']=="True"):
        LoadLevelerRun(options,config)
    else:
        StandardRun(options,config)        
                
if __name__=="__main__":main()

