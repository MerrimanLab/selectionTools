#!/usr/bin/env python
#
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

from optparse import OptionParser
import ConfigParser

## Subprocess import clause required for running commands on the shell##
import subprocess
import logging
logging.basicConfig(format='%(asctime)s %(message)s')
logger = logging.getLogger(__name__)
logger.setLevel(level=logging.INFO)

SUBPROCESS_FAILED_EXIT=10

def parse_arguments():
    parser = OptionParser()
    parser.add_option('-v','--verbose',action="store_true",dest='verbose',help="Print debug messages")
    parser.add_option('-q','--silent',action="store_false",dest='verbose', help="Run Silently")
    parser.add_option('-i','--vcf',dest='vcf_input',help="VCF input file")
    parser.add_option('-o','--out',dest='output_prefix', help="Output file prefix")
    parser.add_option('-c','--chromosome',dest='chromosome',help="Chromosome")
    parser.add_option('-l','--log-fire',dest='log_file',help="Log file for the pipeline process")
    parser.add_option('--maf',dest='maf',help='Minor allele-frequency filter') 
    parser.add_option('--hwe',dest='hwe',help="Hardy-Weinberg Equillibrium filter thread")
    parser.add_option('--remove-missing',dest="remove_missing",help="Remove missing genotypes") 
    parser.add_option('--config-file',dest="config_file", help="Config file")
    parser.add_option('--phased-vcf',action="store_true",dest="phased_vcf",help="Phased vcf file")
    parser.add_option('--imputation',action="store_true", dest="imputation",help="Imputation")
    parser.add_option('--full-process',action="store_true",dest="full_process",help="Run Entire Process")
      
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
    #Optional arguments using sane defaults  
     
    if(options.imputation is None):
        options.imputation = False
    if(options.hwe is None):
        options.hwe = 0.001
    if(options.maf is None):
        options.maf = 0.05
    if(options.remove_missing is None):
        options.remove_missing = 0.99
    if (options.config_file == None):
        options.config_file = "defaults.cfg"
    if (options.phased_vcf == None):
        options.phased_vcf = False
    if (options.full_process == None):
        options.full_process= False
    logger.debug(options.config_file)
    return options 
     
# Calls a subprocess to run vcf tools

def run_vcf_tools(options,config):
    cmd = []
    prefix = options.output_prefix + options.chromosome
    logger.debug("Attempting to call vcf tools to convert to ped/map plink format")
    vcf_tools=config['vcftools']['vcf_tools_executable']
    cmd.append(vcf_tools)
    cmd.extend(['--gzvcf',options.vcf_input, '--plink', '--out',prefix,'--remove-indels'])
    try:
        subprocess.call(cmd) 
    except:
        logger.error("Vcf tools failed to run" + ' '.join(cmd))
        sys.exit(SUBPROCESS_FAILED_EXIT)
     
    logger.debug("Finished vcf tools run")
    return(prefix + '.ped', prefix + '.map')

# Calls a subprocess to run plink

def run_plink(options,config,ped,map):
    logger.debug("Attempting to call vcf tools to convert to ped/map plink format")
    cmd = []
    prefix = ped.split('.')[0]
    plink = config['plink']['plink_executable']
    cmd.append(plink)
    # add standard plink commands #

    cmd.extend(['--noweb','--file',prefix,'--geno',str(options.remove_missing),'--hwe',str(options.hwe),'--maf',str(options.maf),'--recode','--out',prefix])
    try:
        subprocess.call(cmd)
    except:
        logger.error("plink failed to run" + ' '.join(cmd))
        sys.exit(SUBPROCESS_FAILED_EXIT)
    logger.error("Finished plink filtering")
    return(prefix+'.ped',prefix+'.map')

# Calls a subprocess to run shape it

def run_shape_it(options,config,ped,map):
    cmd = []
    prefix = options.output_prefix + options.chromosome + '.phased'
    logger.debug("Attempting to call shape it to phase the data")
    shapeit=config['shapeit']['shapeit_executable']
    cmd.append(shapeit)
    cmd.extend(['--input-ped',ped,map,'-M',config['shapeit']['genetic_map_dir'],'--output-max',prefix,'--thread',config['system']['threads_avaliable']])
    try:
        subprocess.call(cmd)
    except:
        logger.error("shapeit  failed to run" + ' '.join(cmd))
        sys.exit(SUBPROCESS_FAILED_EXIT)
    return(prefix + '.haps')

#Calls a subprocess to run impute   
 
def run_impute2(options,config,gen,sample):
    cmd = []
    prefix = options.output_prefix + options.chromosome + '_impute2'
    logger.debug('Attempting to call impute2 the data')
    impute2=config['impute2']['impute_executable']
    cmd.append(impute2)
#    cmd.append(['--
    
#Calls a subprocess to run the indel filter

def indel_filter(options,config,gen,sample):
    cmd = []    
    output_name= options.output_name + options.chromosome + '_indel_filter.haps'        
    r_executable = config['Rscript']['r_executable']
    indel_filter = config['Rscript']['indel_filter']
    cmd.append(r_executable)
    cmd.append([])
 
def run_aa_annotate(options,config,gen,sample):
    cmd = []
    py_executable = config['python']['python_executable']
    aa_annotate = config['python']['aa_annotate']
    cmd.append(py_executable)
    
#def run_tajimas_d(options,config,gen,sample):

#def run_multi_coreihh(options,config,gen,sample):
    

def parse_config(options):
    config = ConfigParser.ConfigParser()
    config.read(options.config_file)
    config_parsed = {}
    logger.debug(config.sections())
    logger.debug(config.get('system','ram_avaliable'))
    for section in config.sections():
        logger.debug(section)
        opts = config.options(section)
        config_parsed[section] = {}
        for op in opts:
            try:
                config_parsed[section][op] = config.get(section,op)
            except:
                logger.info("exception on {0}".format(op))
                config_parsed[section][op] = None
    return config_parsed

def main():
    options = parse_arguments()
    config = parse_config(options)
    if(options.phased_vcf):
        ancestral_annotation(options,config)
    else:
        (ped,map) = run_vcf_tools(options,config)
        (ped,map) = run_plink(options,config,ped,map)
        (gen) = run_shape_it(options,config,ped,map) 
        if(options.imputation):
            (haps)= run_impute2(options,config,haps)
        haps = indel_filter(options,config,gen)
        #tajimas = run_tajimas_d(options,config,haps)
        #haps = run_aa_annotate(options,config,haps)
                
if __name__=="__main__":main()

