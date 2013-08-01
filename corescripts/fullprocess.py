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
import re

#For matching the file names
import fnmatch

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
    parser.add_option('--population',dest="population",help="Population Code ")
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
    assert options.population is not None, "Population code has not been specified."
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
    logger.debug(config['vcftools']['extra_args'])
    cmd.extend(config['vcftools']['extra_args'].split())
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
    cmd.extend(config['plink']['extra_args'].split())
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
    genetic_map = ''
    for file in os.listdir(config['shapeit']['genetic_map_dir']):
        if fnmatch.fnmatch(file,'genetic_map_chr'+options.chromosome+'*'):
            genetic_map = file
            
    shapeit=config['shapeit']['shapeit_executable']
    cmd.append(shapeit)
    cmd.extend(['--input-ped',ped,map,'-M',os.path.join(config['shapeit']['genetic_map_dir'],genetic_map),'--output-max',prefix,'--thread',config['system']['threads_avaliable']])
    cmd.extend(config['shapeit']['extra_args'].split())
    logger.debug(cmd)
    try:
        subprocess.call(cmd)
    except:
        logger.error("shapeit failed to run" + ' '.join(cmd))
        sys.exit(SUBPROCESS_FAILED_EXIT)
    logger.debug('Shape it phasing has completed')
    return(prefix + '.haps')

#Calls a subprocess to run impute   
 
def run_impute2(options,config,haps):
    cmd = []
    prefix = options.output_prefix + options.chromosome + '_impute2'
    logger.debug('Attempting to call impute2 the data')
    impute2=config['impute2']['impute_executable']
    cmd.append(impute2)
#    cmd.append(['--
    
#Calls a subprocess to run the indel filter

def indel_filter(options,config,haps):
    cmd = []    
    output_name= options.output_prefix + options.chromosome + '_indel_filter.haps'        
    logger.debug('Attempting to run the R indel and maf filter usually reserved for after phasing')
    rscript = config['Rscript']['rscript_executable']
    indel_filter = config['Rscript']['indel_filter']
    cmd.append(rscript)
    cmd.append(indel_filter)
    cmd.extend([haps,str(options.maf),output_name])
    try:
        subprocess.call(cmd)
    except:
        logger.error("maf filter failed to run" + ' '.join(cmd))
        sys.exit(SUBPROCESS_FAILED_EXIT)
    
    return(output_name)
 
def run_aa_annotate_haps(options,config,haps):
    cmd = []
    output_name= options.output_prefix + '_aachanged.haps'
    py_executable = config['python']['python_executable']
    aa_annotate = config['ancestral_allele']['ancestral_allele_script']
    logger.debug('Attempting to run ancestral allele annotation')
    for file in os.listdir(config['ancestral_allele']['ancestral_fasta_dir']):
        if fnmatch.fnmatch(file,config['ancestral_allele']['ancestral_prefix'].replace('chr',options.chromosome)):
            ancestral_fasta = file
    cmd.append(py_executable)
    cmd.append(aa_annotate)
    cmd.extend(['-i',haps ,'-c', options.chromosome, '-o', output_name,'-a',os.path.join(config['ancestral_allele']['ancestral_fasta_dir'],ancestral_fasta)])
    try:
        subprocess.call(cmd)
    except:
        logger.error("ancestral allele annotation failed to run" + ' '.join(cmd))
        sys.exit(SUBPROCESS_FAILED_EXIT)
    return output_name
    
#def run_tajimas_d(options,config,gen,sample):

def run_multi_coreihh(options,config,haps):
    cmd=[]
    output_name= options.output_prefix + '.ihh'
    rscript=config['Rscript']['rscript_executable']
    multicore_ihh=config['multicore_ihh']['multicore_ihh']
    window=config['multicore_ihh']['window']
    overlap=config['multicore_ihh']['overlap']
    cores=config['system']['threads_avaliable']
    logger.debug("Started running multicore iHH (rehh) script")
    # Default offset is 0 as this is the single pc operation something different happens on nesi
    population=options.population
    cmd.append(rscript)
    # Todo look at MAF in rehh
    cmd.extend([multicore_ihh,population,haps,str(options.chromosome),str(window),str(overlap),str(cores),'.','7',str(config['multicore_ihh']['minor_allele_frequency']),options.output_prefix])
    try:
        subprocess.call(cmd)
    except:
        logger.error("multicore ihh (rehh) failed to run" + ' '.join(cmd))
        sys.exit(SUBPROCESS_FAILED_EXIT)
    os.rename(population+'_chr_'+options.chromosome+"_wd_"+'.'+"_.ihh",output_name)
    return output_name

def parse_config(options):
    config = ConfigParser.ConfigParser()
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

def level_the_load(options,arguments):
    #RUN _ON_ NESI
    return 1


def main():
    options = parse_arguments()
    config = parse_config(options)
    if(config['system']['nesi']==True):
        level_the_load(options,config)
    elif(options.phased_vcf):
        haps = ancestral_annotation_vcf(options,config)
        ihh = run_multi_coreihh(options,config,haps)
    else:
        (ped,map) = run_vcf_tools(options,config)
        (ped,map) = run_plink(options,config,ped,map)
        (haps) = run_shape_it(options,config,ped,map) 
        if(options.imputation):
            (haps)= run_impute2(options,config,haps)
        haps = indel_filter(options,config,haps)
        #tajimas = run_tajimas_d(options,config,haps)
        haps = run_aa_annotate_haps(options,config,haps)
        ihh = run_multi_coreihh(options,config,haps)
    logger.info("Pipeline completed successfully")
    logger.info(haps)
    logger.info(ihh)
    logger.info("Goodbye :)")
                
if __name__=="__main__":main()

