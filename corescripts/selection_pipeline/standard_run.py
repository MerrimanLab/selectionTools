
import os
import sys
import re
from run_pipeline import CommandTemplate


#For matching the file names
import fnmatch

from optparse import OptionParser
import ConfigParser

## Subprocess import clause required for running commands on the shell##
import subprocess
import logging
logger = logging.getLogger(__name__)
SUBPROCESS_FAILED_EXIT=10

class StandardRun(CommandTemplate):
    
    def __init__(self,options,config):
        if(options.phased_vcf): 
            haps = self.ancestral_annotation_vcf(options,config)
            ihh = self.run_multi_coreihh(options,config,haps)
        else:
            (ped,map) = self.run_vcf_to_plink(options,config)
            (ped,map) = self.run_plink_filter(options,config,ped,map)
            (haps) = self.run_shape_it(options,config,ped,map) 
        if(options.imputation):
            (haps)= self.run_impute2(options,config,haps)
        haps = self.indel_filter(options,config,haps)
        #tajimas = run_tajimas_d(options,config,haps)
        haps = self.run_aa_annotate_haps(options,config,haps)
        ihh = self.run_multi_coreihh(options,config,haps)
        logger.info("Pipeline completed successfully")
        logger.info(haps)
        logger.info(ihh)
        logger.info("Goodbye :)")
 
    def run_subprocess(self,command,tool):   
        try:
            subprocess.call(command) 
        except:
            logger.error(tool + " failed to run " + ' '.join(command))
            sys.exit(SUBPROCESS_FAILED_EXIT)    
        
            
        logger.error("Finished tool " + tool)

    def run_vcf_to_plink(self,options,config):
        (cmd,prefix) = CommandTemplate.run_vcf_to_plink(self,options,config)
        self.run_subprocess(cmd,'vcftools') 
        return(prefix + '.ped', prefix + '.map')

    # Calls a subprocess to run plink

    def run_plink_filter(self,options,config,ped,map):
        (cmd,prefix) = CommandTemplate.run_plink_filter(self,options,config,ped,map)
        self.run_subprocess(cmd,'plink')
        return(prefix+'.ped',prefix+'.map')

        # Calls a subprocess to run shape it

    def run_shape_it(self,options,config,ped,map):
        (cmd,prefix) = CommandTemplate.run_shape_it(self,options,config,ped,map)
        cmd.extend(['--thread',config['system']['threads_avaliable']])
        self.run_subprocess(cmd,'shapeit')
        return(prefix + '.haps')

    #Calls a subprocess to run impute   
 
    def run_impute2(self,options,config,haps):
        cmd = []
        prefix = options.output_prefix + options.chromosome + '_impute2'
        logger.debug('Attempting to call impute2 the data')
        impute2=config['impute2']['impute_executable']
        cmd.append(impute2)
    #    cmd.append(['--
    
    #Calls a subprocess to run the indel filter

    def indel_filter(self,options,config,haps):
        (cmd,output_name) = CommandTemplate.indel_filter(self,options,config,haps)
        self.run_subprocess(cmd,'indel_filter')
        return(output_name)
 
    def run_aa_annotate_haps(self,options,config,haps):
        (cmd,output_name) = CommandTemplate.indel_filter(self,options,config,haps)
        self.run_subprocess(cmd,'ancestral_annotation')
        return (output_name)
    
#def run_tajimas_d(options,config,gen,sample):

    def run_multi_coreihh(self,options,config,haps):
        (cmd,output_name) = CommandTemplate.run_multi_coreihh(self,options,config,haps)
        cores=config['system']['threads_avaliable']
        cmd.extend(['--cores',cores])
        cmd.extend(['--working_dir','.'])
        with open(haps,'r') as hap_file:
            line = hap_file.readline()
            
   
        self.run_subprocess(cmd,'multcore_ihh')
        os.rename(options.population+'_chr_'+options.chromosome+"_wd_"+'.'+"_.ihh",output_name)
        return output_name
