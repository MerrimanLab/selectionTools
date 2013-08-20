import queue
import os
import sys
import re
from run_pipeline import CommandTemplate


#For matching the file names
import fnmatch

from optparse import OptionParser

## Subprocess import clause required for running commands on the shell##
import subprocess
import logging
logger = logging.getLogger(__name__)
SUBPROCESS_FAILED_EXIT=10
MISSING_EXECUTABLE_ERROR=5

class StandardRun(CommandTemplate):
   
    def is_script(self,fpath):
        return os.path.isfile(fpath)
    def is_exe(self,fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    #Stolen code from 
    #http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    def which(self,program,program_name):
        fpath, fname = os.path.split(program)
        print(program)
        if fpath:
            if self.is_exe(program):
                return program
            elif (self.is_script(program)):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if self.is_exe(exe_file):
                    return exe_file
        logger.error(program_name +" path = " + fpath+" not locatable path or in the directory specified in your config file ")
        return None    
  
    #Ensures the executables specified on the path exist so that the standard run will
    # work on the machine 
    def check_executables_and_scripts_exist(self,options,config):
        executables=['plink','shapeit','impute','Rscript','python','ancestral_allele','indel_filter','multicore_ihh']
        if(self.which(config['plink']['plink_executable'],'plink')is None): 
            return False
        if(self.which(config['shapeit']['shapeit_executable'],'shapeit') is None ):
            return False
        if(self.which(config['ancestral_allele']['ancestral_allele_script'],'ancestral_allele') is None):
            return False
        if(self.which(config['impute']['impute_executable'],'impute2') is None):
            return False
        if(self.which(config['Rscript']['indel_filter'],'indel_filter') is None):
            return False
        if(self.which(config['Rscript']['rscript_executable'],'Rscript') is None):
            return False
        if(self.which(config['multicore_ihh']['multicore_ihh'],'multicore_ihh') is None):
            return False
        return True
        


    def __init__(self,options,config):
        #Perform local executable check.# 
        if( not self.check_executables_and_scripts_exist(options,config)):
            sys.exit(MISSING_EXECUTABLE_ERROR)
        self.threads=config['system']['threads_avaliable']        
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
            exit_code = subprocess.call(command) 
        except:
            logger.error(tool + " failed to run " + ' '.join(command))
            sys.exit(SUBPROCESS_FAILED_EXIT)   

        if(exit_code != 0): 
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
        cmd.extend(['--thread',self.threads])
        self.run_subprocess(cmd,'shapeit')
        return(prefix + '.haps')

    #Calls a subprocess to run impute   

        
    def impute_worker(q):
        while True:
            cmd=q.get()
            self.run_subprocess(cmd,'impute2')
            q.task_done()             
 
    def run_impute2(self,options,config,haps):
        imputeQueue=queue.Queue()
        (cmds,output_prefix) = CommandTemplate.run_impute2(self,options,config,haps)
        threads=self.threads
        no_commands=len(cmds)
        for i in range(thread):
            t = Thread(target=impute_worker,args=imputeQueue)
            t.daemon = True
        for cmd in cmds:
            q.put(cmd)
        impute_queue.join()
        CommandTemplate.join_impute2_files(options,config,output_prefix,no_commands)
        return(output_prefix+'.haps') 
         

    def indel_filter(self,options,config,haps):
        (cmd,output_name) = CommandTemplate.indel_filter(self,options,config,haps)
        self.run_subprocess(cmd,'indel_filter')
        return(output_name)
 
    def run_aa_annotate_haps(self,options,config,haps):
        (cmd,output_name) = CommandTemplate.indel_filter(self,options,config,haps)
        self.run_subprocess(cmd,'ancestral_annotation')
        return (output_name)
    
    #def run_tajimas_d(options,config,gen,sample):
    
    #def fu_and_wus_h(options,config,gen,sample):

    def run_multi_coreihh(self,options,config,haps):
        (cmd,output_name) = CommandTemplate.run_multi_coreihh(self,options,config,haps)
        cores=self.threads
        cmd.extend(['--cores',cores])
        cmd.extend(['--working_dir','.'])
        cmd.extend(['--offset',0])
        with open(haps,'r') as hap_file:
            line = hap_file.readline()
            
   
        self.run_subprocess(cmd,'multcore_ihh')
        os.rename(options.population+'_chr_'+options.chromosome+"_wd_"+'.'+"_.ihh",output_name)
        return output_name
