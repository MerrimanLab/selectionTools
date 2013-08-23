#
# Python Script to perform for running the full process using load leveler
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

from run_pipeline import CommandTemplate

from optparse import OptionParser

## Subprocess import clause required for running commands on the shell##
import subprocess
import logging
#For giving steps unique names
from datetime import datetime
#Import standard run 

logging.basicConfig(format='%(asctime)s %(message)s')
logger = logging.getLogger(__name__)
logger.setLevel(level=logging.INFO)

SUBPROCESS_FAILED_EXIT=10
load_leveler_template="""
        #@ shell = /bin/bash
        #@ group = {0}
        #@ class = {1}
        #@ output = $(jobid).out
        #@ error = $(jobid).err 
"""
load_leveler_step="""
        #@ step_name = {0}
        #@ wall_clock_limit = {1}
        #@ resources = ConsumableMemory({2}gb) ConsumableVirtualMemory({2}gb)
"""

load_leveler_account_no="""
        #@ account_no = {0}
"""

# serial_task_with_threads
load_leveler_serial ="""
        #@ job_type = serial
        #@ parallel_threads = {0}
"""
# mpi_task
load_leveler_mpi = """
        #@ job_type = parallel
        #@ total_tasks_threads = {0}
"""
load_leveler_dependency = """
        #@ dependency = ({0})
""" 
# Queue jobs on load leveler
load_leveler_queue ="""
        #@ queue
    """

# ulimit load leveler ulimit
load_leveler_ulimit ="""
            ulimit -v {0} -m {0}
    """



class LoadLevelerRun(object):
    
    """ Load leveler class takes the pipeline and runs the PHD on the nesi pan 
        cluster.
    """
    #We run this to generate a script for every job
    def write_job_preamble(self,memory_required,cmd):
        ulimit = str(int(memory_required) * 1024 * 1024)
        self.load_leveler_script.write(self.string_to_bytes(load_leveler_ulimit.format(ulimit)))
        self.load_leveler_script.write(self.string_to_bytes(' '.join(cmd)))
        self.load_leveler_script.write(self.string_to_bytes('\n'))
    def queue_ll(self):
        self.load_leveler_script.write(self.string_to_bytes(load_leveler_queue+'\n'))
    def mpi_task(self,threads):
        self.load_leveler_script.write(self.string_to_bytes(load_leveler_mpi.format(threads)))
    def serial_task(self,threads):
        self.load_leveler_script.write(self.string_to_bytes(load_leveler_serial.format(threads)))

    def write_step_preamble(self,memory_required,wall_time,step_name,dependencies=None):
        self.load_leveler_script.write(self.string_to_bytes(self.script_template)) 
        if(dependencies is not None):
            self.load_leveler_script.write(self.string_to_bytes(load_leveler_dependency.format(dependencies)))
        self.load_leveler_script.write(bytes(load_leveler_step.format(step_name,wall_time,memory_required),'UTF-8'))

    def get_dependency(self,dependency):
        return dependency + " >=0"
    #Could potentially fuck out a whole part filtered out of your daat#
    def __init__(self,options,config):
        self.options=options
        self.config=config
        logger.debug('Running the script on nesi')
        self.group=config['nesi']['group'] 
        self.nesi_class=config['nesi']['class']
        self.account_no=config['nesi']['account_no']
        self.load_leveler_script = open('level_the_load.ll','wb')
        #The template for every script
        self.script_template=load_leveler_template.format(self.group,self.nesi_class)
        if(self.account_no is not None):
             self.script_template = self.script_template + (load_leveler_account_no.format(self.account_no))
        self.create_load_leveler_script(options,config)
        self.load_leveler_script.close()

    def create_load_leveler_script(self,options,config):
        if(options.phased_vcf): 
            (dependecies,haps) = self.ancestral_annotation_vcf(options,config)
            (dependencies,haps) = self.run_multi_coreihh(options,config,haps,dependencies)
        else:
            (dependencies,ped,map) = self.run_vcf_to_plink(options,config)
            (dependencies,ped,map) = self.run_plink_filter(options,config,ped,map,dependencies)
            (dependencies,haps) = self.run_shape_it(options,config,ped,map,dependencies) 
        if(options.imputation):
            (dependencies,haps)= self.run_impute2(options,config,haps,dependencies)
        (dependencies,haps) = self.indel_filter(options,config,haps,dependencies)
        #tajimas = run_tajimas_d(options,config,haps)
        (dependencies,haps) = self.run_aa_annotate_haps(options,config,haps,dependencies)
        (dependencies,ihh) = self.run_multi_coreihh(options,config,haps,dependencies)
        logger.info("Pipeline completed successfully")
        logger.info(haps)
        logger.info(ihh)
        logger.info("Goodbye :)")
    # Need to make multicore ihh work using only 
    def run_multi_coreihh(self,options,config,haps,dependencies):
        (cmd,prefix) = CommandTemplate.run_multi_coreihh(self,options,config,haps)
        threads = '10'
        cmd.extend(['--cores',threads])
        cmd.extend(['--working_dir','.'])
        cmd.extend(['--offset','1'])
        memory_required=str(int(threads)*3)
        wall_time="11:59:00"
        step_name = prefix + self.get_date_string()
        self.write_step_preamble(memory_required,wall_time,step_name)
        #10 threads for shapeit
        self.serial_task(10)
        self.queue_ll()
        self.write_step_preamble(memory_required,wall_time,dependencies) 
        self.write_job_preamble(memory_required,cmd)
        return( step_name + ">=0 ",output_name) 

    def run_aa_annotate_haps(self,options,config,haps,dependencies):
        (cmd,prefix) =CommandTemplate.run_aa_annotate_haps(self,options,config,haps)
        memory_required="4"
        wall_time="11:59:00"
        step_name = prefix + self.get_date_string()
        self.write_step_preamble(memory_required,wall_time,step_name)
        #10 threads for shapeit
        self.serial_task(10)
        self.queue_ll()
        self.write_step_preamble(memory_required,wall_time,dependencies) 
        self.write_job_preamble(memory_required,cmd)
        return( step_name + ">=0 ",output_name) 

    def indel_filter(self,options,config,haps,dependencies):
        (cmd,prefix) = CommandTemplate.indel_filter(self,options,config,haps)
        memory_required="4"
        wall_time="11:59:00"
        step_name = prefix + self.get_date_string()
        self.write_step_preamble(memory_required,wall_time,step_name)
        #10 threads for shapeit
        self.serial_task(10)
        self.queue_ll()
        self.write_step_preamble(memory_required,wall_time,dependencies) 
        self.write_job_preamble(memory_required,cmd)
        return(self.get_dependency(step_name),prefix) 

    def run_shape_it(self,options,config,ped,map,dependencies):
        logger.debug("Preparing shapeit for running on nesi")
        (cmd,prefix) = CommandTemplate.run_shape_it(self,options,config,ped,map)
        cmd.extend(['--thread','10'])
        memory_required="24"
        wall_time="11:59:00"
        step_name = prefix + self.get_date_string()
        self.write_step_preamble(memory_required,wall_time,step_name)
        #10 threads for shapeit
        self.serial_task(10)
        self.queue_ll()
        self.write_job_preamble(memory_required,cmd)
        logger.debug("Finished preparing shape it for running on nesi")
        return(step_name +">=0",prefix + '.haps')
    def run_vcf_to_plink(self,options,config):
        logger.debug("Preparing vcf_to_plink for running on pan")
        (cmd,prefix) = CommandTemplate.run_vcf_to_plink(self,options,config)
        memory_required="4"
        wall_time="11:59:00"
        step_name = prefix + self.get_date_string()
        self.write_step_preamble(memory_required,wall_time,step_name)
        self.serial_task(1)
        self.queue_ll()
        self.write_job_preamble(memory_required,cmd)
        #and queue
        logger.debug("Finished preparing vcf_to_plink for running on pan")
        return(step_name + '>= 0',prefix + '.ped', prefix+'.map',) 
    def join_impute2_files(self,output_prefix,no_commands,dependencies):
        step_name = prefix + self.get_date_string()
        haps=[]
        warnings=[]
        info=[]
        for i in range(no_commands):
            haps.append(output_prefix+'_'+str(i)+'.haps')
            warnings.append(output_prefix+'_'+str(i)+'.warnings')
            info.append(output_prefix+'_'+str(i)+'.info')
        command = 'cat {0} > {1}'.format(' '.join(haps)) + '\n'
        command = command +  'cat {0} > {1}'.format(' '.join(haps)) + '\n'
        command = command + 'cat {0} > {1}'.format(' '.join(haps)) + '\n'
        memory_required="4"
        wall_time="11:59:00"
        step_name = prefix + self.get_date_string()
        self.write_step_preamble(memory_required,wall_time,step_name," && ".join(dependencies))
        self.serial_task(1)
        self.queue_ll()
        self.write_job_preamble(memory_required,command)
        return(self.get_dependency(step_name))
    def create_impute_commands(self,prefix,haps,cmd_template,no_jobs):
        cmds=[]
        for i in range(0,no_of_impute_jobs):
            individual_command=cmd_template
            individual_command.extend(['-int',str(i*no_of_impute_jobs),str(i*no_of_impute_jobs+distance)])
            individual_prefix=prefix + '_'+ str(i)
            individual_command.extend(['-o',individual_prefix+'.haps','-w',individual_prefix + '.warnings','-i',individual_prefix +'.info'])
            cmds.append(individual_command)
        return(self.get_dependency(step_name),cmds)       
         

    def run_impute2(self,options,config,haps,dependencies):
        (cmd_template,prefix) = CommandTemplate.run_impute2(self,options,config,haps)
        memory_required="20"
        output_dependencies=[]
        try:
            proc =subprocess.Popen( """tail -1 {0} > awk '{{print $2}}'""".format(join(options.vcf_input),stdout=subprocess.PIPE,shell=True))
        except:
            logger.error("Tail command failed on VCF input")
            sys.exit(SUBPROCESS_FAILED_EXIT)
        distance=int(self.config['impute2']['chromosome_split_size'])*1000000
        no_of_impute_jobs = int(proc.stdout.read())//distance + 1
        (cmds,create_impute_split_step_name)=self.create_impute_commands(prefix,haps,dependencies)
        for cmd in cmds:
            memory_required=""
            wall_time="11:59:00"
            step_name = prefix + self.get_date_string()
            self.write_step_preamble(memory_required,wall_time,step_name,create_impute_split_step_name)
            self.serial_task(1)
            self.queue_ll()
            self.write_job_preamble(memory_required,cmd)
            output_dependencies.append(self.get_dependency(step_name))
        self.join_impute2_files(output_dependencies,no_commands,output_dependencies)
        return(output_dependencies,output_prefix+'.haps')

    def run_plink_filter(self,options, config,ped,map,dependencies):
        logger.debug("Preparing plink filtering for running on pan")
        (cmd,prefix) = CommandTemplate.run_plink_filter(self,options,config,ped,map)
        memory_required="4"
        wall_time="11:59:00"
        step_name = prefix + self.get_date_string()
        self.write_step_preamble(memory_required,wall_time,step_name,dependencies)
        self.serial_task(1)
        self.queue_ll()
        self.write_job_preamble(memory_required,cmd)
        logger.debug("Finished preparing plink filtering for running on pan")
        return(self.get_dependency(step_name),prefix+'.ped',prefix+'.map')
    
    # Utility functions used for setting up the job step names
    def string_to_bytes(self,input):
        return bytes(input,'UTF-8')

    def get_date_string(self):
        return str(datetime.now()).replace(' ','').replace(':','').replace('.','').replace('-','')
    
