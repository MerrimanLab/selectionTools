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

load_leveler_case="""
        case $LOADL_STEP_NAME in
        """
load_leveler_each_case="""
        {0})
            {1} 
        ;;
        """
load_leveler_esac="""
        esac
        """

load_leveler_template="""
        #@ shell = /bin/bash
        #@ group = {0}
        #@ class = {1}
        #@ output = logs/$(step_name).$(jobid).out
        #@ error = logs/$(step_name).$(jobid).err 
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
    # to fix this we need to do a case statement for all the jobs_steps
    
    #We run this to generate a script for every job
    def write_job_preamble(self,memory_required,cmd,job_step,inList=True,modules=None):
        ulimit = str(int(memory_required) * 1024 * 1024)
        command = ''
        command = command + (load_leveler_ulimit.format(ulimit))
        if(modules is not None):
            for module in modules:
                command = command +('module load '+ module) + '\n'
        if(inList):
            command = command +(' '.join(cmd))
        else:
            command = command + (cmd) + '\n'
           
        self.job_commands.append(command)
        self.job_steps.append(job_step)
    def create_multicore_ihh_nesi(self,memory_required,cmd,job_step,inList,module=None):
        return 0 
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
    #Could potentially fuck out a whole part filtered out of your data
    def __init__(self,options,config):
        self.options=options
        self.config=config
        # append to job_steps with all the step names
        # all the commands in a ordered list        
        self.job_commands=[]
        self.job_steps=[]
        logger.debug('Running the script on nesi')
        try:
            os.makedirs('logs')
        except OSError:
            logger.warning("Could not create logs directory folder exists")
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

    def create_case_loadl_case(self):
        self.load_leveler_script.write(self.string_to_bytes(load_leveler_case))
        for step, command in zip(self.job_steps, self.job_commands):
            self.load_leveler_script.write(self.string_to_bytes(load_leveler_each_case.format(step,command)))
        self.load_leveler_script.write(self.string_to_bytes(load_leveler_esac))

    def create_load_leveler_script(self,options,config):
        if(options.phased_vcf): 
            (dependecies,haps) = self.ancestral_annotation_vcf(options,config)
            (dependencies,haps) = self.run_multi_coreihh(options,config,haps,dependencies)
        else:
            (dependencies,ped,map) = self.run_vcf_to_plink(options,config)
            (dependencies,ped,map) = self.run_plink_filter(options,config,ped,map,dependencies)
            print(dependencies)
            (dependencies,haps) = self.run_shape_it(options,config,ped,map,dependencies) 
        if(options.imputation):
            (dependencies,haps)= self.run_impute2(options,config,haps,dependencies)
        (dependencies,haps) = self.indel_filter(options,config,haps,dependencies)
        #tajimas = run_tajimas_d(options,config,haps)
        (dependencies,haps) = self.run_aa_annotate_haps(options,config,haps,dependencies)
        (dependencies,ihh) = self.run_multi_coreihh(options,config,haps,dependencies)
        self.create_case_loadl_case()
        logger.info("Pipeline completed successfully")
        logger.info(haps)
        logger.info(ihh)
        logger.info("Goodbye :)")
    
    # Change multicore ihh so it uses the new arguments
    def run_multi_coreihh(self,options,config,haps,dependencies):
        (cmd,prefix) = CommandTemplate.run_multi_coreihh(self,options,config,haps)
        parralel_cores = '10'
        cmd.extend(['--cores',parralel_cores])
        memory_required=str(3)
        wall_time="00:59:00"
        r_script_path=config['rscript']['rscript_executable']
        #retarded dependecy brought about by my retarded coding
        python_executable='python'
        nesi_script_dir=config['nesi']['nesi_script_dir']
        cmd.extend(['--script_dir',nesi_script_dir,'--python',py_executable,'--rscript',r_script_path])
        step_name = self.get_step_name(prefix)
        self.write_step_preamble(memory_required,wall_time,step_name,dependencies)
        self.serial_task(1)
        self.queue_ll()
        cmd.append('\n' + "mv *.ihh " +prefix)
        self.write_job_preamble(memory_required,cmd,step_name)
        return( self.get_dependency(step_name), prefix + '.haps') 

    def run_aa_annotate_haps(self,options,config,haps,dependencies):
        (cmd,prefix) =CommandTemplate.run_aa_annotate_haps(self,options,config,haps)
        memory_required="4"
        wall_time="11:59:00"
        step_name = self.get_step_name(prefix)
        self.write_step_preamble(memory_required,wall_time,step_name,dependencies)
        #10 threads for shapeit
        self.serial_task(1)
        self.queue_ll()
        modules=[]
        #Make useful. so that its not just a hard coded value#
        modules.append('python/3.3.2')
        self.write_job_preamble(memory_required,cmd,step_name,modules=modules)
        return( step_name + ">=0 ",prefix+'.haps') 

    def indel_filter(self,options,config,haps,dependencies):
        (cmd,prefix) = CommandTemplate.indel_filter(self,options,config,haps)
        memory_required="4"
        wall_time="11:59:00"
        step_name = self.get_step_name(prefix)
        self.write_step_preamble(memory_required,wall_time,step_name,dependencies)
        #10 threads for shapeit
        self.serial_task(1)
        self.queue_ll()
        self.write_job_preamble(memory_required,cmd,step_name)
        return(self.get_dependency(step_name),prefix+'.haps') 

    def run_shape_it(self,options,config,ped,map,dependencies):
        logger.debug("Preparing shapeit for running on nesi")
        (cmd,prefix) = CommandTemplate.run_shape_it(self,options,config,ped,map)
        cmd.extend(['--thread','2'])
        memory_required="10"
        wall_time="00:59:00"
        step_name = self.get_step_name(prefix)
        self.write_step_preamble(memory_required,wall_time,step_name,dependencies)
        #10 threads for shapeit
        self.serial_task(2)
        self.queue_ll()
        self.write_job_preamble(memory_required,cmd,step_name)
        logger.debug("Finished preparing shape it for running on nesi")
        return(self.get_dependency(step_name),prefix + '.haps')
    def run_vcf_to_plink(self,options,config):
        logger.debug("Preparing vcf_to_plink for running on pan")
        (cmd,prefix) = CommandTemplate.run_vcf_to_plink(self,options,config)
        memory_required="4"
        wall_time="11:59:00"
        step_name = self.get_step_name(prefix)
        self.write_step_preamble(memory_required,wall_time,step_name)
        self.serial_task(1)
        self.queue_ll()
        self.write_job_preamble(memory_required,cmd,step_name)
        #and queue
        logger.debug("Finished preparing vcf_to_plink for running on pan")
        return(self.get_dependency(step_name),prefix + '.ped', prefix+'.map',) 
    def join_impute2_files(self,prefix,no_commands,dependencies):
        step_name = self.get_step_name(prefix)
        haps=[]
        warnings=[]
        info=[]
        command=[]
        for i in range(no_commands):
            haps.append(prefix+'_'+str(i)+'.haps')
            warnings.append(prefix+'_'+str(i)+'.warnings')
            info.append(prefix+'_'+str(i)+'.info')
        command = 'cat {0} > {1}'.format(' '.join(haps),prefix + '.haps') + '\n'
        command = command +  'cat {0} > {1}'.format(' '.join(haps), prefix + '.haps') + '\n'
        command = command + 'cat {0} > {1}'.format(' '.join(haps),prefix +'.haps') + '\n'
        memory_required="4"
        wall_time="11:59:00"
        step_name = self.get_step_name(prefix)
        self.write_step_preamble(memory_required,wall_time,step_name," && ".join(dependencies))
        self.serial_task(1)
        self.queue_ll()
        self.write_job_preamble(memory_required,command,step_name,inList=False)
        return(self.get_dependency(step_name))
    def create_impute_commands(self,prefix,haps,cmd_template,no_of_impute_jobs,distance):
        cmds=[]
        print(cmd_template)
        for i in range(0,no_of_impute_jobs):
                        
            individual_command =[]
            individual_command.extend(cmd_template)
            individual_command.extend(['-int',str(i*no_of_impute_jobs),str(i*no_of_impute_jobs+distance)])
            individual_prefix=prefix + '_'+ str(i)
            individual_command.extend(['-o',individual_prefix+'.haps','-w',individual_prefix + '.warnings','-i',individual_prefix +'.info'])
            print(cmd_template)
            cmds.append(individual_command)
        return(cmds)       
         

    def run_impute2(self,options,config,haps,dependencies):
        (cmd_template,prefix) = CommandTemplate.run_impute2(self,options,config,haps)
        memory_required="20"
        output_dependencies=[]
        try:
            proc =subprocess.Popen( """tail -1 {0} | awk '{{print $2}}'""".format(options.vcf_input),stdout=subprocess.PIPE,shell=True)
        except:
            logger.error("Tail command failed on VCF input")
            sys.exit(SUBPROCESS_FAILED_EXIT)
        distance=int(self.config['impute2']['chromosome_split_size'])*1000000
        no_of_impute_jobs = int(proc.stdout.read())//distance + 1
        (cmds)=self.create_impute_commands(prefix,haps,cmd_template,no_of_impute_jobs,distance)
        for cmd in cmds:
            memory_required="20"
            wall_time="11:59:00"
            step_name = self.get_step_name(prefix)
            self.write_step_preamble(memory_required,wall_time,step_name,dependencies)
            self.serial_task(1)
            self.queue_ll()
            self.write_job_preamble(memory_required,cmd,step_name)
            output_dependencies.append(self.get_dependency(step_name))
        step_name=self.join_impute2_files(prefix,no_of_impute_jobs,output_dependencies)
        return(step_name,prefix+'.haps')

    def run_plink_filter(self,options, config,ped,map,dependencies):
        logger.debug("Preparing plink filtering for running on pan")
        (cmd,prefix) = CommandTemplate.run_plink_filter(self,options,config,ped,map)
        memory_required="4"
        wall_time="11:59:00"
        step_name = self.get_step_name(prefix)
        self.write_step_preamble(memory_required,wall_time,step_name,dependencies)
        self.serial_task(1)
        self.queue_ll()
        self.write_job_preamble(memory_required,cmd,step_name)
        logger.debug("Finished preparing plink filtering for running on pan")
        return(self.get_dependency(step_name),prefix+'.ped',prefix+'.map')
    
    # Utility functions used for setting up the job step names
    def string_to_bytes(self,input):
        return bytes(input,'UTF-8')

    def get_step_name(self,prefix):
        return (prefix + str( datetime.now())).replace(' ','').replace(':','').replace('.','').replace('-','')
    
