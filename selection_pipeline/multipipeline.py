#
#
# Multipopulation script calls the selection
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

from collections import OrderedDict
import math
import sys
import os
import subprocess
from optparse import OptionParser
import ConfigParser
import logging
import re
from .environment import set_environment
from .standard_run_utilities import *
logger = logging.getLogger(__name__)


SUBPROCESS_FAILED_EXIT = 10
CANNOT_FIND_EXECUTABLE = 20
CANNOT_FIND_CONFIG = 30



# generate RSB after we have calculated ihs

def rsb(config,options,populations):
    rscript = config['Rscript']['rscript_executable']
    generate_rsb =config['Rscript']['generate_rsb']
    directory = 'rsb'
    if not os.path.exists(directory):
        os.mkdir(directory)
    pops = list(populations.keys())
    orig_dir = os.getcwd()
    os.chdir(directory)
    for i in range(0,len(pops)-1):
        cmd=[]
        pop1 = pops[i]
        cmd.append(rscript)
        pop1_ihh_file = os.path.join(orig_dir,pop1,'results', pop1 + 'chr' +options.chromosome + '.ihh')
        cmd.extend([generate_rsb,'--pop1',pop1,'--pop1file',pop1_ihh_file])
        for j in range(i+1,len(pops)):
            tmp_cmd=[]
            tmp_cmd.extend(cmd)
            pop2=pops[j] 
            pop2_ihh_file = os.path.join(orig_dir,pop2,'results', pop2 + 'chr' +options.chromosome + '.ihh')
            tmp_cmd.extend(['--pop2',pop2,'--pop2file',pop2_ihh_file])
            tmp_cmd.extend(['--chr',options.chromosome]) 
            run_subprocess(tmp_cmd,'rsb_generation') 
    os.chdir(orig_dir)     



def get_populations(populations):
    pops = {}
    for pop in populations:
        with open(pop, 'r') as f:
            for i, line in enumerate(f):
                line = line.strip()
                if ( i == 0 ):
                    pop_name = line
                    pops[pop_name]=[]
                else:
                    pops[pop_name].append(line)
    return pops
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

def check_executables_and_scripts_exist(options,config): 
        executables=['vcf-subset','selection_pipeline']
        if(which(config['vcftools']['vcf_subset_executable'],'vcf-subset')is None):
            return False
        if(which(config['selection_pipeline']['selection_pipeline_executable'],'selection_pipeline') is None):
            return False
        return True
    
def __concat__vcfs__(vcf_inputs,config,populations):
    return 0
def subset_vcf(vcf_input,config,populations):
    vcf_outputs = []
    vcf_dict = {}
    no_pops = len(populations)
    threads=config['system']['cores_avaliable']
    # Take threads and divide by the number of jobs.
    threads_per_job = math.ceil(threads / no_pops)
    # get vcf line count
    line_count = get_vcf_line_count(vcf_input)
    # split length is the size of each chunk
    split_length = line_count // threads
    # split positions
    split_positions = [split_length* i for i in range(1,threads+1)]
    remainder_length = line_count % threads 
    split_positions[len(split_positions) - 1] += remainder_length
    vcf_inputs = split_vcf(vcf_input,split_positions)
    cmds=[]
    stdouts=[]
    for i, vcf in enumerate(vcf_inputs):
        for key, value in populations.items():
            cmd = []
            output_file= key +str(i) +'.vcf'  
                vcf_dict[key].append(vcf_output)
            except KeyError:
                vcf_dict[key]=[vcf_output]
            comma_list_ids = ','.join(value)
            vcf_subset_executable=config['vcftools']['vcf_subset_executable']
            cmd.append(vcf_subset_executable)
            cmd.extend(['-f','-c',comma_list_ids,vcf])
            stdouts.append(output_file)
            #run_subprocess(cmd,'vcf-subset',stdout=vcf_output)
            cmds.append(list(cmd))
    queue_jobs(cmds,config['system']['cores_avaliable'],stdouts=stdouts)
    cmds=[]
    for key, value in vcf_dict.items():
        # generate the commands for vcf concat for each output file generated
        cmd=[]
        vcf_concat_executable=config['vcf_tools']['vcf_concat_executable']
        cmd.append(vcf_concat_executable)
        cmd.extend(value)
    queue_jobs(cmds,config['system']['cores_avaliable'],stdouts=vcf_outputs)
    # call the queue jobs to run vcf-subset 
    # return the population concatenated vcf file
    return vcf_outputs 

def run_selection_pipeline(output_vcfs,options,populations,config):
    orig_dir = os.getcwd()
    if(options.extra_args is not None):
        extra_args=options.extra_args
    else:
        extra_args='' 
    # Run the selection pipeline for a single run job #
    selection_pipeline_executable=config['selection_pipeline']['selection_pipeline_executable']
    for vcf, population_name in zip(output_vcfs, populations):
        directory=population_name
        # Create directory for each sub population to run in
        if not os.path.exists(directory):
            os.mkdir(directory)
        
        cmd=[]
        cmd.append(selection_pipeline_executable) 
        cmd.extend(['-c',options.chromosome,'-i',os.path.abspath(vcf),'-o',population_name,'--population',population_name,'--config-file',os.path.abspath(options.config_file)])
        cmd.append(extra_args)  
        os.chdir(directory)
        run_subprocess(cmd,'selection_pipeline')
        #running_log.close()
        os.chdir(orig_dir)
def fst_vcf(input_vcf,config,options,populations):
    vcf_tools =config['vcftools']['vcf_tools_executable']
    directory = 'fst'
    if not os.path.exists(directory):
        os.mkdir(directory)
    pops = list(populations.keys())
    orig_dir = os.getcwd()
    os.chdir(directory)
    for i in range(0,len(pops)-1):
        p = pops[i]
        cmd=[]
        cmd.append(vcf_tools)
        first_pop_name = open('first_pop.tmp','w')
        first_pop_name.write('\n'.join(populations[p]))
        first_pop_name.close()
        cmd.extend(['--fst-window-size',options.fst_window_size,'--fst-window-step',options.fst_window_step,'--weir-fst-pop','first_pop.tmp','--vcf',input_vcf])
        for j in range(i+1,len(pops)):
            s = pops[j]
            tmp_cmd = []
            tmp_cmd.extend(cmd)
            tmp_cmd.extend(['--weir-fst-pop','second_pop.tmp']) 
            second_pop_name = open('second_pop.tmp','w')
            second_pop_name.write('\n'.join(populations[s]))
            second_pop_name.close()
            run_subprocess(tmp_cmd,'fst_calculation')
            os.rename('out.windowed.weir.fst',options.chromosome + p + s + '.fst')
    os.remove('second_pop.tmp')
    os.remove('first_pop.tmp')        
    os.remove('out.log') 
    os.chdir(orig_dir)
def main():
    parser=OptionParser()
    parser.add_option('-p','--population',action='append',dest="populations",help='population_files')
    parser.add_option('-a','--arguments-selection-pipelines',dest="extra_args",help='Arguments to the selection pipeline script')
    parser.add_option('-l','--log-file',dest="log_file",help="Log file")
    parser.add_option('-i','--vcf-input-file',dest="vcf_input",help="VCF Input File")
    parser.add_option('-c','--chromosome',dest="chromosome",help="Chromosome label doesn't actually have to correspond to the real chromosome but is required to determine what output files to make")
    parser.add_option('--config-file',dest='config_file',help='Configuration File')
    parser.add_option('--fst-window-size',dest="fst_window_size",help="FST window size")
    parser.add_option('--fst-window-step',dest="fst_window_step",help="FST window step size")
    parser.add_option('--no-clean-up',dest="no_clean_up",action="store_true",help="Do not clean up intermediate datafiles")
    
    (options,args) = parser.parse_args()
    assert options.vcf_input is not None, "No VCF file has been specified as input"
    assert options.chromosome is not None, "No chromosome has been specified to the script"
    if options.config_file is None:
	options.config_file = 'defaults.cfg'
    config = parse_config(options)
    if not (check_executables_and_scripts_exist(options,config)):
        sys.exit(CANNOT_FIND_EXECUTABLE)
    
    if options.no_clean_up is None:
        options.clean_up_files = False
    if options.fst_window_step is None:
        options.fst_window_step = str(1000)
    else:
        options.fst_window_step = str(options.fst_window_step)
    if options.log_file is None:
        options.log_file = 'multi_population.log'   
    if options.fst_window_size is None:
        options.fst_window_size = str(1000)
    else:
        options.fst_window_size = str(options.fst_window_size) 
    logging.basicConfig(format='%(asctime)s %(message)s',filename=options.log_file,filemode='w',level=logging.INFO)
    set_environment(config['environment'])
    options.vcf_input = os.path.abspath(options.vcf_input)
    populations=get_populations(options.populations)
    populations=OrderedDict(sorted(populations.items(), key=lambda t: t[0]))
    output_fst =  fst_vcf(options.vcf_input,config,options,populations)
    output_vcfs = subset_vcf(options.vcf_input,config,populations)
    run_selection_pipeline(output_vcfs,options,populations,config)
    rsb(config,options,populations)
    if not os.path.exists('logs'):
        os.mkdir('logs')    
    os.rename(options.log_file,'logs/'+options.log_file)
    if not options.no_clean_up:
        clean_folder('.')
    logger.info("Multi_population Complete")
    logger.info("Goodbye :")
if __name__=="__main__":main()
