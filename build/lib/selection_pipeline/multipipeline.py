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
import subprocess
from optparse import OptionParser
import configparser
import logging
logging.basicConfig(format='%(asctime)s %(message)s')
logger = logging.getLogger(__name__)
logger.setLevel(level=logging.INFO)


SUBPROCESS_FAILED_EXIT = 10
CANNOT_FIND_EXECUTABLE = 20
CONFIG_FILE_NOT_FOUND = 30

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
    if not os.path.isfile(options.config_file):
        logger.error("Config file not found please check arguments")
        sys.exit(CONFIG_FILE_NOT_FOUND)  
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

def run_subprocess(command,tool,stdout=None):  
        try:
            if(stdout is None):
                exit_code = subprocess.call(command) 
            else:
                exit_code = subprocess.call(command,stdout=stdout)
        except:
            logger.error(tool + " failed to run " + ' '.join(command))
            sys.exit(SUBPROCESS_FAILED_EXIT)   

        if(exit_code != 0): 
            sys.exit(SUBPROCESS_FAILED_EXIT)   
        logger.error("Finished tool " + tool)


def check_executables_and_scripts_exist(options,config): 
        executables=['vcf-subset','selection_pipeline']
        if(which(config['vcftools']['vcf_subset_executable'],'vcf-subset')is None):
            return False
        if(which(config['selection_pipeline']['selection_pipeline_executable'],'selection_pipeline') is None):
            return False
        return True


# Copied from standard run to check whether the programs exist. So much copy pasta.

def is_script(fpath):
    return os.path.isfile(fpath)
def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
#Stolen code from 
#http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
def which(program,program_name):
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
        elif (is_script(program)):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    logger.error(program_name +" path = " + fpath+" not locatable path or in the directory specified in your config file ")
    return None

def subset_vcf(vcf_input,config,populations):
    print(populations)
    vcf_outputs = []
    for key, value in populations.items():
        cmd = []
#        vcf_output = open(key + '.vcf','w')
        population = key
        comma_list_ids = ','.join(value)
        vcf_merge_exec=config['vcftools']['vcf_subset_executable']
        cmd.append(vcf_merge_exec)
        cmd.extend(['-f','-c',comma_list_ids,vcf_input])
  #      run_subprocess(cmd,'vcf-merge',stdout=vcf_output)
        vcf_outputs.append(key + '.vcf')
 #       vcf_output.close()
    return vcf_outputs 

def run_selection_pipeline(output_vcfs,options,populations,config):
    if(options.extra_args is not None):
        extra_args=options.extra_args
    else:
        extra_args='' 
    # Run the selection pipeline for a single run job #
    selection_pipeline_executable=config['selection_pipeline']['selection_pipeline_executable']
    for vcf, population_name in zip(output_vcfs, populations):
        directory=population_name
        # Create directory for each sub population to run in
        print(population_name)
        if not os.path.exists(directory):
            os.mkdir(directory)
        running_log= open(os.path.join(directory,population_name+'.log'),'w')
        
        cmd=[]
        cmd.append(selection_pipeline_executable) 
        cmd.extend(['-c',options.chromosome,'-i',os.path.abspath(vcf),'-o',population_name,'--population',population_name,'--config-file',os.path.abspath(options.config_file)])
        cmd.append(extra_args)  
        os.chdir(directory)
        run_subprocess(cmd,'selection_pipeline')
        run_subprocess.close()

def main():
    parser=OptionParser()
    parser.add_option('-p','--population',action='append',dest="populations",help='population_files')
    parser.add_option('-a','--arguments-selection-pipelines',dest="extra_args",help='Arguments to the selection pipeline script')
    parser.add_option('-i','--vcf-input-file',dest="vcf_input",help="VCF Input File")
    parser.add_option('-c','--chromosome',dest="chromosome",help="Chromosome label doesn't actually have to correspond to the real chromosome but is required to determine what output files to make")
    parser.add_option('-C','--config-file',dest='config_file',help='Configuration File')
    parser.add_option('--fst-window-size',dest="fst_window_size",help="FST window size")
    parser.add_option('--fst-window-step',dest="fst_window_step",help="FST window step size")
    (options,args) = parser.parse_args()
    assert options.vcf_input is not None, "No VCF file has been specified as input"
    assert options.chromosome is not None, "No chromosome has been specified to the script"
    assert options.config_file is not None, "No config file has been specified for the program"
    config = parse_config(options)
    if not (check_executables_and_scripts_exist(options,config)):
        sys.exit(CANNOT_FIND_EXECUTABLE)
     
    populations=get_populations(options.populations)
    output_vcfs = subset_vcf(options.vcf_input,config,populations)
    print(output_vcfs) 
    run_selection_pipeline(output_vcfs,options,populations,config)

if __name__=="__main__":main()
