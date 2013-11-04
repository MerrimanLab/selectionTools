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

import sys
import os
import subprocess
from optparse import OptionParser
import configparser
import logging
from .environment import set_environment
logging.basicConfig(format='%(asctime)s %(message)s')
logger = logging.getLogger(__name__)
logger.setLevel(level=logging.INFO)


SUBPROCESS_FAILED_EXIT = 10
CANNOT_FIND_EXECUTABLE = 20
CANNOT_FIND_CONFIG = 30



# using fey and wus h
# generate the statistics for fey and wus h
def fey_and_wus_h():
    return 0
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
    for key, value in list(populations.items()):
        cmd = []
        vcf_output = open(key + '.vcf','w')
        population = key
        comma_list_ids = ','.join(value)
        vcf_merge_exec=config['vcftools']['vcf_subset_executable']
        cmd.append(vcf_merge_exec)
        cmd.extend(['-f','-c',comma_list_ids,vcf_input])
        run_subprocess(cmd,'vcf-merge',stdout=vcf_output)
        vcf_outputs.append(key + '.vcf')
        vcf_output.close()
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
        #run_subprocess.close()
        running_log.close()
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
            print(tmp_cmd)
            run_subprocess(tmp_cmd,'fst_calculation')
            os.rename('out.windowed.weir.fst',options.chromosome + p + s + '.fst')
    os.remove('second_pop.tmp')
    os.remove('first_pop.tmp')        
 
    os.chdir(orig_dir)
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
    if not os.path.isfile(options.config_file):
        logger.error("Cannot find config file specified at path {0}".format(options.config_file))
        sys.exit(CANNOT_FIND_CONFIG)
    config = parse_config(options)
    if not (check_executables_and_scripts_exist(options,config)):
        sys.exit(CANNOT_FIND_EXECUTABLE)
    if options.fst_window_step is None:
        options.fst_window_step = str(10)
    else:
        options.fst_window_step = str(options.fst_window_step)
    
    if options.fst_window_size is None:
        options.fst_window_size = str(1000)
    else:
        options.fst_window_size = str(options.fst_window_size) 
    set_environment(config['environment'])
    options.vcf_input = os.path.abspath(options.vcf_input)
    populations=get_populations(options.populations)
    output_fst =  fst_vcf(options.vcf_input,config,options,populations)
    output_vcfs = subset_vcf(options.vcf_input,config,populations)
    run_selection_pipeline(output_vcfs,options,populations,config)
if __name__=="__main__":main()
