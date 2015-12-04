# Multipopulation script calls the selection

# pipeline for each population that we need
# to do then zips up and runs a script to p# each of the cross population
# statistics once
# it has all finished.
# institution: University of Otago
# author: James Boocock
#
#
# requires that the selection pipeline
# installed.
#

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict
import math
import sys
import os
from optparse import OptionParser
import configparser
import logging
from .environment import set_environment
from .standard_run_utilities import *
logger = logging.getLogger(__name__)

SUBPROCESS_FAILED_EXIT = 10
CANNOT_FIND_EXECUTABLE = 20
CANNOT_FIND_CONFIG = 30

# generate RSB after we have calculated ihs


def rsb(config, options, populations):
    """ Runs RSB script

        Reads config options and populations and generates
        rsb statistics for each population pair.
    """
    rscript = config['Rscript']['rscript_executable']
    generate_rsb = config['Rscript']['generate_rsb']
    directory = 'rsb'
    if not os.path.exists(directory):
        os.mkdir(directory)
    pops = list(populations.keys())
    orig_dir = os.getcwd()
    os.chdir(directory)
    for i in range(0, len(pops)-1):
        cmd = []
        pop1 = pops[i]
        cmd.append(rscript)
        pop1_ihh_file = os.path.join(orig_dir, pop1, 'results',
                                     pop1 + 'chr' + options.chromosome +
                                     '.ihh')
        cmd.extend([generate_rsb, '--pop1', pop1, '--pop1file', pop1_ihh_file])
        for j in range(i+1, len(pops)):
            tmp_cmd = []
            tmp_cmd.extend(cmd)
            pop2 = pops[j]
            pop2_ihh_file = os.path.join(orig_dir, pop2, 'results', pop2 +
                                         'chr' + options.chromosome + '.ihh')
            tmp_cmd.extend(['--pop2', pop2, '--pop2file', pop2_ihh_file])
            tmp_cmd.extend(['--chr', options.chromosome])
            run_subprocess(tmp_cmd, 'rsb_generation')
    os.chdir(orig_dir)


def get_populations(populations):
    """ Returns the populations and ids

        Imports data from the input population files and creates a dictionary
        of the populations.
    """
    pops = {}
    for pop in populations:
        with open(pop, 'r') as f:
            for i, line in enumerate(f):
                line = line.strip()
                if (i == 0):
                    pop_name = line
                    pops[pop_name] = []
                else:
                    pops[pop_name].append(line)
    return pops


def parse_config(options):
    """ Parse config file

        Read the config file and save the results in a dictionary
    """
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
                config_parsed[section][op] = config.get(section, op)
            except:
                logger.info("exception on {0}".format(op))
                config_parsed[section][op] = None
    return config_parsed


def check_executables_and_scripts_exist(options, config):
    """ Check the executables actually exist where specified.

        Uses the config file to determine whether the executable
        are at the location as expeceted.
    """
    if(which(config['vcftools']['vcf_subset_executable'],
             'vcf-subset') is None):
        return False
    if(which(config['selection_pipeline']['selection_pipeline_executable'],
             'selection_pipeline') is None):
        return False
    return True


def subset_vcf(vcf_input, config, populations):
    """ Run subset VCF to break the VCF file into populations

        Uses the VCF input file and the population dictionary
        to run VCF-subset uses parallelizations as specified
        by the number of cores.
    """
    vcf_outputs = []
    vcf_dict = {}
    no_pops = len(populations)
    threads = int(config['system']['cores_avaliable'])
    threads_per_job = int(math.ceil(threads / float(no_pops)))
    line_count = get_vcf_line_count(vcf_input)
    split_length = line_count // threads_per_job
    split_positions = [split_length * i for i in range(1, threads_per_job+1)]
    remainder_length = line_count % threads
    split_positions[len(split_positions) - 1] += remainder_length
    vcf_inputs = split_vcf(vcf_input, split_positions)
    cmds = []
    stdouts = []
    for i, vcf in enumerate(vcf_inputs):
        for key, value in list(populations.items()):
            cmd = []
            output_file = key + str(i) + '.vcf'
            try:
                vcf_dict[key].append(output_file)
            except KeyError:
                vcf_dict[key] = [output_file]
            comma_list_ids = ','.join(value)
            vcf_subset_executable = config['vcftools']['vcf_subset_executable']
            cmd.append(vcf_subset_executable)
            cmd.extend(['-f', '-c', comma_list_ids, vcf])
            stdouts.append(output_file)
            cmds.append(list(cmd))
    queue_jobs(cmds, 'vcf-subset',
               config['system']['cores_avaliable'], stdouts=stdouts)
    cmds = []
    for key, value in list(vcf_dict.items()):
        # generate the commands for vcf concat for each output file generated
        cmd = []
        output_file = key + '.vcf'
        # Append to vcf_outputs
        vcf_outputs.append(output_file)
        if(len(value) == 1):
            os.rename(value[0], output_file)
        else:
            vcf_concat_executable = config['vcftools']['vcf_concat_executable']
            cmd.append(vcf_concat_executable)
            cmd.extend(value)
            cmds.append(list(cmd))
    if(len(cmds) != 0):
        queue_jobs(cmds, 'vcf-concat', config['system']['cores_avaliable'],
                   stdouts=vcf_outputs)
  # call the queue jobs to run vcf-subset
    # return the population concatenated vcf file
    return vcf_outputs


def run_selection_pipeline(output_vcfs, options, populations, config):
    """ Runs the selection_pipeline script for each population

        Uses the population dictionary and the output vcfs from the subset
        process to run the selection pipeline on each population.
    """
    cores = config['system']['cores_avaliable']
    parralelise_populations = False
    # Arbitrary cut off for parralelising each population
    # 4 at the moment could be calculated given the amount
    # of parralelisation needed in each run.
    if(len(populations) >= 4 and int(cores) >= 4):
        parralelise_populations = True
        cores_per_run = str(int(cores) // len(populations))
    else:
        cores_per_run = cores
    orig_dir = os.getcwd()
    if(options.extra_args is not None):
        extra_args = options.extra_args
    else:
        extra_args = ''
    if options.cores is not None:
        extra_args += ' --cores ' + cores_per_run
    # Run the selection pipeline for a single run job #
    selection_pipeline_executable = \
        config['selection_pipeline']['selection_pipeline_executable']
    cmds = []
    # check whether we should disable rsb given that iHS could have potentially
    # been disabled and if that is the case we cannot perform rsb calculation
    if options.extra_args is not None:
        if '--no-ihs' in options.extra_args:
            options.no_rsb = True
    if parralelise_populations:
        folder_names = []
    for vcf, population_name in zip(sorted(output_vcfs), sorted(populations)):
        cmd = []
        cmd.append(selection_pipeline_executable)
        cmd.extend(['-c', options.chromosome, '-i', os.path.abspath(vcf),
                   '--population', population_name,
                   '--config-file', os.path.abspath(options.config_file)])
        cmd.extend(extra_args.split())
        cmds.append(cmd)
        directory = population_name
        if not os.path.exists(directory):
            os.mkdir(directory)
        if parralelise_populations:
            folder_names.append(directory)
        else:
        # Create directory for each sub population to run in
            os.chdir(directory)
            run_subprocess(cmd, 'selection_pipeline')
            os.chdir(orig_dir)
    if parralelise_populations:
        queue_jobs(cmds, 'selection_pipeline',
                   cores, folder_names=folder_names)


def fst_vcf(input_vcf, config, options, populations):
    """ Generates FST statistics for every pair of populations

        Uses the population dictionary to generate a weir-fst statistics
        using VCF-TOOLS.
    """
    vcf_tools = config['vcftools']['vcf_tools_executable']
    directory = 'fst'
    if not os.path.exists(directory):
        os.mkdir(directory)
    pops = list(populations.keys())
    orig_dir = os.getcwd()
    os.chdir(directory)
    for i in range(0, len(pops)-1):
        p = pops[i]
        cmd = []
        cmd_hapmap = [] 
        cmd.append(vcf_tools)
        first_pop_name = open('first_pop.tmp', 'w')
        first_pop_name.write('\n'.join(populations[p]))
        first_pop_name.close()
        cmd.extend(['--fst-window-size', options.fst_window_size, 
                    '--fst-window-step', options.fst_window_step,
                    '--vcf',input_vcf])
        cmd_hapmap.extend(cmd)
        cmd.extend(['--weir-fst-pop', 'first_pop.tmp'])
        cmd_hapmap.extend(['--hapmap-fst-pop','first_pop.tmp'])
        for j in range(i+1, len(pops)):
            s = pops[j]
            tmp_cmd_hapmap = []
            tmp_cmd = []
            tmp_cmd_hapmap.extend(cmd_hapmap)
            tmp_cmd.extend(cmd)
            tmp_cmd.extend(['--weir-fst-pop', 'second_pop.tmp'])
            tmp_cmd_hapmap.extend(['--hapmap-fst-pop','second_pop.tmp'])
            second_pop_name = open('second_pop.tmp', 'w')
            second_pop_name.write('\n'.join(populations[s]))
            second_pop_name.close()
            run_subprocess(tmp_cmd, 'fst_calculation_weir')
            run_subprocess(tmp_cmd_hapmap,'fst_calculation_hapmap')
            os.rename('out.windowed.weir.fst',
                      options.chromosome + p + s + '.weir.fst')
            os.rename('out.windowed.hapmap.fst',
                      options.chromosome + p + s +'.hapmap.fst')
    os.remove('second_pop.tmp')
    os.remove('first_pop.tmp')
    os.remove('out.log')
    os.chdir(orig_dir)


def main():
    """ Main function for multi_population

        Reads config and options and runs the multi_population
        pipeline
    """
    parser = OptionParser()
    parser.add_option('-p', '--population', action='append',
                      dest="populations", help='population_files')
    parser.add_option('-a', '--arguments-selection-pipelines',
                      dest="extra_args", help=('Arguments to the selection'
                                               'pipeline script'))
    parser.add_option('-l', '--log-file', dest="log_file", help="Log file")
    parser.add_option('-i', '--vcf-input-file', dest="vcf_input",
                      help="VCF Input File")
    parser.add_option('-c', '--chromosome', dest="chromosome",
                      help=("Chromosome label doesn't actually have to"
                            "correspond to the real chromosome but is required"
                            " to determine what output files to make"))
    parser.add_option('--config-file', dest='config_file',
                      help='Configuration File')
    parser.add_option('--fst-window-size', dest="fst_window_size",
                      help="FST window size (kb)")
    parser.add_option('--fst-window-step', dest="fst_window_step",
                      help="FST window step size (kb)")
    parser.add_option('--no-clean-up', dest="no_clean_up",
                      action="store_true",
                      help="Do not clean up intermediate datafiles")
    parser.add_option('--cores', dest="cores", help=("Overrides number of "
                      "cores avaliable as provided in the config file"))
    parser.add_option('--no-rsb',dest="no_rsb", action="store_true",
                      help="Do not calculate RSB")
    (options, args) = parser.parse_args()
    print((options.extra_args))
    assert options.vcf_input is not None, \
        "no VCF file has been specified as input"
    assert os.path.isfile(options.vcf_input), \
        "Cannot locate vcf file at path = {0)".format(options.vcf_input)
    assert options.chromosome is not None, \
        "no chromosome has been specified to the script"
    assert options.populations is not None and \
        len(options.populations) >= 2, \
        "At least two population files are required"
    if options.config_file is None:
        options.config_file = 'defaults.cfg'
        if not(os.path.isfile(options.config_file)):
            raise Exception("Cannot find config file")
    elif not(os.path.isfile(options.config_file)):
        raise Exception("Cannot find config file")
    config = parse_config(options)
    if options.log_file is None:
        options.log_file = 'multi_population.log'
    logging.basicConfig(format='%(asctime)s %(message)s',
                        filename=options.log_file, filemode='w',
                        level=logging.INFO)
    if not (check_executables_and_scripts_exist(options, config)):
        sys.exit(CANNOT_FIND_EXECUTABLE)
    if options.no_clean_up is None:
        options.clean_up_files = False
    if options.fst_window_step is None:
        options.fst_window_step = str(1000)
    else:
        options.fst_window_step = str(
                float(options.fst_window_step) * 1e3)
    if options.fst_window_size is None:
        options.fst_window_size = str(1000)
    else:
        options.fst_window_size = str(
                float(options.fst_window_size) * 1e3)
    if options.no_rsb is None:
        options.no_rsb = False
    if options.cores is not None:
        config['system']['cores_avaliable'] = options.cores
    set_environment(config['environment'])
    options.vcf_input = os.path.abspath(options.vcf_input)
    populations = get_populations(options.populations)
    populations = OrderedDict(sorted(list(populations.items()), key=lambda t: t[0]))
    fst_vcf(options.vcf_input, config, options, populations)
    output_vcfs = subset_vcf(options.vcf_input, config, populations)
    run_selection_pipeline(output_vcfs, options, populations, config)
    # TODO move FST to here on filtered dataset
    if not (options.no_rsb):
        rsb(config, options, populations)
    if not os.path.exists('logs'):
        os.mkdir('logs')
    os.rename(options.log_file, 'logs/' + options.log_file)
    if not options.no_clean_up:
        keep = [os.path.basename(options.vcf_input),os.path.basename(options.config_file)]
        keep.extend(options.populations)
        clean_folder('.', keep=keep)
    logger.info("Multi_population Complete")
    logger.info("Goodbye :")
    print("Multi-population selection pipeline completed successfully !:)")
if __name__ == "__main__":
    main()
