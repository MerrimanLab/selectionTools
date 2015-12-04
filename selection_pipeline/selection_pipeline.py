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
# edited Aug 2015 (MC) to add seed setting and threshold flags
from optparse import OptionParser
from optparse import OptionGroup
import configparser
import logging
import os
import sys
from .standard_run import StandardRun
from .environment import set_environment
from ._version import __version__
logger = logging.getLogger(__name__)
SUBPROCESS_FAILED_EXIT = 10


def parse_config(options):
    """ Parse config file

        Reads a config and parses the
        arguments into a dictionary.
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


def parse_arguments():
    """ Parse the comand line arguments

        read the arguments and set sensible
        default values for the program
    """
    parser = OptionParser()
    debug_options = OptionGroup(parser, "Debug Options")
    faw_options = OptionGroup(parser, "Fay and Wu's Options")
    ihs_options = OptionGroup(parser, "iHS Options")
    filter_options = OptionGroup(parser, "Filtering Options")
    tajimas_options = OptionGroup(parser, "Tajima's D Options")
    impute_options = OptionGroup(parser, "Impute2 Options")
    debug_options.add_option('-v', '--debug',
                      action="store_true", dest='debug',
                      help="Print debug messages")
    debug_options.add_option('-q', '--silent', action="store_false",
                      dest='verbose', help="Run Silently")
    req_options = OptionGroup(parser, "Required Options")
    req_options.add_option('-i', '--vcf',
                      dest='vcf_input', help="VCF input file")
    req_options.add_option('-c', '--chromosome',
                      dest='chromosome', help="Chromosome")
    debug_options.add_option('-l', '--log-fire', dest='log_file',
                      help="Log file for the pipeline process")
    filter_options.add_option('--maf', dest='maf',
                      help='Minor allele-frequency filter')
    filter_options.add_option('--hwe', dest='hwe',
                      help="Hardy-Weinberg Equillibrium filter proportion")
    filter_options.add_option('--remove-missing', dest="remove_missing",
                      help="Remove missing genotypes")
    req_options.add_option('--config-file', dest="config_file",
                      help="Config file")
    parser.add_option('--phased-vcf', action="store_true",
                      dest="phased_vcf", help="Phased vcf file")
    req_options.add_option('--population', dest="population",
                      help="Population Code ")
    impute_options.add_option('--imputation', action="store_true",
                      dest="imputation", help="Imputation")
    parser.add_option('--full-process', action="store_true",
                      dest="full_process", help="Run Entire Process")
    parser.add_option('--gzvcf', action="store_true",
                      dest="vcf_gz", help="VCF input is in GZ file (optional)")
    tajimas_options.add_option('--TajimaD', dest='tajimas_d',
                      help="Output Tajima's D statistic in bins of size (bp)")
    faw_options.add_option('--fay-Window-Width', dest='fayandWuWindowWidth',
                      help="Sliding window width for Fay and Wu's H (kb)")
    faw_options.add_option('--fay-Window-Jump', dest="fayandWuWindowJump",
                      help=("Window Jump for Fay and Wus ( if fay-Window-Width"
                            " = fay-Window-Jump non-overlapping windows "
                            "are used (kb)"))
    parser.add_option('--no-clean-up', dest="no_clean_up", action="store_true",
                      help="Do not clean up intermediate datafiles")
    impute_options.add_option('--impute-split-size', dest='impute_split_size',
                      help="impute2 split size (Mb)")
    ihs_options.add_option('--ehh-window-size', dest="multi_window_size",
                      help="Multicore window size (Mp)")
    ihs_options.add_option('--ehh-overlap', dest="ehh_overlap",
                      help="EHH window overlap (Mb)")
    filter_options.add_option('--daf', dest='daf',
                      help="Derived Allele Frequency filter proportion")
    ihs_options.add_option('--big-gap', dest="big_gap",
                      help=("Gap size for not calculating iHH if "
                            "core SNP spans this gap (kb)"))
    ihs_options.add_option('--small-gap', dest='small_gap',
                      help=("Gap size for applying a penalty to "
                            "the area calculated by iHH (kb)"))
    ihs_options.add_option('--small-gap-penalty', dest="small_gap_penalty",
                      help=("Penalty multiplier for intergration steps"
                            "in iHH see manual for formula, usually the "
                            "same as small-gap"))
    parser.add_option('--cores', dest='cores',
                      help="Override cores avaliable setting")
    ihs_options.add_option('--no-ihs',dest='no_ihs',action="store_true"
                      , help='Disable iHS and iHH calculation')
    parser.add_option('--haps', dest='haps',
                        help="Shapeit haps file")
    parser.add_option('--sample', dest='sample',
                        help='Corresponding sample file to accompany haps')
    parser.add_option('--beagle',dest='beagle',action='store_true',
                      help="Use beagle to phase")
    parser.add_option('--no-gmap',dest="no_genetic_map",action="store_true",
                      help="Do not use a genetic map for the analysis")
    ihs_options.add_option('--physical-ihs',dest="physical_ihs",help="Use physical map for calculating iHS",action="store_true")
    parser.add_option("--no-plots" , dest="no_plots", action="store_true",
                      help="Do not create rudimentary plots")
    debug_options.add_option('--version', dest = "ver", action="store_true",
                      help="Print version info")
    parser.add_option('--set-shapeit-seed',dest="shapeitSeed", action="store_true",
                      help="NOT IMPLEMENTED YET seed to start shapeit with")
    impute_options.add_option('--set-impute-seed',dest='imputeSeed', action='store_true',
                      help='NOT IMPLEMENTED YET seed to start impute2 with')
    parser.add_option_group(req_options)
    parser.add_option_group(filter_options)
    parser.add_option_group(tajimas_options)
    parser.add_option_group(faw_options)
    parser.add_option_group(ihs_options)
    parser.add_option_group(impute_options)
    parser.add_option_group(debug_options)
    (options, args) = parser.parse_args()
    if(options.verbose is not None):
        if(options.debug):
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.ERROR)
    if(options.ver is True):
        print(("Version: {0}".format(__version__)))
        sys.exit(1)

    # Obligatory arguments
    assert options.vcf_input or (options.haps and options.sample) is not None, \
        "No VCF or haps/sample file has been specified as input"
    assert options.chromosome is not None, \
        "No chromosome has been specified to the script"
    assert options.population is not None, \
        "Population code has not been specified."
    assert options.config_file is not None, \
        "Config file has not been specified."
    if(options.haps and options.sample):
        assert os.path.isfile(options.haps), \
                "Cannot locate haps file path = {0}".format(options.haps)
        assert os.path.isfile(options.sample), \
                "Cannot locate sample file path = {0}".format(options.sample)
    elif(options.vcf_input):
        assert os.path.isfile(options.vcf_input), \
                "Cannot locate vcf input file path = {0}".format(options.vcf_input)
    if(options.fayandWuWindowJump is None):
        options.fayandWuWindowJump = str(5000)
    else:
        options.fayandWuWindowJump = str(
            int(float(options.fayandWuWindowJump) * 1e3))
    if(options.fayandWuWindowWidth is None):
        options.fayandWuWindowWidth = str(5000)
    else:
        options.fayandWuWindowWidth = str(
            int(float(options.fayandWuWindowWidth) * 1e3))
    if(options.no_clean_up is None):
        options.no_clean_up = False
    if(options.tajimas_d is None):
        options.tajimas_d = str(5000)
    else:
        options.tajimas_d = str(
            int(float(options.tajimas_d) * 1e3))
    if(options.imputation is None):
        options.imputation = False
    if(options.hwe is None):
        options.hwe = str(0.0001)
    if(options.maf is None):
        options.maf = str(0.01)
    if(options.daf is None):
        options.daf = str(0.00)
    if(options.remove_missing is None):
        options.remove_missing = str(0.99)
    if (options.phased_vcf is None):
        options.phased_vcf = False
    if (options.full_process is None):
        options.full_process = False
    if (options.vcf_gz is None):
        options.vcf_gz = False
    if (options.no_ihs is None):
        options.no_ihs = False
    if(options.log_file is None):
        options.log_file = options.population + \
            options.chromosome + "_selection_pipeline.log"
    if (options.impute_split_size is None):
        options.impute_split_size = str(5000000)
    else:
        options.impute_split_size = str(
            int(float(options.impute_split_size) * 1e6))
    if (options.multi_window_size is None):
        options.multi_window_size = str(int(5*1e6))
    else:
        options.multi_window_size = str(
            int(float(options.multi_window_size) * 1e6))
    if (options.ehh_overlap is None):
        options.ehh_overlap = str(int(2*1e6))
    else:
        options.ehh_overlap = str(
            int(float(options.ehh_overlap) * 1e6))
    if (options.big_gap is None):
        options.big_gap = str(0)
    else:
        options.big_gap = str(
            int(float(options.big_gap) * 1e3))
    if (options.small_gap is None):
        options.small_gap = str(0)
    else:
         options.small_gap = str(
            int(float(options.small_gap) * 1e3))
    if (options.small_gap_penalty is None):
        options.small_gap_penalty = str(0)
    else:
        options.small_gap_penalty = str(
            int(float(options.small_gap_penalty) * 1e3))
    if (options.no_genetic_map):
        # Must set beagle to true becasue shapeit will not
        # Work without a genetic map
        options.beagle = True
    if (options.no_plots is None):
        options.no_plots = False
    if (options.physical_ihs is None):
        options.physical_ihs = False
    return options


def main():
    """ The main function

        Runs the selection pipeline.
    """
    options = parse_arguments()
    config = parse_config(options)
    set_environment(config['environment'])
    if options.cores is not None:
        config['system']['cores_avaliable'] = options.cores
    logging.basicConfig(format='%(asctime)s     %(message)s',
                        filename=options.log_file, filemode='w',
                        level=logging.INFO)
    s = StandardRun(options, config=config)
    s.run_pipeline()
    print("Selection Pipeline Completed Successfully :)!")

if __name__ == "__main__":
    main()
