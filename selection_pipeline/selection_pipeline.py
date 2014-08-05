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
from optparse import OptionParser
import ConfigParser
import logging
from .standard_run import StandardRun
from .environment import set_environment
logger = logging.getLogger(__name__)
SUBPROCESS_FAILED_EXIT = 10


def parse_config(options):
    """ Parse config file

        Reads a config and parses the
        arguments into a dictionary.
    """
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
    parser.add_option('-v', '--debug',
                      action="store_true", dest='debug',
                      help="Print debug messages")
    parser.add_option('-q', '--silent', action="store_false",
                      dest='verbose', help="Run Silently")
    parser.add_option('-i', '--vcf',
                      dest='vcf_input', help="VCF input file")
    parser.add_option('-c', '--chromosome',
                      dest='chromosome', help="Chromosome")
    parser.add_option('-l', '--log-fire', dest='log_file',
                      help="Log file for the pipeline process")
    parser.add_option('--maf', dest='maf',
                      help='Minor allele-frequency filter')
    parser.add_option('--hwe', dest='hwe',
                      help="Hardy-Weinberg Equillibrium filter proportion")
    parser.add_option('--remove-missing', dest="remove_missing",
                      help="Remove missing genotypes")
    parser.add_option('--config-file', dest="config_file",
                      help="Config file")
    parser.add_option('--phased-vcf', action="store_true",
                      dest="phased_vcf", help="Phased vcf file")
    parser.add_option('--population', dest="population",
                      help="Population Code ")
    parser.add_option('--imputation', action="store_true",
                      dest="imputation", help="Imputation")
    parser.add_option('--full-process', action="store_true",
                      dest="full_process", help="Run Entire Process")
    parser.add_option('--gzvcf', action="store_true",
                      dest="vcf_gz", help="VCF input is in GZ file (optional)")
    parser.add_option('--TajimaD', dest='tajimas_d',
                      help="Output Tajima's D statistic in bins of size (bp)")
    parser.add_option('--fay-Window-Width', dest='fayandWuWindowWidth',
                      help="Sliding window width for Fay and Wu's H (kb)")
    parser.add_option('--fay-Window-Jump', dest="fayandWuWindowJump",
                      help=("Window Jump for Fay and Wus ( if fay-Window-Width"
                            " = fay-Window-Jump non-overlapping windows "
                            "are used (kb)"))
    parser.add_option('--no-clean-up', dest="no_clean_up", action="store_true",
                      help="Do not clean up intermediate datafiles")
    parser.add_option('--impute-split-size', dest='impute_split_size',
                      help="impute2 split size (Mb)")
    parser.add_option('--ehh-window-size', dest="multi_window_size",
                      help="Multicore window size (Mp)")
    parser.add_option('--ehh-overlap', dest="ehh_overlap",
                      help="EHH window overlap (Mb)")
    parser.add_option('--daf', dest='daf',
                      help="Derived Allele Frequency filter proportion")
    parser.add_option('--big-gap', dest="big_gap",
                      help=("Gap size for not calculating iHH if "
                            "core SNP spans this gap (kb)"))
    parser.add_option('--small-gap', dest='small_gap',
                      help=("Gap size for applying a penalty to "
                            "the area calculated by iHH (kb)"))
    parser.add_option('--small-gap-penalty', dest="small_gap_penalty",
                      help=("Penalty multiplier for intergration steps"
                            "in iHH see manual for formula, usually the "
                            "same as small-gap"))
    parser.add_option('--cores', dest='cores',
                      help="Override cores avaliable setting")
    parser.add_option('--no-ihs',dest='no_ihs',action="store_true"
                      , help='Disable iHS and iHH calculation')
    parser.add_option('--haps', dest='haps',
                        help="Shapeit haps file")
    parser.add_option('--sample', dest='sample',
                        help='Corresponding sample file to accompany haps')
    parser.add_option('--beagle',dest='beagle',action='store_true',
                      help="Use beagle to phase")
    parser.add_option('--no-gmap',dest="no_genetic_map",action="store_true",
                      help="Do not use a genetic map for the analysis")
    parser.add_option('--physical-ihs',dest="physical_ihs",help="Use physical map for calculating iHS",action="store_true")
    parser.add_option("--no-plots" , dest="no_plots", action="store_true",
                      help="Do not create rudimentary plots")
    (options, args) = parser.parse_args()
    if(options.verbose is not None):
        if(options.debug):
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.ERROR)
    # Obligatory arguments
    assert options.vcf_input or (options.haps and options.sample) is not None, \
        "No VCF or haps/sample file has been specified as input"
    assert options.chromosome is not None, \
        "No chromosome has been specified to the script"
    assert options.population is not None, \
        "Population code has not been specified."
    assert options.config_file is not None, \
        "Config file has not been specified."
    assert os.path.isfile(options.config_file), \
        "Config file cannot be found on path {0}".format(options.config_file)
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
