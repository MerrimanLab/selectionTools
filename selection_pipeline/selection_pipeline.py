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
    parser.add_option('-v', '--verbose',
                      action="store_true", dest='verbose',
                      help="Print debug messages")
    parser.add_option('-q', '--silent', action="store_false",
                      dest='verbose', help="Run Silently")
    parser.add_option('-i', '--vcf',
                      dest='vcf_input', help="VCF input file")
    parser.add_option('-o', '--out',
                      dest='output_prefix', help="Output file prefix")
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
                      help="Output Tajima's D statistic in bins of size <int>")
    parser.add_option('--fay-Window-Width', dest='fayandWuWindowWidth',
                      help="Sliding window width for Fay and Wu's H")
    parser.add_option('--fay-Window-Jump', dest="fayandWuWindowJump",
                      help=("Window Jump for Fay and Wus ( if fay-Window-Width"
                            " = fay-Window-Jump non-overlapping windows "
                            "are used"))
    parser.add_option('--no-clean-up', dest="no_clean_up", action="store_true",
                      help="Do not clean up intermediate datafiles")
    parser.add_option('--impute-split-size', dest='impute_split_size',
                      help="impute2 split size (Mb)")
    parser.add_option('--ehh-window-size', dest="multi_window_size",
                      help="Multicore window size (bp)")
    parser.add_option('--ehh-overlap', dest="ehh_overlap",
                      help="EHH window overlap")
    parser.add_option('--daf', dest='daf',
                      help="Derived Allele Frequency filter proportion")
    parser.add_option('--big-gap', dest="big_gap",
                      help=("Gap size in kb for not calculating iHH if "
                            "core SNP spans this gap"))
    parser.add_option('--small-gap', dest='small_gap',
                      help=("Gap size in kb for applying a penalty to "
                            "the area calculated by iHH"))
    parser.add_option('--small-gap-penalty', dest="small_gap_penalty",
                      help=("Penalty multiplier for intergration steps"
                            "in iHH see manual for formula, usually the "
                            "same as small-gap"))
    parser.add_option('--cores', dest='cores',
                      help="Override cores avaliable setting")
    (options, args) = parser.parse_args()
    if(options.verbose is not None):
        if(options.verbose):
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.ERROR)
    # Obligatory arguments
    assert options.vcf_input is not None, \
        "No VCF file has been specified as input"
    assert options.chromosome is not None, \
        "No chromosome has been specified to the script"
    assert options.output_prefix is not None, \
        "Output file prefix has not been specified."
    assert options.population is not None, \
        "Population code has not been specified."
    if(options.fayandWuWindowJump is None):
        options.fayandWuWindowJump = str(5000)
    if(options.fayandWuWindowWidth is None):
        options.fayandWuWindowWidth = str(5000)
    if(options.no_clean_up is None):
        options.no_clean_up = False
    if(options.tajimas_d is None):
        options.tajimas_d = str(5000)
    else:
        options.tajimas_d = str(options.tajimas_d)
    if(options.imputation is None):
        options.imputation = False
    if(options.hwe is None):
        options.hwe = str(0.001)
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
    if(options.log_file is None):
        options.log_file = options.population + \
            options.chromosome + "_selection_pipeline.log"
    if (options.impute_split_size is None):
        options.impute_split_size = str(5)
    # TODO give sensible defaults kilobases ?
    if (options.multi_window_size is None):
        options.multi_window_size = str(5000000)
    if (options.ehh_overlap is None):
        options.ehh_overlap = str(2000000)
    if (options.big_gap is None):
        options.big_gap = str(0)
    if (options.small_gap is None):
        options.small_gap = str(0)
    if (options.small_gap_penalty is None):
        options.small_gap_penalty = str(0)
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


if __name__ == "__main__":
    main()
