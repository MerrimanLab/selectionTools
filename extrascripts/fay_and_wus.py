import argparse
from selection_pipeline import *
import logging

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
parser = argparse.ArgumentParser(description='Generate fay and wus data using variscan')
parser.add_argument('--fay-window-width',dest='fayandWuWindowWidth',help="Sliding window width for Fay and Wu's")
parser.add_argument('--fay-window-jump',dest='fayandWuWindowJump',help="Window JJump for Fay and Wus(if fay-window-width = fay-window-jump non-overlapping windows are used")
parser.add_argument('-i',dest="haps_file",help="Haps input file")
parser.add_argument('-s',dest="sample_input_file",help="Sample input file")
parser.add_argument('--config-file',dest='config_file')
parser.add_argument('-o',dest="output_prefix",help="Prefix for output files")
parser.add_argument('-c',dest="chromosome",help="Chromosome to perform fay and wus on")
args=parser.parse_args()
if(args.fayandWuWindowWidth is None):
    args.fayandWuWindowWidth = str(5000)
if(args.fayandWuWindowJump is None):
    args.fayandWuWindowJump = str(5000);
if(args.config_file is None):
    args.config_file = 'defaults.cfg'
if(args.output_prefix is None):
    args.output_prefix = 'FAY'

# Uses the pipeline to run just fay and wus.
config  = parse_config(args)
s=StandardRun(args,config,full_run=False)
haps = s.run_aa_annotate_haps(args.haps_file)
new_sample_file = s.fix_sample_file(args.sample_input_file)
haps2 = s.prepare_haps_for_variscan(haps,new_sample_file)
fayandwus= s.variscan_fayandwus(haps2)
logging.info("Fay and Wus Done")
logging.info(fayandwus)
logging.info("Done :)")


