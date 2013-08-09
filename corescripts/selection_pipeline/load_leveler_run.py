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

from optparse import OptionParser
import ConfigParser

## Subprocess import clause required for running commands on the shell##
import subprocess
import logging

#Import standard run 

logging.basicConfig(format='%(asctime)s %(message)s')
logger = logging.getLogger(__name__)
logger.setLevel(level=logging.INFO)

SUBPROCESS_FAILED_EXIT=10

load_leveler_template="""
        #@ shell = /bin/bash
        #@ 
"""


class LoadLevelerRun(object):
    
    """ Load leveler class takes the pipeline and runs the PHD on the nesi pan 
        cluster.
    """
    def __init__(self):
        logger.debug('Running the script on nesi')   
     
        if(options.phased_vcf): 
            haps = self.ancestral_annotation_vcf(options,config)
            ihh = self.run_multi_coreihh(options,config,haps)
        else:
            (ped,map) = self.run_vcf_tools(options,config)
            (ped,map) = self.run_plink(options,config,ped,map)
            (haps) = self.run_shape_it(options,config,ped,map) 
        if(options.imputation):
            (haps)= self.run_impute2(options,config,haps)
        haps = self.indel_filter(options,config,haps)
        #tajimas = run_tajimas_d(options,config,haps)
        haps = self.run_aa_annotate_haps(options,config,haps)
        ihh = self.run_multi_coreihh(options,config,haps)
        logger.info("Pipeline completed successfully")
        logger.info(haps)
        logger.info(ihh)
        logger.info("Goodbye :)")
    def ancestral_annotation_vcf(self,options,config):
        return 1   
    def run_vcf_tools(self,options,config):
        return 1
    def run_plink(self,options, config):
        return 1
