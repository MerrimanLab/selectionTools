
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
logging.basicConfig(format='%(asctime)s %(message)s')
logger = logging.getLogger(__name__)
logger.setLevel(level=logging.INFO)

class StandardRun(object):
    
    def __init__(options,config):
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
          
 
    def run_vcf_tools(self,config):
        cmd = []
        prefix = options.output_prefix + options.chromosome
        logger.debug("Attempting to call vcf tools to convert to ped/map plink format")
        vcf_tools=config['vcftools']['vcf_tools_executable']
        cmd.append(vcf_tools)
        cmd.extend(['--gzvcf',options.vcf_input, '--plink', '--out',prefix,'--remove-indels'])
        logger.debug(config['vcftools']['extra_args'])
        cmd.extend(config['vcftools']['extra_args'].split())
        try:
            subprocess.call(cmd) 
        except:
            logger.error("Vcf tools failed to run" + ' '.join(cmd))
            sys.exit(SUBPROCESS_FAILED_EXIT)
     
        logger.debug("Finished vcf tools run")
        return(prefix + '.ped', prefix + '.map')

    # Calls a subprocess to run plink

    def run_plink(self,options,config,ped,map):
        logger.debug("Attempting to call vcf tools to convert to ped/map plink format")
        cmd = []
        prefix = ped.split('.')[0]
        plink = config['plink']['plink_executable']
        cmd.append(plink)
        # add standard plink commands #

        cmd.extend(['--noweb','--file',prefix,'--geno',str(options.remove_missing),'--hwe',str(options.hwe),'--maf',str(options.maf),'--recode','--out',prefix])
        cmd.extend(config['plink']['extra_args'].split())
        try:
            subprocess.call(cmd)
        except:
            logger.error("plink failed to run" + ' '.join(cmd))
            sys.exit(SUBPROCESS_FAILED_EXIT)
        logger.error("Finished plink filtering")
        return(prefix+'.ped',prefix+'.map')

        # Calls a subprocess to run shape it

    def run_shape_it(self,options,config,ped,map):
        cmd = []
        prefix = options.output_prefix + options.chromosome + '.phased'
        logger.debug("Attempting to call shape it to phase the data")
        genetic_map = ''
        for file in os.listdir(config['shapeit']['genetic_map_dir']):
            if fnmatch.fnmatch(file,'genetic_map_chr'+options.chromosome+'*'):
                genetic_map = file
            
        shapeit=config['shapeit']['shapeit_executable']
        cmd.append(shapeit)
        cmd.extend(['--input-ped',ped,map,'-M',os.path.join(config['shapeit']['genetic_map_dir'],genetic_map),'--output-max',prefix,'--thread',config['system']['threads_avaliable']])
        cmd.extend(config['shapeit']['extra_args'].split())
        logger.debug(cmd)
        try:
            subprocess.call(cmd)
        except:
            logger.error("shapeit failed to run" + ' '.join(cmd))
        sys.exit(SUBPROCESS_FAILED_EXIT)
        logger.debug('Shape it phasing has completed')
        return(prefix + '.haps')

    #Calls a subprocess to run impute   
 
    def run_impute2(self,options,config,haps):
        cmd = []
        prefix = options.output_prefix + options.chromosome + '_impute2'
        logger.debug('Attempting to call impute2 the data')
        impute2=config['impute2']['impute_executable']
        cmd.append(impute2)
    #    cmd.append(['--
    
    #Calls a subprocess to run the indel filter

    def indel_filter(self,options,config,haps):
        cmd = []    
        output_name= options.output_prefix + options.chromosome + '_indel_filter.haps'        
        logger.debug('Attempting to run the R indel and maf filter usually reserved for after phasing')
        rscript = config['Rscript']['rscript_executable']
        indel_filter = config['Rscript']['indel_filter']
        cmd.append(rscript)
        cmd.append(indel_filter)
        cmd.extend([haps,str(options.maf),output_name])
        try:
            subprocess.call(cmd)
        except:
            logger.error("maf filter failed to run" + ' '.join(cmd))
            sys.exit(SUBPROCESS_FAILED_EXIT)
    
        return(output_name)
 
    def run_aa_annotate_haps(self,options,config,haps):
        cmd = []
        output_name= options.output_prefix + '_aachanged.haps'
        py_executable = config['python']['python_executable']
        aa_annotate = config['ancestral_allele']['ancestral_allele_script']
        logger.debug('Attempting to run ancestral allele annotation')
        for file in os.listdir(config['ancestral_allele']['ancestral_fasta_dir']):
            if fnmatch.fnmatch(file,config['ancestral_allele']['ancestral_prefix'].replace('?',options.chromosome)):
                ancestral_fasta = file
        cmd.append(py_executable)
        cmd.append(aa_annotate)
        cmd.extend(['-i',haps ,'-c', options.chromosome, '-o', output_name,'-a',os.path.join(config['ancestral_allele']['ancestral_fasta_dir'],ancestral_fasta)])
        try:
            subprocess.call(cmd)
        except:
            logger.error("ancestral allele annotation failed to run" + ' '.join(cmd))
            sys.exit(SUBPROCESS_FAILED_EXIT)
        return output_name
    
#def run_tajimas_d(options,config,gen,sample):

    def run_multi_coreihh(self,options,config,haps):
        cmd=[]
        output_name= options.output_prefix + '.ihh'
        rscript=config['Rscript']['rscript_executable']
        multicore_ihh=config['multicore_ihh']['multicore_ihh']
        window=config['multicore_ihh']['window']
        overlap=config['multicore_ihh']['overlap']
        cores=config['system']['threads_avaliable']
        logger.debug("Started running multicore iHH (rehh) script")
        # Default offset is 0 as this is the single pc operation something different happens on nesi
        population=options.population
        cmd.append(rscript)
        # Todo look at MAF in rehh
        cmd.extend([multicore_ihh,population,haps,str(options.chromosome),str(window),str(overlap),str(cores),'.','4',str(config['multicore_ihh']['minor_allele_frequency']),options.output_prefix])
        try:
            subprocess.call(cmd)
        except:
            logger.error("multicore ihh (rehh) failed to run" + ' '.join(cmd))
            sys.exit(SUBPROCESS_FAILED_EXIT)
        os.rename(population+'_chr_'+options.chromosome+"_wd_"+'.'+"_.ihh",output_name)
        return output_name
