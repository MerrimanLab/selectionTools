import re
import os
import fnmatch

from optparse import OptionParser

import logging
logger = logging.getLogger(__name__)

class CommandTemplate(object):

    def run_vcf_to_plink(self,options,config):
        cmd = []
        prefix = options.output_prefix + options.chromosome
        logger.debug("Attempting to call vcf tools to convert to ped/map plink format")
        vcf_tools=config['vcftools']['vcf_tools_executable']
        cmd.append(vcf_tools)
        cmd.extend(['--gzvcf',options.vcf_input, '--plink', '--out',prefix,'--remove-indels'])
        logger.debug(config['vcftools']['extra_args'])
        cmd.extend(config['vcftools']['extra_args'].split())
        return (cmd,prefix)

    def run_plink_filter(self,options,config,ped,map):
        
        logger.debug("Attempting to call vcf tools to convert to ped/map plink format")
        cmd = []
        prefix = ped.split('.')[0]
        plink = config['plink']['plink_executable']
        cmd.append(plink)
        # add standard plink commands #

        cmd.extend(['--noweb','--file',prefix,'--geno',str(options.remove_missing),'--hwe',str(options.hwe),'--maf',str(options.maf),'--recode','--out',prefix])
        cmd.extend(config['plink']['extra_args'].split())
        return(cmd,prefix)
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
        cmd.extend(['--input-ped',ped,map,'-M',os.path.join(config['shapeit']['genetic_map_dir'],genetic_map),'--output-max',prefix])
        cmd.extend(config['shapeit']['extra_args'].split())
        return(cmd,prefix)


    def indel_filter(self,options,config,haps):
        cmd = []    
        output_name= options.output_prefix + options.chromosome + '_indel_filter.haps'        
        logger.debug('Attempting to run the R indel and maf filter usually reserved for after phasing')
        rscript = config['Rscript']['rscript_executable']
        indel_filter = config['Rscript']['indel_filter']
        cmd.append(rscript)
        cmd.append(indel_filter)
        cmd.extend([haps,str(options.maf),output_name])
    
        return(cmd,output_name)



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
        return (cmd,output_name)

    def run_multi_coreihh(self,options,config,haps):
        cmd=[]
        output_name= options.output_prefix + '.ihh'
        rscript=config['Rscript']['rscript_executable']
        multicore_ihh=config['multicore_ihh']['multicore_ihh']
        window=config['multicore_ihh']['window']
        overlap=config['multicore_ihh']['overlap']
        # Default offset is 0 as this is the single pc operation something different happens on nesi
        population=options.population
        cmd.append(rscript)
        # Todo look at MAF in rehh
        cmd.extend([multicore_ihh,'-p',population,'-i',haps,'-c',str(options.chromosome),'--window',str(window),'--overlap',str(overlap),'--maf',str(config['multicore_ihh']['derived_allele_frequency'])])
        return (cmd,output_name)

