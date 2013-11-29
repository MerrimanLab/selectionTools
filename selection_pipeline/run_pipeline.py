import re
import os
import fnmatch
import sys
import subprocess

from optparse import OptionParser

import logging
logger = logging.getLogger(__name__)
SUBPROCESS_FAILED_EXIT=10
class CommandTemplate(object):

    def run_vcf_to_plink(self,options,config):
        cmd = []
        prefix = options.output_prefix + options.chromosome
        logger.debug("Attempting to call vcf tools to convert to ped/map plink format")
        vcf_tools=config['vcftools']['vcf_tools_executable']
        cmd.append(vcf_tools)
        if(options.vcf_gz):
            cmd.append('--gzvcf')
        else:
            cmd.append('--vcf')
            cmd.extend([options.vcf_input, '--plink', '--out',prefix,'--remove-indels'])
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
            if fnmatch.fnmatch(file,config['shapeit']['genetic_map_prefix'].replace('?',options.chromosome)):
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

    def join_impute2_files(options,config,output_prefix,no_commands):
        output_haps=open(output_prefix+'.haps','w')
        output_warnings=open(output_prefix+'.warnings','w')
        output_info=open(output_prefix+'.info','w')
        for i in range(no_commands):
            with open(output_prefix+'_'+str(i)+'.haps_haps','r') as h:
                with open(output_prefix + '_'+str(i)+'.warnings','r') as w:
                    with open(output_prefix + '_'+str(i) + '.info','r')as f:
                        output_haps.write(h.read())
                        output_warnings.write(w.read())
                        output_info.write(f.read())
        output_haps.close()
        output_warnings.close()
        output_info.close()
                         
    def run_impute2(self,options,config,haps):
        prefix = options.output_prefix + options.chromosome + '_impute2'
        logger.debug('Attempting to call impute2 the data')
        impute2=config['impute2']['impute_executable']
        genetic_map_dir=config['impute2']['impute_map_dir']
        genetic_map=''
        for file in os.listdir(config['impute2']['impute_map_dir']):
            if fnmatch.fnmatch(file,config['impute2']['impute_map_prefix'].replace('?',options.chromosome)):
                genetic_map = os.path.join(config['impute2']['impute_map_dir'],file)
        
        legend_file =''
        for file in os.listdir(config['impute2']['impute_reference_dir']):
            if fnmatch.fnmatch(file,config['impute2']['impute_reference_prefix'].replace('?',options.chromosome)+'.legend'):
                legend_file = os.path.join(config['impute2']['impute_reference_dir'],file)
        
        hap_file = ''
        for file in os.listdir(config['impute2']['impute_reference_dir']):
            if fnmatch.fnmatch(file,config['impute2']['impute_reference_prefix'].replace('?',options.chromosome)+'.hap'):
                hap_file = os.path.join(config['impute2']['impute_reference_dir'],file)
        #create the command template
        cmd_template=[]
        cmd_template.append(impute2)
        cmd_template.extend(['-m',genetic_map,'-h',hap_file,'-l',legend_file,'-known_haps_g',haps,'-phase'])
        return (cmd_template,prefix)
    def get_ancestral_fasta(self,options,config):     
        if('reference_fasta' in config['ancestral_allele'].keys()):
            cmd.append('--ref-fasta')
            ancestral_fasta=config['ancestral_allele']['reference_fasta']
        else:
            for file in os.listdir(config['ancestral_allele']['ancestral_fasta_dir']):
                if fnmatch.fnmatch(file,config['ancestral_allele']['ancestral_prefix'].replace('?',options.chromosome)):
                    ancestral_fasta = os.path.join(config['ancestral_allele']['ancestral_fasta_dir'],file)
        return ancestral_fasta

    def run_aa_annotate_haps(self,options,config,in_file,vcf=False):
        cmd = []
        output_haps =options.output_prefix.split('.haps')[0] +'_aachanged.haps'
        if(vcf):
            output_sample= options.output_prefix.split('.haps')[0]+'_aachanged.sample' 
        py_executable = config['python']['python_executable']
        aa_annotate = config['ancestral_allele']['ancestral_allele_script']
        cmd.append(py_executable)
        cmd.append(aa_annotate)
        ancestral_fasta = self.get_ancestral_fasta(options,config)
        cmd.extend(['-c', options.chromosome, '-o', output_haps,'-a',ancestral_fasta])
        if(vcf):
            cmd.extend(['-v',in_file,'-s',output_sample])
            return(cmd,output_haps,output_sample)
        else:
            cmd.extend(['-i',in_file])
            return(cmd,output_haps)

    def run_multi_coreihh(self,options,config,haps):
        cmd=[]
        output_name= options.output_prefix + 'chr' + options.chromosome+ '.ihh'
        rscript=config['Rscript']['rscript_executable']
        multicore_ihh=config['multicore_ihh']['multicore_ihh']
        window=options.multi_window_size
        overlap=options.ehh_overlap
        # Default offset is 0 as this is the single pc operation something different happens on nesi
        population=options.population
        cmd.append(rscript)
        # Todo look at MAF in rehh
        cmd.extend([multicore_ihh,'-p',population,'-i',haps,'-c',str(options.chromosome),'--window',str(window),'--overlap',str(overlap),'--maf',options.daf])
        cmd.extend(['--big_gap',options.big_gap,'--small_gap',options.small_gap,'--small_gap_penalty',options.small_gap_penalty])
        return (cmd,output_name)
    
    def fix_sample_file(self,options,config,sample_file):
        cmd=[]
        cmd.extend(['cut', '-d',' ', '-f', '1-6', sample_file])
        sample_file=sample_file.split('.sample')[0] + '_fixed.sample'
        return(cmd,sample_file)

    def haps_to_vcf(self,options,config,haps,new_sample_file):
        cmd=[]
        output_name=options.output_prefix + options.chromosome + '.vcf'
        qctool_executable=config['qctool']['qctool_executable']
        cmd.append(qctool_executable)
        ancestral_fasta = self.get_ancestral_fasta(options,config)
        cmd.extend(['-filetype','shapeit_haplotypes','-g',haps,'-s',new_sample_file,'-og',output_name])
        return (cmd,output_name)

    def fix_vcf_qctool(self,options,config,vcf):
        cmd=[]
        output_name=vcf.split('.vcf')[0] + '_fixed.vcf'
        cmd.extend(['sed', 's/^NA/{0}/g'.format(options.chromosome),vcf])
        return(cmd,output_name)  
           
    def vcf_to_tajimas_d(self,options,config,vcf):
        cmd=[]
        output_name = 'out.Tajima.D'   
        vcftools_executable=config['vcftools']['vcf_tools_executable']
        cmd.append(vcftools_executable)
        cmd.extend(['--TajimaD',options.tajimas_d,'--vcf',vcf])
        return(cmd,output_name)
    def prepare_haps_for_variscan(self,options,config,haps,sample):
        cmd=[]
        output_name = options.output_prefix + options.chromosome + '.hapmap'
        haps_executable=config['variscan']['haps_to_hapmap_executable']
        ancestral_fasta = self.get_ancestral_fasta(options,config)
        cmd.append(haps_executable)
        cmd.extend(['-i',haps,'-s',sample,'-o', output_name, '--id','ANCESTOR','-a',ancestral_fasta,'-c',options.chromosome])
        return(cmd,output_name)
    def variscan_fayandwus(self,options,config,hap2):
        cmd=[]
        v_config_name = 'variscan.conf'
        output_name = options.output_prefix + options.chromosome + '.faw' 
        variscan_config = open(v_config_name,'w')
        variscan_executable = config['variscan']['variscan_executable']
        cmd.append(variscan_executable)
        cmd.extend([hap2,v_config_name])
        # generate default config file for variscan
        config_string = 'RefPos = 0 \n'
        config_string += 'RefSeq = 1 \n'
        config_string += 'BlockDataFile = none \n'
        config_string += 'SeqChoice = all \n'
        config_string += 'OutGroup = last \n'
        config_string += 'RunMode = 22 \n'
        config_string += 'IndivNames = \n'
        config_string += 'UseMuts = 1 \n'
        config_string += 'CompleteDeletion = 0 \n'
        config_string += 'FixNum = 0 \n'
        config_string += 'NumNuc = 4 \n'
        config_string += 'SlidingWindow = 1 \n'
        config_string += 'WidthSW = {0} \n'.format(options.fayandWuWindowWidth)
        config_string += 'JumpSW = 5000 \n'.format(options.fayandWuWindowJump)
        config_string += 'WindowType = 0 \n'
        config_string += 'UseLDSinglets = 0 \n'
        variscan_config.write(config_string)
        variscan_config.close()
        return(cmd,output_name,v_config_name)
        
 
