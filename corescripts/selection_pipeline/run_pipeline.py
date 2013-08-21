import re
import os
import fnmatch

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
        output_haps=open(output_prefix+'.haps','wb')
        output_warnings=open(output_prefix+'.warnings','wb')
        output_info=open(output_prefix+'.info','wb')
        for i in range(no_commands):
            with open(output_prefix+'_'+str(i)+'haps','r') as h:
                with open(output_prefix + '_'+str(i)+'.warnings','r') as w:
                    with open(output_prefix + '_'+str(i) + '.info','r')as f:
                        output_haps.write(h.read())
                        output_warnings.write(w.read())
                        output_info.write(f.read())
        output_haps.close()
        output_warnings.close()
        output_info.close()
                         
    #def run_tajimas_d(options,config,gen,sample):
    
    #def fu_and_wus_h(options,config,gen,sample):

     

    def run_impute2(self,options,config,haps):
        cmds=[]
        prefix = options.output_prefix + options.chromosome + '_impute2'
        logger.debug('Attempting to call impute2 the data')
        impute2=config['impute2']['impute_executable']
        distance=config['impute2']['chromosome_split_size']
        genetic_map_dir=config['impute2']['genetic_map_dir']
        genetic_map=''
        for file in os.listdir(config['impute2']['genetic_map_dir']):
            if fnmatch.fnmatch(file,config['impute2']['genetic_map_chr'].replace('?',options.chromosome)):
                genetic_map = file
        
        legend_file =''
        for file in os.listdir(config['impute2']['impute_reference_dir']):
            if fnmatch.fnmatch(file,config['impute2']['impute_reference_prefix'].replace('?',options.chromosome+'.legend')):
                legend_file = file
        
        hap_file = ''
        for file in os.listdir(config['impute2']['impute_reference_dir']):
            if fnmatch.fnmatch(file,config['impute2']['impute_reference_prefix'].replace('?',options.chromosome+'.hap')):
                hap_file = file
        distance=int(distance) * 1000000
        # Break files into 5 megabase regions.
        try:
            subprocess.call("`tail -1 ${0}| awk '{print $3}'`".format(haps),stdout=subprocess.PIPE) 
        except:
            logger.error("Tail command failed on haps file")
        #Get the max position from your haps file# 
        max_position=int(proc.stdout.read().strip())
        no_of_impute_jobs=max_position//distance + 1
        
        #create the command template
        cmd_template=[]
        cmd_template.append(impute2)
        cmd_template.extend(['-m',genetic_map,'-h',hap_file,'-l',legend_file,'-known_haps_g',haps])
        for i in range(0,no_of_impute_jobs):
            individual_command=cmd_template
            individual_command.extend(['-int',str(i*no_of_impute_jobs),str(i*no_of_impute_jobs+distance)])
            individual_prefix=prefix + '_'+ str(i)
            individual_command.extend('-o',individual_prefix+'.haps','-w',individual_prefix + '.warninigs','-i',individual_prefix +'.info')
            cmds.append(individual_command)
        return (cmds,prefix)
             

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
        output_name= options.output_prefix + 'chr' + options.chromosome+ '.ihh'
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

