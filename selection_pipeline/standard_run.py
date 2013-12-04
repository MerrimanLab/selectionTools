import os
import sys
from .standard_run_utilities import *
from .run_pipeline import CommandTemplate

## Subprocess import clause required for running commands on the shell##
import subprocess
import logging
logger = logging.getLogger(__name__)


class StandardRun(CommandTemplate):
    def is_script(self, fpath):
        return os.path.isfile(fpath)

    def is_exe(self, fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    #http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python

    def which(self, program, program_name):
        fpath, fname = os.path.split(program)
        if fpath:
            if self.is_exe(program):
                return program
            elif (self.is_script(program)):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if self.is_exe(exe_file):
                    return exe_file
        logger.error(program_name + " path = " + fpath +
                     " not locatable path or in the directory"
                     " specified in your self.config file")
        return None

    def check_executables_and_scripts_exist(self):
        if(self.which(
                self.config['plink']['plink_executable'],
                'plink')is None):
            return False
        if(self.which(
                self.config['shapeit']['shapeit_executable'],
                'shapeit') is None):
            return False
        if(self.which(
                self.config['ancestral_allele']['ancestral_allele_script'],
                'ancestral_allele'
                ) is None):
            return False
        if(self.which(
                self.config['impute2']['impute_executable'],
                'impute2'
                ) is None):
            return False
        if(self.which(
                self.config['Rscript']['indel_filter'],
                'indel_filter'
                ) is None):
            return False
        if(self.which(
                self.config['Rscript']['rscript_executable'],
                'Rscript'
                ) is None):
            return False
        if(self.which(
                self.config['multicore_ihh']['multicore_ihh'],
                'multicore_ihh'
                ) is None):
            return False
        if(self.which(
                self.config['qctool']['qctool_executable'],
                'qctool'
                ) is None):
            return False
        return True

    def __init__(self, options, config, full_run=True):
        #Perform local executable check.#
        self.options = options
        self.config = config
        if(full_run):
            if(not self.check_executables_and_scripts_exist()):
                sys.exit(MISSING_EXECUTABLE_ERROR)
            self.threads = self.config['system']['cores_avaliable']

    def run_pipeline(self):
        if(self.options.phased_vcf):
            (haps, sample) = self.run_aa_annotate_haps(
                self.options.vcf_input, vcf=True)
            haps = self.indel_filter(haps)
            new_sample_file = self.fix_sample_file(sample)
            ihh = self.run_multi_coreihh(haps)
            tajimaSD = self.vcf_to_tajimas_d(vcf_input)
            haps2_haps = self.prepare_haps_for_variscan(new_sample_file)
            fayandwus = self.variscan_fayandwus(haps2_haps)
        else:
            (ped, map) = self.run_vcf_to_plink()
            (ped, map) = self.run_plink_filter(ped, map)
            (haps, sample) = self.run_shape_it(ped, map)
            if(self.options.imputation):
                (haps) = self.run_impute2(haps)
            haps = self.indel_filter(haps)
            new_sample_file = self.fix_sample_file(sample)
            vcf = self.haps_to_vcf(new_sample_file)
            vcf = self.fix_vcf_qctool(vcf)
            haps = self.run_aa_annotate_haps(haps)
            haps2_haps = self.prepare_haps_for_variscan(new_sample_file)
            fayandwus = self.variscan_fayandwus(self.config, haps2_haps)
            tajimaSD = self.vcf_to_tajimas_d(vcf)
        ihh = self.run_multi_coreihh(haps)
        ihs_file = ihh.split('.ihh')[0]+'.ihs'
        haplo_hh = ihh.split('.ihh')[0] +'.haphh'
        if not os.path.exists('results'):
            os.mkdir('results')
        os.rename(tajimaSD, 'results/' + tajimaSD)
        os.rename(vcf, 'results/' + vcf)
        os.rename(ihh, 'results/' + ihh)
        os.rename(ihs_file, 'results/'+ihs_file)
        os.rename(haps, 'results/' + haps)
        os.rename(fayandwus, 'results/' + fayandwus)
        if not os.path.exists('log'):
            os.mkdir('log')
        logger.info(self.options.log_file)
        os.rename(self.options.log_file, 'log/' + self.options.log_file)
        if not self.options.no_clean_up:
            clean_folder('.')
        logger.info(tajimaSD)
        logger.info(vcf)
        logger.info(haps)
        logger.info(ihh)
        logger.info(fayandwus)
        logger.info("Pipeline completed successfully")
        logger.info("Goodbye :)")

    def run_vcf_to_plink(self):
        (cmd, prefix) = CommandTemplate.run_vcf_to_plink(
            self, self.options, self.config)
        run_subprocess(cmd, 'vcftools') 
        return(prefix + '.ped', prefix + '.map')

    def run_plink_filter(self, ped, map):
        (cmd,prefix) = CommandTemplate.run_plink_filter(
            self, self.options, self.config, ped, map)
        run_subprocess(cmd,'plink')
        return(prefix+'.ped',prefix+'.map')

    # Calls a subprocess to run shape it

    def run_shape_it(self,ped,map):
        (cmd,prefix) = CommandTemplate.run_shape_it(
            self,self.options,self.config,ped,map)
        cmd.extend(['--thread',self.threads])
        run_subprocess(cmd,'shapeit')
        return(prefix + '.haps',prefix + '.sample')

    def haps_to_vcf(self,haps,new_sample_file):
        (cmd,output_name) = CommandTemplate.haps_to_vcf(
            self,self.options,self.config,haps,new_sample_file)
        run_subprocess(cmd,'hapstovcf')
        return(output_name) 
    
             
    def join_impute2_files(self,output_prefix,no_commands):
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
    def run_impute2(self,haps):
        (cmd_template,output_prefix) = CommandTemplate.run_impute2(
               self,self.options,self.config,haps)

        # change from megabases to bp which is what is
        # expected by the impute2 command line self.options
        distance=int(self.options.impute_split_size) * 1000000
        # Break files into 5 megabase regions.
        try:
            proc = subprocess.Popen("""tail -1 {0}| awk '{{print $3}}'""".format(haps),stdout=subprocess.PIPE,shell=True) 
        except:
            logger.error("Tail command failed on haps file")
            sys.exit(SUBPROCESS_FAILED_EXIT)
        # get the start of the haps file 
        try:
            head = subprocess.Popen("""head -1 {0}| awk '{{print $3}}'""".format(haps),stdout=subprocess.PIPE,shell=True) 
        except:
            logger.error("Head command failed on haps file")
            sys.exit(SUBPROCESS_FAILED_EXIT)
        #Get the max position from your haps file# 
        start_position = int(head.stdout.read())
        no_of_impute_jobs = (int(proc.stdout.read())-int(start_position))//distance + 1
        #create the command template
        #Get the max position from your haps file# 

        # get the start of the first window
        # 
        first_window = start_position // distance 
        cmds = []
        for i in range(0,no_of_impute_jobs):
            individual_command=list(cmd_template)
            individual_command.extend(['-int',str((i+first_window)*distance),str((i+first_window+1)*distance)])
            individual_prefix=output_prefix + '_'+ str(i)
            individual_command.extend(['-o',individual_prefix+'.haps','-w',individual_prefix + '.warnings','-i',individual_prefix +'.info'])
            cmds.append(list(individual_command))
        queue_jobs(cmds,'impute2',self.config['system']['cores_avaliable'])
        join_impute2_files(self.options,self.config,output_prefix,no_of_impute_jobs)
        return(output_prefix+'.haps') 
         

    def indel_filter(self,haps):
        (cmd,output_name) = CommandTemplate.indel_filter(self,self.options,self.config,haps)
        run_subprocess(cmd,'indel_filter')
        return(output_name)
    
    def run_aa_annotate_haps(self,haps,vcf=False):
        if(vcf):
            (cmd,output_name,sample_name) = CommandTemplate.run_aa_annotate_haps(self,self.options,self.config,haps,vcf)
            run_subprocess(cmd,'ancestral_annotation')
            return(output_name,sample_name)
        else:
            (cmd,output_name) = CommandTemplate.run_aa_annotate_haps(self,self.options,self.config,haps)
            run_subprocess(cmd,'ancestral_annotation')
            return(output_name)
    
    def run_multi_coreihh(self,haps):
        (cmd,output_name) = CommandTemplate.run_multi_coreihh(self,self.options,self.config,haps)
        cores=self.threads
        ihs_output=output_name.split('.ihh')[0]+'.ihs'
        cmd.extend(['--cores',cores])
        cmd.extend(['--working_dir','.'])
        cmd.extend(['--offset','1'])
        cmd.extend(['--ihs'])
        run_subprocess(cmd,'multcore_ihh')
        os.rename(self.options.population+'_chr_'+self.options.chromosome+"_wd_"+'.'+"_.ihh",output_name)
        os.rename(self.options.population+'_chr_'+self.options.chromosome+'_wd_'+'.'+"_.ihs",ihs_output)
        return output_name
    def fix_sample_file(self,sample_file):
        (cmd,output_name) = CommandTemplate.fix_sample_file(self,self.options,self.config,sample_file)
        new_sample_file=open(output_name,'w')
        run_subprocess(cmd,'fix sample file',stdout=new_sample_file)
        new_sample_file.close()
        return(output_name) 
    def vcf_to_tajimas_d(self,vcf):
        (cmd,output_name) = CommandTemplate.vcf_to_tajimas_d(self,self.options,self.config,vcf)
        run_subprocess(cmd,'tajimas_d')
        taj_file =self.options.population + self.options.chromosome + '.taj_d'
        os.rename(output_name,taj_file)
        return(taj_file)
    def fix_vcf_qctool(self,vcf):
        (cmd,output_name) = CommandTemplate.fix_vcf_qctool(self,self.options,self.config,vcf)
        fixed_vcf = open(output_name,'w')
        run_subprocess(cmd,'fix vcf qctool',stdout=fixed_vcf)
        fixed_vcf.close()
        return(output_name)  
    def prepare_haps_for_variscan(self,haps,sample):   
        (cmd,output_name) = CommandTemplate.prepare_haps_for_variscan(self,self.options,self.config,haps,sample)
        run_subprocess(cmd,'haps to hapmap') 
        return(output_name)
    def variscan_fayandwus(self,hap2):
        (cmd,output_name,varscan_conf) = CommandTemplate.variscan_fayandwus(self,self.options,self.config,hap2) 
        output_variscan = open(output_name,'w')
         
        # get the start and the end of the hapmap2 file
        #
        varscan_config = open(varscan_conf ,'a')
        try:
            proc = subprocess.Popen("""tail -1 {0}| awk '{{print $4}}'""".format(hap2),stdout=subprocess.PIPE,shell=True) 
        except:
            logger.error("Tail command failed on haps file")
            sys.exit(SUBPROCESS_FAILED_EXIT)
        # get the start of the haps file 
        try:
            head = subprocess.Popen("""head -2 {0} | tail -1 | awk '{{print $4}}'""".format(hap2),stdout=subprocess.PIPE,shell=True) 
        except:
            logger.error("Head command failed on haps file")
            sys.exit(SUBPROCESS_FAILED_EXIT)
        start_pos = head.stdout.read()
        start_position = int(start_pos)
        end_position = int(proc.stdout.read())
        
        # Write the start and end of the Fey and Wu's into the self.config file #
        varscan_config.write('StartPos = ' + str(start_position) + "\n")
        varscan_config.write('EndPos = ' + str(end_position)+ '\n')
        varscan_config.close() 
        run_subprocess(cmd,'variscan',stdout=output_variscan)
        output_variscan.close()
        return(output_name)
