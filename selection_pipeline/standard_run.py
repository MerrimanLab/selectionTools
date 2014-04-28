import os
import sys
from .standard_run_utilities import *
from .run_pipeline import CommandTemplate

## Subprocess import clause required for running commands on the shell##
import subprocess
import logging
logger = logging.getLogger(__name__)

MISSING_EXECUTABLE = 50


class StandardRun(CommandTemplate):
    def is_script(self, fpath):
        """ Determines if the fpath is a file

        """
        return os.path.isfile(fpath)

    def is_exe(self, fpath):
        """ Determines if the fpath is a executable file

        """
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
#http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python

    def which(self, program, program_name):
        """ Determines if the specified program or script exist on the path

        """
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
        """ Checks to ensure all the scripts and executables exist.

        Uses the config dictionary to determine whether all the required
        executables and scripts exists where they have been specified.
        """
        if(self.which(
                self.config['plink']['plink_executable'],
                'plink')is None):
            logger.error("plink not found check config file")
            return False
        if(self.which(
                self.config['shapeit']['shapeit_executable'],
                'shapeit') is None):
            logger.error("shapeit not found check config file")
            return False
        if(self.which(
                self.config['ancestral_allele']['ancestral_allele_script'],
                'ancestral_allele'
                ) is None):
            logger.error("ancestral_allele not found check config file")
            return False
        if(self.which(
                self.config['impute2']['impute_executable'],
                'impute2'
                ) is None):
            logger.error("impute2 not found check config file")
            return False
        if(self.which(
                self.config['Rscript']['indel_filter'],
                'indel_filter'
                ) is None):
            logger.error("indel_filter not found check config file")
            return False
        if(self.which(
                self.config['Rscript']['rscript_executable'],
                'Rscript'
                ) is None):
            logger.error("Rscript not found check config file")
            return False
        if(self.which(
                self.config['multicore_ihh']['multicore_ihh'],
                'multicore_ihh'
                ) is None):
            logger.error("multicore_ihh not found check config file")
            return False
        if(self.which(
                self.config['qctool']['qctool_executable'],
                'qctool'
                ) is None):
            logger.error("qctool not found check config file")
            return False
        if(self.which(
                self.config['haps_scripts']['haps_to_hapmap_script'],
                'haps_to_hapmap'
                ) is None):
            logger.error("haps_to_hapmap not found check config file")
            return False
        if(self.which(
                self.config['haps_scripts']['haps_filter_script'],
                'haps_filter'
                ) is None):
            logger.error('haps_filter not found check config file')
            return False
        return True

    def __init__(self, options, config, full_run=True):
        """ Constructor for StandardRun class

        """
        super(StandardRun, self).__init__(options, config)
        self.options = options
        self.config = config
        if(full_run):
            if(not self.check_executables_and_scripts_exist()):
                sys.exit(MISSING_EXECUTABLE)
            self.threads = self.config['system']['cores_avaliable']

    def run_pipeline(self):
        """ Run pipeline runs the pipeline for a standard run

        """
        vcf = self.run_remove_indels_from_vcf()
        if(self.options.phased_vcf):
            (haps, sample) = self.vcf_to_haps(vcf)
        elif(self.options.haps and self.options.sample):
            haps = self.options.haps
            sample = self.options.sample
        else:
            (ped, map) = self.run_vcf_to_plink()
            (ped, map) = self.run_plink_filter(ped, map)
            (haps, sample) = self.run_shape_it(ped, map)
        if(self.options.imputation):
            (haps) = self.run_impute2(haps)
            haps = self.indel_filter(haps)
        haps = self.haps_filter(haps)
        new_sample_file = self.fix_sample_file(sample)
        haps2_haps = self.prepare_haps_for_variscan(haps, new_sample_file)
        fayandwus = self.variscan_fayandwus(haps2_haps)
        vcf = self.haps_to_vcf(haps, new_sample_file)
        vcf = self.fix_vcf_qctool(vcf)
        haps = self.run_aa_annotate_haps(haps)
        tajimaSD = self.vcf_to_tajimas_d(vcf)
        if (not self.options.no_ihs):
            ihh = self.run_multi_coreihh(haps)
        ihs_file = ihh.split('.ihh')[0] + '.ihs'
        haplo_hh = ihh.split('.ihh')[0] + '.RData'
        if not os.path.exists('results'):
            os.mkdir('results')
        os.rename(tajimaSD, 'results/' + tajimaSD)
        if (not self.options.no_ihs):
            os.rename(haplo_hh, 'results/' + haplo_hh)
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
            keep = ['selection_stderr.tmp', 'selection_stdout.tmp']
            clean_folder('.', keep=keep)
        logger.info(tajimaSD)
        logger.info(vcf)
        logger.info(haps)
        logger.info(ihh)
        logger.info(fayandwus)
        logger.info("Pipeline completed successfully")
        logger.info("Goodbye :)")

    def run_remove_indels_from_vcf(self):
        """ Run remove indels from vcf using subprocess

        """
        (cmd, output_name) = \
            super(StandardRun, self).run_remove_indels_from_vcf()
        run_subprocess(cmd, 'remove indels')
        return(output_name)

    def vcf_to_haps(self, vcf):
        """ Run vcf to haps usng subprocess

        """
        (cmd, haps, sample) = super(StandardRun, self).vcf_to_haps(vcf)
        run_subprocess(cmd, 'vcf to haps')
        return(haps, sample)

    def run_vcf_to_plink(self):
        """ Run vcf to plink using subprocess

        """
        (cmd, prefix) = super(StandardRun, self).run_vcf_to_plink()
        run_subprocess(cmd, 'vcftools')
        return(prefix + '.ped', prefix + '.map')

    def run_plink_filter(self, ped, map):
        """ Run plink filtering using subprocess

        """
        (cmd, prefix) = super(StandardRun, self).run_plink_filter(ped, map)
        run_subprocess(cmd, 'plink')
        return(prefix+'.ped', prefix+'.map')

    def run_shape_it(self, ped, map):
        """ Run shape it phasing using subprocess

        """
        (cmd, prefix) = super(StandardRun, self).run_shape_it(
            ped, map)
        cmd.extend(['--thread', self.threads])
        run_subprocess(cmd, 'shapeit')
        return(prefix + '.haps', prefix + '.sample')

    def haps_to_vcf(self, haps, new_sample_file):
        """ Run haps to vcf conversion using subprocess

        """
        (cmd, output_name) = super(StandardRun, self).haps_to_vcf(
            haps, new_sample_file)
        run_subprocess(cmd, 'hapstovcf')
        return(output_name)

    def haps_filter(self, haps):
        """ Runs haps filter using subprocess

        """
        (cmd, output_name) = super(StandardRun, self).haps_filter(haps)
        run_subprocess(cmd, 'haps filter')
        return(output_name)

    def join_impute2_files(self, output_prefix, no_commands):
        """ Join the impute2 output files

            Opens all the haps, warnings and info files
            for impute2
        """
        output_haps = open(output_prefix+'.haps', 'w')
        output_warnings = open(output_prefix+'.warnings', 'w')
        output_info = open(output_prefix+'.info', 'w')
        for i in range(no_commands):
            with open(output_prefix+'_'+str(i)+'.haps_haps', 'r') as h:
                with open(output_prefix + '_'+str(i)+'.warnings', 'r') as w:
                    with open(output_prefix + '_'+str(i) + '.info', 'r')as f:
                        output_haps.write(h.read())
                        output_warnings.write(w.read())
                        output_info.write(f.read())
        output_haps.close()
        output_warnings.close()
        output_info.close()

    def run_impute2(self, haps):
        """ Run impute2 using subprocess

            Impute2 splits the files into the region specified in
            the options by default it is the recommended 5 megabases.

            These chunks are then run in parralel if possible given
            the users cores setting.
        """
        (cmd_template, output_prefix) = super(StandardRun, self).run_impute2(
            haps)
        distance = int(self.options.impute_split_size)
        try:
            proc = subprocess.Popen("""tail -1 {0}| awk '{{print $3}}'"""
                                    .format(haps), stdout=subprocess.PIPE,
                                    shell=True)
        except:
            logger.error("Tail command failed on haps file")
            sys.exit(SUBPROCESS_FAILED_EXIT)
        try:
            head = subprocess.Popen("""head -1 {0}| awk '{{print $3}}'"""
                                    .format(haps), stdout=subprocess.PIPE,
                                    shell=True)
        except:
            logger.error("Head command failed on haps file")
            sys.exit(SUBPROCESS_FAILED_EXIT)
        start_position = int(head.stdout.read())
        no_of_impute_jobs = ((int(proc.stdout.read())-int(start_position))
                             // distance + 1)
        first_window = start_position // distance
        cmds = []
        for i in range(0, no_of_impute_jobs):
            individual_command = list(cmd_template)
            individual_command.extend(['-int', str((i+first_window)*distance),
                                      str((i+first_window+1)*distance)])
            individual_prefix = output_prefix + '_' + str(i)
            individual_command.extend(['-o', individual_prefix+'.haps',
                                      '-w', individual_prefix + '.warnings',
                                      '-i', individual_prefix + '.info'])
            cmds.append(list(individual_command))
        queue_jobs(cmds, 'impute2', self.config['system']['cores_avaliable'])
        self.join_impute2_files(output_prefix, no_of_impute_jobs)
        return(output_prefix+'.haps')

    def indel_filter(self, haps):
        """ Run the indel filter script using subprocess

        """
        (cmd, output_name) = super(StandardRun, self).indel_filter(haps)
        run_subprocess(cmd, 'indel_filter')
        return(output_name)

    def run_aa_annotate_haps(self, haps, vcf=False):
        """ Run ancestral annotation using subprocess

        """
        if(vcf):
            (cmd, output_name, sample_name) = (super(StandardRun, self)
                                               .run_aa_annotate_haps(
                                                   haps, vcf))
            run_subprocess(cmd, 'ancestral_annotation')
            return(output_name, sample_name)
        else:
            (cmd, output_name) = super(StandardRun, self).run_aa_annotate_haps(
                haps)
            run_subprocess(cmd, 'ancestral_annotation')
            return(output_name)

    def run_multi_coreihh(self, haps):
        """ Run multicore ihh using subprocess

            Uses the number of cores to specify the parallelisation
            to multi_core ihh
        """
        (cmd, output_name) = super(StandardRun, self).run_multi_coreihh(haps)
        cores = self.threads
        ihs_output = output_name.split('.ihh')[0]+'.ihs'
        rdata_output = output_name.split('.ihh')[0]+'.RData'
        cmd.extend(['--cores', cores])
        cmd.extend(['--working_dir', '.'])
        cmd.extend(['--offset', '1'])
        cmd.extend(['--ihs'])
        run_subprocess(cmd, 'multcore_ihh')
        os.rename(self.options.population+'_chr_'+self.options.chromosome +
                  "_wd_"+'.'+"_.ihh", output_name)
        os.rename(self.options.population+'_chr_'+self.options.chromosome +
                  '_wd_'+'.'+"_.ihs", ihs_output)
        os.rename(self.options.population+'_chr_'+self.options.chromosome +
                  '_wd_'+'.'+"_.RData", rdata_output)
        return output_name

    def fix_sample_file(self, sample_file):
        """ Run fix sample file using subprocess

        """
        (cmd, output_name) = super(StandardRun, self).fix_sample_file(
            sample_file)
        new_sample_file = open(output_name, 'w')
        run_subprocess(cmd, 'fix sample file', stdout=new_sample_file)
        new_sample_file.close()
        return(output_name)

    def vcf_to_tajimas_d(self, vcf):
        """ Run vcf to tajimas D using subprocess

        """
        (cmd, output_name) = super(StandardRun, self).vcf_to_tajimas_d(vcf)
        run_subprocess(cmd, 'tajimas_d')
        taj_file = self.options.population + self.options.chromosome + '.taj_d'
        os.rename(output_name, taj_file)
        return(taj_file)

    def fix_vcf_qctool(self, vcf):
        """ Run qctool to vcf fix using subprocess

        """
        (cmd, output_name) = super(StandardRun, self).fix_vcf_qctool(vcf)
        fixed_vcf = open(output_name, 'w')
        run_subprocess(cmd, 'fix vcf qctool', stdout=fixed_vcf)
        fixed_vcf.close()
        return(output_name)

    def prepare_haps_for_variscan(self, haps, sample):
        """ Run prepare haps file for variscan using subprocess

        """
        (cmd, output_name) = (super(StandardRun, self).
                              prepare_haps_for_variscan(haps, sample))
        run_subprocess(cmd, 'haps to hapmap')
        return(output_name)

    def variscan_fayandwus(self, hap2):
        """ Run fay and wus variscan

            Gets the start and the end of the haps file
            to set the StartPos and EndPos variables in the
            variscan input file
        """
        (cmd, output_name, varscan_conf) = (super(StandardRun, self)
                                            .variscan_fayandwus(hap2))
        output_variscan = open(output_name, 'w')
        varscan_config = open(varscan_conf, 'a')
        try:
            proc = subprocess.Popen(
                """tail -1 {0}| awk '{{print $4}}'""".
                format(hap2), stdout=subprocess.PIPE,
                shell=True)
        except:
            logger.error("Tail command failed on haps file")
            sys.exit(SUBPROCESS_FAILED_EXIT)
        try:
            head = subprocess.Popen(
                """head -2 {0} | tail -1 | awk '{{print $4}}'""".
                format(hap2), stdout=subprocess.PIPE, shell=True)
        except:
            logger.error("Head command failed on haps file")
            sys.exit(SUBPROCESS_FAILED_EXIT)
        start_pos = head.stdout.read()
        start_position = int(start_pos)
        end_position = int(proc.stdout.read())
        varscan_config.write('StartPos = ' + str(start_position) + "\n")
        varscan_config.write('EndPos = ' + str(end_position) + '\n')
        varscan_config.close()
        run_subprocess(cmd, 'variscan', stdout=output_variscan)
        output_variscan.close()
        return(output_name)
