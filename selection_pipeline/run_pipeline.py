import os
import fnmatch
import logging
logger = logging.getLogger(__name__)
SUBPROCESS_FAILED_EXIT = 10


class CommandTemplate(object):
    """ Represents a full selection pipeline run.

        Creates command_templates for each of the tools
        in the selection pipeline. Any new system that
        uses the selection pipeline should inherit from Command
        Template. Any commands you want to use on a new system
        should call the parent class first to setup the template
        to built on for the specifiec architecture.
        e.g. Standard linux box, Load Leveler or another
        cluster interface.
    """

    def __init__(self, options, config):
        """ Initialises the class variables self.config and self.options.

        """
        self.config = config
        self.options = options

    def run_vcf_to_plink(self):
        """ Template for vcf to plink conversion

            Uses vcf tools to convert a VCF file
            to ped/map plink format.
        """
        cmd = []
        prefix = self.options.output_prefix + self.options.chromosome
        vcf_tools = self.config['vcftools']['vcf_tools_executable']
        cmd.append(vcf_tools)
        if(self.options.vcf_gz):
            cmd.append('--gzvcf')
        else:
            cmd.append('--vcf')
            cmd.extend([self.options.vcf_input, '--plink', '--out',
                       prefix, '--remove-indels'])
            cmd.extend(self.config['vcftools']['extra_args'].split())
        return (cmd, prefix)

    def run_remove_indels_from_vcf(self):
        """ Template for running remove indels from vcf

        """
        cmd = []
        output_name = \
            self.options.vcf_input.split('.vcf')[0]
        vcftools = self.config['vcftools']['vcf_tools_executable']
        cmd.append(vcftools)
        cmd.extend(['--vcf', self.options.vcf_input, '--remove-indels',
                   '--out', output_name, '--recode'])
        return(cmd, output_name + '.recode.vcf')

    def run_plink_filter(self, ped, map):
        """ Template for running the plink filter

            Uses PLINK to filter the ped/map file
            filters HWE, MAF and missing Genotype
            information.
        """
        cmd = []
        prefix = ped.split('.')[0]
        plink = self.config['plink']['plink_executable']
        cmd.append(plink)
        # add standard plink commands #

        cmd.extend(['--noweb', '--file', prefix, '--geno',
                   str(self.options.remove_missing), '--hwe',
                   str(self.options.hwe), '--maf', str(self.options.maf),
                   '--recode', '--out', prefix])
        cmd.extend(self.config['plink']['extra_args'].split())
        return(cmd, prefix)

    def run_shape_it(self, ped, map):
        """ Template for running shapeit

            Sets up the default command for running shapeit
            for phasing genotypes.

            Reads the config file to find the location of the
            genetic_map
        """
        cmd = []
        prefix = self.options.output_prefix + \
            self.options.chromosome + '.phased'
        genetic_map = ''
        for file in os.listdir(self.config['shapeit']['genetic_map_dir']):
            if fnmatch.fnmatch(
                file, self.config['shapeit']['genetic_map_prefix'].replace(
                    '?', self.options.chromosome)):
                genetic_map = file
        shapeit = self.config['shapeit']['shapeit_executable']
        cmd.append(shapeit)
        cmd.extend(['--input-ped', ped, map, '-M',
                   os.path.join(self.config['shapeit']['genetic_map_dir'],
                                genetic_map), '--output-max', prefix])
        cmd.extend(self.config['shapeit']['extra_args'].split())
        return(cmd, prefix)

    def indel_filter(self, haps):
        """ Return a template for running the indel filter

            Sets up the default command for running the Rscript
            indel_filter.
        """
        cmd = []
        output_name = self.options.output_prefix + \
            self.options.chromosome + '_indel_filter.haps'
        rscript = self.config['Rscript']['rscript_executable']
        indel_filter = self.config['Rscript']['indel_filter']
        cmd.append(rscript)
        cmd.append(indel_filter)
        cmd.extend([haps, str(self.options.maf), output_name])
        return(cmd, output_name)

    def run_impute2(self, haps):
        """ Return a template for running impute2

            Sets up a command for impute2 searches the paths
            specified in the config file for the genetic_map
            and the known haps and legend files.
        """
        prefix = self.options.output_prefix + self.options.chromosome + \
            '_impute2'
        impute2 = self.config['impute2']['impute_executable']
        genetic_map = ''
        for file in os.listdir(self.config['impute2']['impute_map_dir']):
            if fnmatch.fnmatch(file, (
                self.config['impute2']['impute_map_prefix'].replace(
                    '?', self.options.chromosome))):
                genetic_map = os.path.join(
                    self.config['impute2']['impute_map_dir'], file)
        legend_file = ''
        for file in os.listdir(self.config['impute2']['impute_reference_dir']):
            if fnmatch.fnmatch(file, (
                self.config['impute2']['impute_reference_prefix'].replace(
                    '?', self.options.chromosome) + '.legend')):
                legend_file = os.path.join(
                    self.config['impute2']['impute_reference_dir'], file)
        hap_file = ''
        for file in os.listdir(self.config['impute2']['impute_reference_dir']):
            if fnmatch.fnmatch(file, (
                self.config['impute2']['impute_reference_prefix'].replace(
                    '?', self.options.chromosome) + '.hap')):
                hap_file = os.path.join(
                    self.config['impute2']['impute_reference_dir'], file)
        #create the command template
        cmd_template = []
        cmd_template.append(impute2)
        cmd_template.extend(['-m', genetic_map, '-h', hap_file, '-l',
                            legend_file, '-known_haps_g', haps, '-phase'])
        return (cmd_template, prefix)

    def get_ancestral_fasta(self):
        """ Get the ancestral fasta file for the pipeline

            Reads the config file and gets the ancestral fasta file
            to be used for ancestral annotation.
        """
        if('reference_fasta' in self.config['ancestral_allele'].keys()):
            ancestral_fasta = \
                self.config['ancestral_allele']['reference_fasta']
        else:
            for file in os.listdir(
                    self.config['ancestral_allele']['ancestral_fasta_dir']):
                if fnmatch.fnmatch(
                        file,
                        self.config['ancestral_allele']['ancestral_prefix'].
                        replace('?', self.options.chromosome)):
                    ancestral_fasta = os.path.join(
                        self.config['ancestral_allele']['ancestral_fasta_dir'],
                        file)
        return ancestral_fasta

    def run_aa_annotate_haps(self, in_file, vcf=False):
        """ Return the template for running ancestral annotation

            runs the ancestral annotation python script to convert
            the haps file to derived / ancestral alleles
        """
        cmd = []
        output_haps = self.options.output_prefix.split('.haps')[0] + \
            '_aachanged.haps'
        if(vcf):
            output_sample = self.options.output_prefix.split('.haps')[0] + \
                '_aachanged.sample'
        py_executable = self.config['python']['python_executable']
        aa_annotate = \
            self.config['ancestral_allele']['ancestral_allele_script']
        cmd.append(py_executable)
        cmd.append(aa_annotate)
        ancestral_fasta = self.get_ancestral_fasta()
        cmd.extend(['-c', self.options.chromosome, '-o',
                   output_haps, '-a', ancestral_fasta])
        if('reference_fasta' in self.config['ancestral_allele'].keys()):
            cmd.append('--ref-fasta')
        if(vcf):
            cmd.extend(['-v', in_file, '-s', output_sample])
            return(cmd, output_haps, output_sample)
        else:
            cmd.extend(['-i', in_file])
            return(cmd, output_haps)

    def run_multi_coreihh(self, haps):
        """ Return the template for running multi_coren ihh

        """
        cmd = []
        output_name = self.options.output_prefix + 'chr' + \
            self.options.chromosome + '.ihh'
        rscript = self.config['Rscript']['rscript_executable']
        multicore_ihh = self.config['multicore_ihh']['multicore_ihh']
        window = self.options.multi_window_size
        overlap = self.options.ehh_overlap
        population = self.options.population
        cmd.append(rscript)
        cmd.extend([multicore_ihh, '-p', population, '-i',
                   haps, '-c', str(self.options.chromosome),
                   '--window', str(window), '--overlap', str(overlap),
                   '--maf', self.options.daf])
        cmd.extend(['--big_gap', self.options.big_gap, '--small_gap',
                   self.options.small_gap, '--small_gap_penalty',
                   self.options.small_gap_penalty, '--haplo_hh'])
        return (cmd, output_name)

    def fix_sample_file(self, sample_file):
        """ Return the template for running fix sample file

            Command just cuts the extra columns from the file
            to conform with the standard sample file input
        """
        cmd = []
        cmd.extend(['cut', '-d', ' ', '-f', '1-6', sample_file])
        sample_file = sample_file.split('.sample')[0] + '_fixed.sample'
        return(cmd, sample_file)

    def haps_to_vcf(self, haps, new_sample_file):
        """ Return the template for running haps to vcf

        """
        cmd = []
        output_name = self.options.output_prefix + \
            self.options.chromosome + '.vcf'
        qctool_executable = self.config['qctool']['qctool_executable']
        cmd.append(qctool_executable)
        cmd.extend(['-filetype', 'shapeit_haplotypes', '-g',
                   haps, '-s', new_sample_file, '-og', output_name])
        return (cmd, output_name)

    def fix_vcf_qctool(self, vcf):
        """ Return the template for running fix vcf qctool

        """
        cmd = []
        output_name = vcf.split('.vcf')[0] + '_fixed.vcf'
        cmd.extend(['sed', 's/^NA/{0}/g'.format(self.options.chromosome), vcf])
        return(cmd, output_name)

    def vcf_to_tajimas_d(self, vcf):
        """ Return the template for running vcf to tajima's D

        """
        cmd = []
        output_name = 'out.Tajima.D'
        vcftools_executable = self.config['vcftools']['vcf_tools_executable']
        cmd.append(vcftools_executable)
        cmd.extend(['--TajimaD', self.options.tajimas_d, '--vcf', vcf])
        return(cmd, output_name)

    def prepare_haps_for_variscan(self, haps, sample):
        """ Return the template for running haps to variscan

        """
        cmd = []
        output_name = self.options.output_prefix + self.options.chromosome + \
            '.hapmap'
        haps_executable = self.config['variscan']['haps_to_hapmap_executable']
        (ancestral_fasta, regex) = self.get_ancestral_fasta()
        cmd.append(haps_executable)
        cmd.extend(['-i', haps, '-s', sample, '-o', output_name, '--id',
                   'ANCESTOR', '-a', ancestral_fasta, '-c',
                   self.options.chromosome])
        return(cmd, output_name)

    def variscan_fayandwus(self, hap2):
        """ Return the template for running variscan fay and wus

        """
        cmd = []
        v_config_name = 'variscan.conf'
        output_name = self.options.output_prefix + self.options.chromosome + \
            '.faw'
        variscan_config = open(v_config_name, 'w')
        variscan_executable = self.config['variscan']['variscan_executable']
        cmd.append(variscan_executable)
        cmd.extend([hap2, v_config_name])
        # generate default self.config file for variscan
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
        config_string += 'WidthSW = {0} \n'.format(
            self.options.fayandWuWindowWidth)
        config_string += 'JumpSW = 5000 \n'.format(
            self.options.fayandWuWindowJump)
        config_string += 'WindowType = 0 \n'
        config_string += 'UseLDSinglets = 0 \n'
        variscan_config.write(config_string)
        variscan_config.close()
        return(cmd, output_name, v_config_name)
