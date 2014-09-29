from selection_pipeline.haps_filters import\
    hardy_weinberg_asymptotic, filter_haps_file, stats
import doctest
import unittest
import selection_pipeline
import os
from selection_pipeline.selection_pipeline import parse_config
from selection_pipeline.run_pipeline import CommandTemplate
from selection_pipeline.aa_annotate import aa_seq, write_sample_file, \
    get_haps_line, aa_check

import vcf

suite = doctest.DocTestSuite(selection_pipeline)

def get_file(fname):
    return os.path.join(os.path.dirname(__file__),fname)

class Args:
    """ Represents command line arguments
    """
    pass    
class TestHapsFilter(unittest.TestCase):
    
    def test_hardy_exact(self):
        p_value = hardy_weinberg_asymptotic(138,1469,5)

    def test_hardy_asymptotic(self):
        p_value = hardy_weinberg_asymptotic(138,1469,5)
        assert p_value == 0.34263933319103679
    def test_filter_hap_file(self):
        args = Args()
        args.haps = get_file('filter.haps')
        args.output = get_file('output.haps')
        args.maf = 0.05
        args.missing = 0.05
        args.hwe = .05
        args.chi_square = True
        filter_haps_file(args)
        with open(args.output) as f:
            lines = sum(1 for line in f)
        assert lines == 1
        with open(args.output) as f:
            line = f.readline()
            line = line.split()
            assert line[0] == '2'
            assert line[1] == 'rs4662641'
            assert line[2] == '130000272' 
        os.remove(args.output)

    def test_remove_triallelic(self):
        args = Args()
        args.haps = get_file('triallelic_haps.haps')
        args.output = get_file('output.haps')
        args.maf = 0.05
        args.missing = 0.80
        args.hwe = .05
        args.chi_square = True
        filter_haps_file(args)
        with open(args.output) as f:
            lines = sum(1 for line in f) 
        assert lines == 1
        os.remove(args.output)

class TestRunPipeline(unittest.TestCase):
   
    def setUp(self): 
        self.options = Args()
        self.options.config_file = get_file('defaults.cfg')
        self.config = parse_config(self.options)
        self.options.output_prefix = 'CEU'
        self.options.chromosome = '5'
        self.options.vcf_input = 'testcase.vcf'
        self.options.population = 'CEU'
        self.template = CommandTemplate(self.options,self.config)

    #def test_run_impute2(self):
    #    (cmd_template, prefix) = self.template.run_impute2('test.haps')
    #    assert prefix == \
    #        self.options.output_prefix + self.options.chromosome +\
    #        '_impute2'
    #    assert cmd_template[0] == '/home/smilefreak/selectionTools/bin/impute2'
    #    assert len(cmd_template) == 10

    def test_remove_indels_vcf(self):
        (cmd,output_name) = self.template.run_remove_indels_from_vcf()
        assert output_name == 'testcase.recode.vcf'
        assert len(cmd) == 7
    
class TestAncestralAnnotation(unittest.TestCase):

    def __verify_sequence__(self,aaSeq):
        assert ''.join(aaSeq) == 'GCCG'

    def test_aa_seq_single_chromosome(self):
        options = Args()
        options.ancestralfasta = get_file('ancestor.fa')
        options.single_chromosome = True
        aaSeq = aa_seq(options)
        self.__verify_sequence__(aaSeq)

    def test_aa_seq_header_regex(self):
        options = Args()
        options.ancestralfasta = get_file('ancestor.fa') 
        options.header = 'ANCESTOR_?_FA'
        options.chromosome = '2'
        options.single_chromosome = False
        aaSeq = aa_seq(options)
        self.__verify_sequence__(aaSeq)

    def test_write_sample_file(self):
        options = Args()
        options.sample_file = get_file('test_sample.sample')
        vcf_reader = vcf.Reader(filename=get_file('CEU_test.vcf'))
        write_sample_file(options,vcf_reader)
        sample_names = []
        with open(options.sample_file) as f:
            for i, sample in enumerate(f):
                if( i < 2):
                    continue
                else:
                    sample_name = sample.split()[0]
                    sample_names.append(sample_name)
        for s1, s2 in zip(sorted(sample_names),sorted(vcf_reader.samples)):
            assert s1 == s2
        os.remove(options.sample_file)

    def test_get_haps_line(self):
        options = Args()
        vcf_reader = vcf.Reader(filename=get_file("CEU_test.vcf"))
        record = vcf_reader.next()
        line = get_haps_line(options,record)
        line = line.split()
        assert line[0] == "rs147096179"
        assert line[2] == "130000004"
        assert line[3] == 'C'
        assert line[4] == 'T' 
    
    def test_aa_check(self):
        options = Args()
        realAA= 'G'
        ref = 'C'
        alt = 'T'
        format = "lower"
        line = '2 rs1000 1 C T 0 1'
        # Test ancestral allele != ref
        new_line = aa_check(realAA,ref,alt,format,line)
        for item in new_line.split()[5:]:
            assert item == '1' 
        # Test ancestral allele == ref
        realAA = 'C'
        new_line = aa_check(realAA,ref,alt,format,line)
        assert new_line == '2 rs1000 1 C T 0 1'
        realAA = 'T'
        new_line = aa_check(realAA,ref,alt,format,line)
        assert new_line == '2 rs1000 1 T C 1 0' 
        
            
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestHapsFilter))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestRunPipeline))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestAncestralAnnotation))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestHapsFilter))
