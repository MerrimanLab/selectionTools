from selection_pipeline.haps_filters import\
    hardy_weinberg_asymptotic, filter_haps_file, stats
import doctest
import unittest
import selection_pipeline
import os
from selection_pipeline.selection_pipeline import parse_config
from selection_pipeline.run_pipeline import CommandTemplate


suite = doctest.DocTestSuite(selection_pipeline)

def get_file(fname):
    return os.path.join(os.path.dirname(__file__),fname)

class Args:
    """ Represents command line arguments
    """
    pass    
class TestHapsFilter(unittest.TestCase):

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
        args.asymptotic = True
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

class TestRunPipeline(unittest.TestCase):
   
    def setUp(self): 
        self.options = Args()
        self.options.config_file = get_file('defaults.cfg')
        self.config = parse_config(self.options)
        self.options.output_prefix = 'CEU'
        self.options.chromosome = '5'
        self.options.vcf_input = 'testcase.vcf'
        self.template = CommandTemplate(self.options,self.config)

    def test_run_impute2(self):
        (cmd_template, prefix) = self.template.run_impute2('test.haps')
        assert prefix == \
            self.options.output_prefix + self.options.chromosome +\
            '_impute2'
        assert cmd_template[0] == '/home/sfk/selectionTools/bin/impute2'
        assert len(cmd_template) == 10

    def test_remove_indels_vcf(self):
        (cmd,output_name) = self.template.run_remove_indels_from_vcf()
        assert output_name == 'testcase.recode.vcf'
        assert len(cmd) == 7
        assert cmd[0] == '/home/sfk/selectionTools/bin/vcftools'

suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestHapsFilter))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestRunPipeline))
