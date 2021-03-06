#
# James Boocock and Murray Cadzow
# July 2013 / Sept 2014
# University of Otago
#
import re
import vcf
from optparse import OptionParser
from pyfasta import Fasta
#
# Command Line Arguments
#
# --haps phased haps
# --aa ancestral allele annotation
# --chr Chromosome
# --output Output file
# --format |High|Low|
# optional reference chromosome fasta can be used also
#
# 10000 genomesfasta Ancestral allele file can be downloaded from
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.tar.bz2
#
# Takes a FASTA File Containing the Ancestral Alleles.
# Using the same format as the 1000 genomes ancestral alleles file.
# ATCG high confidence call
# actg low confidence call
# N : failure
# - : the exant species contains a insertion at this position.
# . : no coverage at this alignment.
# annotate haps file.
#Annotate a phased_vcf file that is a file from 1000 genomes
#with the ancestral alleles
# Could potentially take a population list of
#ids also but for now just creates it from the
# vcf file it is given


def aa_seq(options):
    """ Gets the ancestral sequence from a Fasta file

    """
    f = Fasta(options.ancestralfasta)
    keyz = (f.keys())
    match = ''
    if (options.single_chromosome):
        # Single chromosome fasta should only have one sequence.
        # that sequence should be the sequence of interest.
        keyz = list(keyz)
        key = keyz[0]
    else:
        get_chromosome_from_header = options.header
        get_chromosome_from_header = \
            get_chromosome_from_header.replace('?', options.chromosome)
        for key in keyz:
            if(re.match(get_chromosome_from_header, key) is not None):
                match = key
        if(match is ''):
            raise Exception("No match possible is something wrong with the"
                            " regex specified to the program as"
                            "--header-regex")
    aaSeq = f[key]
    return(aaSeq)


def write_sample_file(options, vcf_reader):
    if(options.sample_file is not None):
        sample_file = open(options.sample_file, 'w')
        sample_header = ("ID_1 ID_2 missing father mother sex plink_pheno"
                         "\n0 0 0 D D D B\n")
        sample_file.write(sample_header)
        for sample in vcf_reader.samples:
            sample_file.write(sample + ' ' + sample + ' 0 0 0 0 -9 ' + '\n')
        sample_file.close()


def get_haps_line(options, record):
    if(record.ID is not None):
        line = (record.ID + ' ' + record.ID + ' ' + str(record.POS) +
                ' ' + str(record.REF) + ' ' + str(record.ALT[0]))
    else:
        id = options.chromosome + ":" + str(record.POS)
        line = (id + ' ' + id + ' ' + str(record.POS) + ' ' +
                str(record.REF) + ' ' + str(record.ALT[0]))

    for samples in record.samples:
        gt = samples['GT']
        # Need to skip any snps that have any missing phase data to
        # increase certainty of our results.
        # If every snp will indeed be phased
            
        if(gt != None and '|' in gt):
            gtSplit = gt.split('|')
            if(gtSplit[0] == '.'):
                gtSplit[0] = options.missing_data_code
            if(gtSplit[1] == '.'):
                gtSplit[1] = options.missing_data_code
            line = line + ' ' + gtSplit[0] + ' ' + gtSplit[1]
        else:
            line = line + ' '+options.missing_data_code +  ' '+ options.missing_data_code
    return line


def write_hap_line(options,output_line,output=None):
    """ Writes a haps file out

        Either writes a haps file to stdout or stderr.
    """
    if(output_line is not None):
        if (options.output is not None):
                output.write(output_line + "\n")
        else:
                print(output_line)

def close_files(options,output=None):
    """ Close output haps file.

    """
    if(options.output is not None):
        output.close()


def vcf_to_haps(options):
    """ Converts a VCF file to haps format

    """
    if(options.output is not None):
        output = open(options.output, 'w')
    else:
        output = None
    vcf_reader = vcf.Reader(filename=options.vcf_file)
    write_sample_file(options, vcf_reader)
    for record in vcf_reader:
        write_hap_line(options, get_haps_line(options, record), output)
    close_files(options, output)


def annotate_vcf(options):
    if(options.output is not None):
        output = open(options.output, 'w')
        if(options.output_af is not None):
            output_af = open(options.output,'a')
    else:
        output = None
    vcf_reader = vcf.Reader(filename=options.vcf_file)
    write_sample_file(options, vcf_reader)
    aaSeq = aa_seq(options)
    for record in vcf_reader:
        line = get_haps_line(options, record)
        if(line is not None):
            output_line = aa_check(aaSeq[record.POS-1], record.REF,
                                   record.ALT, options.format, line)
            if(options.output_af is not None):
                output_af.write(str(pos) + "\t" + ref + "\t" + alt + "\t" + tempSeq + "\t" + allele_freq(line.split()[5:],outputLine.split()[5:]) + "\n")
            write_hap_line(options, output_line, output)
    close_files(options, output)
    if(options.output_af is not None):
        close_files(options,output_af)


def aa_check(realAA, ref, alt, format, line):
    if(re.match('[ACTGactg]', realAA)):
        if(realAA.islower() and format == "upper"):
            return None
        else:
            if(realAA.upper() == ref.upper()):
                return line.strip()
            elif(realAA.upper() == alt.upper()):
                newLine = line.split()
                newLine[3] = alt
                newLine[4] = ref
                for i in range(5, len(newLine)):
                    if((newLine[i]) == "1"):
                        newLine[i] = '0'
                    elif((newLine[i]) == "0"):
                        newLine[i] = '1'
                    # DO nothing leave it as missing
            else:
                newLine = line.split()
                newLine[3] = realAA
                newLine[4] = ref
                for i in range(5, len(newLine)):
                        newLine[i] = '1'
            return ' '.join(newLine)
    else:
        return None

def allele_freq(ref, ances):
    ref = ' '.join(ref)
    ances = ' '.join(ances)
    p1 = ref.count('0')
    q1 = ref.count('1')
    p2 = ances.count('0')
    q2 = ances.count('1')
    maf = q1 / float(p1 + q1)
    daf = q2 /float(p2 + q2)
#    print(str(p1) +"\t" + str(q1) + "\t" +str(p2) + "\t" +str(q2))
    return str(maf) + "\t" + str(daf)


def annotate_haps(options):
    aaSeq = aa_seq(options)
    output = None
    if(options.output is not None):
        output = open(options.output, 'w')
        if(options.output_af is not None):
            output_af = open(options.output_af, 'a')
    with open(options.haps, 'r') as haps:
        for line in haps:
            lineSplit = line.split()
            pos = int(lineSplit[2])
            ref = lineSplit[3]
            alt = lineSplit[4]
            tempSeq = aaSeq[pos-1]
            outputLine = aa_check(tempSeq, ref, alt, options.format, line)
            if(outputLine is not None):
                if(options.output is not None):
                    output.write(outputLine + "\n")
                    if(options.output_af is not None):
                        output_af.write(str(pos) + "\t" + ref + "\t" + alt + "\t" + tempSeq + "\t" + allele_freq(line.split()[5:],outputLine.split()[5:]) + "\n")
                else:
                    print(outputLine)
    if(options.output is not None):
        output.close()
        if(options.output_af is not None):
            output_af.close()


def main():
    parser = OptionParser()
    parser.add_option('-i', '--haps', dest='haps',
                      help="Haplotype File (.haps)")
    parser.add_option('-a', '--aa', dest='ancestralfasta',
                      help="Ancestral Allele Fasta file")
    parser.add_option('--ref-fasta', action='store_true',
                      dest='ref_fasta',
                      help=('Use reference fasta which does not split'
                            'by chromosome'))
    parser.add_option('-c', '--chr', dest="chromosome", help="Chromosome")
    parser.add_option('-o', '--output', dest="output",
                      help="Output File (optional)")
    parser.add_option('-f', '--format', dest="format",
                      help=("Format use upper case or upper "
                            "and lower case bases"))
    parser.add_option('-v', '--phased-vcf', dest="vcf_file",
                      help="Phased VCF file (.vcf)")
    parser.add_option('-s', '--sample-file', dest="sample_file",
                      help="Output sample_file")
    parser.add_option('--header-regex', dest="header",
                      help=("To determine which chromosome to extract "
                            "is a regex with a ? for the chromosome number"))
    parser.add_option('--single-chromosome', action='store_true',
                      dest='single_chromosome')
    parser.add_option('--no-annotation', action="store_true",
                      dest="no_annotation",
                      help=("No annotation of VCF file just"
                            " convert to haps"))
    parser.add_option('--missing-code',dest='missing_data_code',
                      help='Missing code for output file')
    parser.add_option('--af',dest='output_af',help="filename for file with minor and derived allele frequencies")
    (options, args) = parser.parse_args()
    if(options.missing_data_code is None):
        options.missing_data_code = '.'
    if(options.format is None):
        options.format = 'lower'
    # Will annotate the haps file with exactly what is required
    # More options could be added later covering a wider range of file types
    # andy maybe different input ancestral alleles.
    assert options.haps is not None or options.vcf_file is not None,\
        "Haps or VCF input file required to run ancestral annotation."
    if(options.output_af is None):
        options.output_af = options.output.split('.haps')[0] + ".af"
    if(options.output_af is not None):
        f =  open(options.output_af, 'w')
        f.write("Pos\tRef\tAlt\tAnc\tMAF\tDAF\n")
        f.close()
    if(options.haps is not None):
        annotate_haps(options)
    elif(options.vcf_file is not None):
        if(options.no_annotation is None):
            annotate_vcf(options)
        else:
            vcf_to_haps(options)
    if(options.single_chromosome is None):
        options.single_chromosome = False
        assert options.header is None, \
            "Option header_regex required if the fasta file is"\
            "split by chromosome"

if __name__ == "__main__":
    main()
