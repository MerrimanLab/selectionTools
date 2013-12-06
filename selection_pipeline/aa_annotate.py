#
# James Boocock
# July 2013
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
    f = Fasta(options.ancestralfasta)
    keyz = (f.keys())
    chroms = {}
    for key in keyz:
        if(options.ref_fasta is not None):
            chroms[(key.split(' '))[0]] = key
        else:
            chroms[(key.split(':'))[2]] = key
    aaSeq = f[chroms[options.chromosome]]
    return(aaSeq)


def annotate_vcf(options):
    if(options.output is not None):
        output = open(options.output, 'w')
    #TODO MAKE WORK WITH REFERENCE FASTA
    aaSeq = aa_seq(options)
    vcf_reader = vcf.Reader(filename=options.vcf_file)
    if(options.sample_file is not None):
        sample_file = open(options.sample_file, 'w')
        sample_header = ("ID_1 ID_2 missing father mother sex plink_pheno"
                         "\n 0 0 0 D D D B\n")
        sample_file.write(sample_header)
        print(options.sample_file)
        for sample in vcf_reader.samples:
            sample_file.write(sample + ' ' + sample + ' 0 0 0 0 -9 ' + '\n')
    for record in vcf_reader:
        if(record.ID is not None):
            line = (record.ID + ' ' + record.ID + ' ' + str(record.POS) +
                    ' ' + str(record.REF) + ' ' + str(record.ALT[0]))
        else:
            id = options.chromosome + ":" + record.POS
            line = (id + ' ' + id + ' ' + str(record.POS) + ' ' +
                    str(record.REF) + ' ' + str(record.ALT[0]))
        for samples in record.samples:
            gt = samples['GT']
            # Need to skip any snps that have any missing phase data to
            # increase certainty of our results.
            # If every snp will indeed be phased
            if('|' in gt):
                gtSplit = gt.split('|')
                line = line + ' ' + gtSplit[0] + ' ' + gtSplit[1]
            else:
                line = line + '- -'
                break
        if(line is not None):
            output_line = aa_check(aaSeq[record.POS], record.REF,
                                   record.ALT, options.format, line)
            if(output_line is not None):
                if (options.output is not None):
                        output.write(output_line + "\n")
                else:
                        print(output_line)
    if(options.output is not None):
        output.close()
    if(options.sample_file is not None):
        sample_file.close()


def aa_check(realAA, ref, alt, format, line):
    if(re.match('[ACTGactg]', realAA)):
        if(realAA.islower() and format == "high"):
            return None
        else:
            if(realAA == ref):
                return line.strip()
            elif(realAA == alt):
                newLine = line.split()
                newLine[3] = alt
                newLine[4] = ref
                for i in range(5, len(newLine)):
                    if((newLine[i]) == "1"):
                        newLine[i] = '0'
                    else:
                        newLine[i] = '1'
            else:
                newLine = line.split()
                newLine[3] = realAA
                newLine[4] = ref
                for i in range(5, len(newLine)):
                        newLine[i] = '1'
            return ' '.join(newLine)
    else:
        return None


def annotate_haps(options):
    aaSeq = aa_seq(options)
    output = None
    if(options.output is not None):
        output = open(options.output, 'w')
    with open(options.haps, 'r') as haps:
        for line in haps:
            lineSplit = line.split()
            pos = int(lineSplit[2])
            ref = lineSplit[3]
            alt = lineSplit[4]
            tempSeq = aaSeq[pos].decode()
            outputLine = aa_check(tempSeq, ref, alt, options.format, line)
            if(outputLine is not None):
                if(options.output is not None):
                    output.write(outputLine + "\n")
                else:
                    print(outputLine)

    if(options.output is not None):
        output.close()


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
                      help="Format us High or use Low & High")
    parser.add_option('-v', '--phased-vcf', dest="vcf_file",
                      help="Phased VCF file (.vcf)")
    parser.add_option('-s', '--sample-file', dest="sample_file",
                      help="Output sample_file")
    (options, args) = parser.parse_args()
    #print(options)
    if(options.format is None):
        options.format = 'high'
    # Will annotate the haps file with exactly what is required
    # More options could be added later covering a wider range of file types
    # andy maybe different input ancestral alleles.
    if(options.haps is not None):
        annotate_haps(options)
    elif(options.vcf_file is not None):
        annotate_vcf(options)


if __name__ == "__main__":
    main()
