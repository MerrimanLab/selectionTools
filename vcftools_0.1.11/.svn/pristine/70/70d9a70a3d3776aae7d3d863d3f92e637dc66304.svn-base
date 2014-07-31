
<div class="item">
<h1>vcftools: Usage and Options</h1>
<p>The vcftools program is intended for analysis of <b>diploid SNP data</b> in VCF format. 
The program is run from the command line, and the interface is
inspired by <a href="http://pngu.mgh.harvard.edu/%7Epurcell/plink/">PLINK</a>,
and so should be largely familiar to users of that package. Commands
take the following form:</p>
<p class="codebox">vcftools --vcf file1.vcf --chr 20 --freq
</p>
<p>The above command tells vcftools to read in the file <i>file1.vcf</i>, extract
sites on chromosome 20, and calculate the allele frequency at each
site. The resulting allele frequency estimates
are stored in the output file, <i>out.freq</i>. As in the
above example, output from vcftools is mainly sent
to output files, as opposed to being shown on the screen.
</p>
<p>
If you want to output a new VCF file, you need the --recode option. For example, the command may look something like:
</p>
<p class="codebox">vcftools --vcf file1.vcf --chr 20 --out chr20 --recode
</p>
<p>This command tells vcftools to extract chromosome 20 from <i>file1.vcf</i> and generate a new VCF file, which will be 
called <i>chr20.recode.vcf</i>. </p> 
<p>A description of each of the available options follows. <b>Note that some commands may only be 
available in the latest version of vcftools. To obtain the latest version, you should use SVN to 
checkout the latest code, as described on the <a href="http://vcftools.sourceforge.net/index.html">home page</a>.</b>
<br>
<br>
<a href="#basic">Basic&nbsp;Options</a><br>
<a href="#site_filter">Site Filter Options</a><br>
<a href="#indv_filter">Individual Filters</a><br>
<a href="#geno_filter">Genotype Filters</a><br>
<a href="#stats">Output Statistics</a><br>
<a href="#other_format">Output in Other Formats</a><br>
<a href="#misc">Miscellaneous</a><br>
<a href="#diff">File Comparison</a><br>
<a href="#development">Options still in development</a><br>
</p>
</div>
<div class="item">
<h1><a name="basic" class="Q">Basic&nbsp;Options</a></h1>
<ul>
<li>--vcf &lt;filename&gt;</li>
<br>
This option defines the VCF file to be processed. The files need to be decompressed prior to use with vcftools. 
<b>vcftools expects files in VCF format v4.0, a specification of which can be found <a href="specs.html">here</a></b>. Version
v4.1 is also supported, with some small limitations. 
<br>
<br>
</li>
<li>--gzvcf &lt;filename&gt;</li>
<br>
This option can be used in place of the --vcf option to read compressed (gzipped) VCF files directly. 
<br>
<br>
</li>
<li>--bcf &lt;filename&gt;</li>
<br>
This option can be used in place of the --vcf option to read BCF files directly. 
You do not need to specify if this file is compressed with BGZF encoding.
Due to differences in current BCF formatting, you may need to also include the --gatk tag to signify that your BCF file came from that program.
<br>
<br>
</li>
<li>--out &lt;prefix&gt;<br>
<br>
This option defines the output filename prefix for all files generated
by vcftools. For example, if &lt;prefix&gt; is set to&nbsp;<span style="font-style: italic;">output_filename</span>,
then all output files will be of the form <span style="font-style: italic;">output_filename.*** </span>.
If this option is omitted, all output files will have the prefix 'out.'.<br>
<br>
</li>
</ul>
</div>
<div class="item">
<h1><a name="site_filter" class="Q">Site Filter Options</a></h1>
<ul>
<li>--chr &lt;chromosome&gt;<br>
<li>--not-chr &lt;chromosome&gt;<br>
<br>
Include or exclude chromosomes with identifiers matching
&lt;chromosome&gt;<br>. These options can be used more than once to include or exclude multiple chromosomes.
<br>
</li>
<li>--from-bp &lt;integer&gt;</li>
<li>--to-bp &lt;integer&gt;<br>
<br>
These options define the&nbsp;physical range of sites will be
processed. Sites outside of this range will be excluded. These options
can only be used in conjuction with --chr.<br>
<br>
</li>
<li>--snp &lt;string&gt;<br>
<br>
Include SNP(s) with matching ID (e.g. a dbSNP rsID). This command can be used multiple
times in order to include more than one SNP.<br>
<br>
</li>
<li>--snps &lt;filename&gt;<br>
<br>
Include a list of SNPs given in a file. The file should contain a list
of SNP IDs (e.g. dbSNP rsIDs), with one ID per line.<br>
<br>
</li>
<li>--exclude &lt;filename&gt;<br>
<br>
Exclude a list of SNPs given in a file. The file should contain a list
of SNP IDs, with one ID per line.<br>
<br>
</li>
<li>--positions &lt;filename&gt;</li>
<li>--exclude-positions &lt;filename&gt;<br>
<br>
Include/exclude a set of sites on the basis of a 
list of positions in a file. Each line of the input file should contain a (tab-separated) chromosome and position.
The file can have comment lines that start with a '#'. <br>
<br>
</li>
<li>--keep-only-indels</li>
<li>--remove-indels<br>
<br>
Include or exclude sites that contain an indel. For this option 'indel' means any variant that 
alters the length of the REF allele.<br>
<br>
</li>
<li>--bed &lt;filename&gt;<br>
<li>--exclude-bed &lt;filename&gt;<br>
<br>
Include or exclude a set of sites on the basis of a 
<a href="http://genome.ucsc.edu/FAQ/FAQformat#format1">BED</a> file. Only 
the first three columns (chrom, chromStart and chromEnd) are required. 
The BED file should have a header line. <br>
<br>
</li>
<li>--remove-filtered-all</li>
<li>--remove-filtered &lt;string&gt;<br>
<li>--keep-filtered &lt;string&gt;<br>
<br>
These options are used to filter sites on the basis of their FILTER flag. 
The first option removes all sites with a FILTER flag. 
The second option can be used to exclude sites with a specific filter flag. 
The third option can be used to select sites on the basis of specific filter flags. The second and third options can be used multiple times to specify multiple FILTERs. 
The --keep-filtered option is applied before the --remove-filtered option.
NOTE: The PASS tag does not count as a filter, it is just the absence of a filter. If you would like to keep only the sites that pass all filters use the --remove-filtered-all option.
<br>
<br>
</li>
<li>--remove-INFO &lt;string&gt;<br>
<li>--keep-INFO &lt;string&gt;<br>
<br>
These options are used to filter sites on the basis of INFO field flags. 
The first option can be used to exclude sites with a specific INFO flag. 
The second option can be used to select sites on the basis of specific INFO flags. These options can be used multiple times to specify multiple INFO flags. 
The --keep-INFO option is applied before the --remove-INFO option. Note that <b>only INFO flags</b> can currently be used as filters 
(i.e. there is currently no support for filtering by INFO field values). <br>
<br>
</li>
<li>--minQ &lt;float&gt;<br>
<br>
Include only sites with Quality above this threshold.<br>
<br>
</li>
<li>--min-meanDP &lt;float&gt;</li>
<li>--max-meanDP &lt;float&gt;<br>
<br>
Include sites with mean Depth within the thresholds defined by these
options.<br>
<br>
</li>
<li>--maf &lt;float&gt;</li>
<li>--max-maf &lt;float&gt;<br>
<br>
Include only sites with Minor Allele Frequency within the specified
range.<br>
<br>
</li>
<li>--non-ref-af &lt;float&gt;</li>
<li>--max-non-ref-af &lt;float&gt;<br>
<br>
Include only sites with all Non-Reference Allele Frequencies within the specified
range.<br>
<br>
</li>
<li>--mac &lt;int&gt;</li>
<li>--max-mac &lt;int&gt;<br>
<br>
Include only sites with Minor Allele Count within the specified
range.<br>
<br>
</li>
<li>--non-ref-ac &lt;float&gt;</li>
<li>--max-non-ref-ac &lt;float&gt;<br>
<br>
Include only sites with all Non-Reference Allele Counts within the specified
range.<br>
<br>
</li>
<li>--hwe &lt;float&gt;<br>
<br>
Assesses sites for Hardy-Weinberg Equilibrium using an exact test, as
defined by Wigginton, Cutler and Abecasis (2005). Sites with a p-value
below the threshold defined by this option are taken to be out of HWE,
and therefore excluded.<br>
<br>
</li>
<li>--geno &lt;float&gt;<br>
<br>
Exclude sites on the basis of the proportion of missing data (defined
to be between 0 and 1, where 1 indicates no missing data allowed).<br>
<br>
</li>
<li>--max-missing-count &lt;int&gt;<br>
<br>
Exclude sites with more than this number of missing chromosomes.<br>
<br>
</li>
<li>--min-alleles &lt;int&gt;</li>
<li>--max-alleles &lt;int&gt;<br>
<br>
Include only sites with a number of alleles within the specified range.
For example, to include only bi-allelic sites, one could use:<br>
<br>
<p class="codebox">vcftools --vcf file1.vcf --min-alleles 2 --max-alleles 2
</p>
<br><br>
</li>
<li>--thin &lt;int&gt;<br>
<br>
Thin sites so that no two sites are within the specified distance.<br>
<br>
</li>

<li>--mask &lt;filename&gt;<br>
<li>--invert-mask &lt;filename&gt;<br>
<li>--mask-min &lt;int&gt;<br>
<br>
Include sites on the basis of a FASTA-like file. The provided file contains a sequence of integer digits (between 0 and 9) for each position on a chromosome that specify if a site at that position 
should be filtered or not. An example mask file would look like:<br>
<p class="codebox">&gt;1<br>0000011111222... </p>
In this example, sites in the VCF file located within the first 5 bases of the start of chromosome 1 would be kept, whereas sites at position 6 onwards would be filtered out.
The threshold integer that determines if sites are filtered or not is set using the --mask-min option, which defaults to 0. 
The chromosomes contained in the mask file must be sorted in the same order as the VCF file. 
The --mask option is used to specify the mask file to be used, whereas the --invert-mask option can be used to specify a mask file that will be inverted before being applied.
<br>
</li>

</ul>
</div>
<div class="item">
<h1><a name="indv_filter" class="Q">Individual Filters</a></h1>
<ul>
<li>--indv &lt;string&gt;<br>
<br>
Specify an individual to be kept in the analysis. This option can be used multiple times
to specify multiple individuals.
<br><br>
</li>
<li>--keep &lt;filename&gt;<br>
<br>
Provide a file containing a list of individuals to include in
subsequent analysis. Each individual ID (as defined in the VCF
headerline) should be included on a separate line.<br>
<br>
</li>
<li>--remove-indv &lt;string&gt;<br>
<br>
Specify an individual to be removed from the analysis. This option can be used multiple times
to specify multiple individuals. If the --indv option is also specified, then the --indv option is
executed before the --remove-indv option.
<br><br>
</li>
<li>--remove &lt;filename&gt;<br>
<br>
Provide a file containing a list of individuals to exclude in
subsequent analysis. Each individual ID (as defined in the VCF
headerline) should be included on a separate line. If both the --keep
and the --remove options are used, then the --keep option is execute before
the --remove option.<br>
<br>
</li>
<li>--min-indv-meanDP &lt;float&gt;</li>
<li>--max-indv-mean-DP &lt;float&gt;<br>
<br>
Calculate the mean coverage on a per-individual basis. Only individuals
with coverage within the range specified by these options are included
in subsequent analyses.<br>
<br>
</li>
<li>--mind &lt;float&gt;<br>
<br>
Specify the minimum call rate threshold for each individual. <br>
<br>
</li>
<li>--phased<br>
<br>
First excludes all individuals having all genotypes unphased, and <span style="text-decoration: underline;">subsequently excludes
all sites with unphased genotypes</span>. The remaining
data&nbsp;therefore consists of&nbsp;phased data only.
<br><br>
</li>
<li>--max-indv &lt;int&gt;<br>
<br>
Randomly thins individuals so that only the specified number are retained.
</li>
</ul>
</div>
<div class="item">
<h1><a name="geno_filter" class="Q">Genotype Filters</a></h1>
<ul>
<li>--remove-filtered-geno-all</li>
<li>--remove-filtered-geno &lt;string&gt;<br>
<br>
The first option removes all genotypes with a FILTER flag not equal to PASS or '.'. The second option can be used to exclude genotypes with a specific filter flag.<br>
<br>
</li>
<li>--minGQ &lt;float&gt;<br>
<br>
Exclude all genotypes with a quality below the threshold specified by
this option (GQ).<br>
<br>
</li>
<li>--minDP &lt;float&gt;<br>
<li>--maxDP &lt;float&gt;<br>
<br>
Exclude all genotypes with a sequencing depth (DP) outside the range defined by these options.</li>
</ul>
</div>
<div class="item">
<h1><a name="stats" class="Q">Output Statistics</a></h1>
<ul>
<li>--freq</li>
<li>--counts</li>
<li>--freq2</li>
<li>--counts2<br>
<br>
Output per-site frequency information. The --freq outputs the allele frequency in a 
file with the suffix '.frq'. The --counts option outputs a similar file with the suffix '.frq.count', 
that contains the raw allele counts at each site.<br>
The --freq2 and --count2 options are used to suppress allele information in the output file. In this case, 
the order of the freqs/counts depends on the numbering in the VCF file.<br>
With any of these options, adding the --derived command will re-order the output file columns so that the ancestral allele 
appears first. The --derived option relies on the ancestral allele being specified in the VCF file using the AA tag in the INFO field.<br>
<br>
</li>
<li>--depth<br>
<br>
Generates a file containing the mean depth per individual. This file
has the suffix '.idepth'.<br>
<br>
</li>
<li>--site-depth</li>
<li>--site-mean-depth<br>
<br>
Generates a file containing the depth per site. The --site-depth option
outputs the depth for each site summed across individuals. This file
has the suffix '.ldepth'. Likewise, the --site-mean-depth outputs the
mean depth for each site, and the output file has the suffix
'.ldepth.mean'.<br>
<br>
</li>
<li>--geno-depth<br>
<br>
Generates
a (possibly very large) file containing the depth for each genotype in
the VCF file. Missing entries are given the value -1. The file has the suffix '.gdepth'.<br>
<br>
</li>
<li>--site-quality<br>
<br>
Generates a file containing the per-site SNP quality, as found in the QUAL column of the VCF file. This file
has the suffix '.lqual'.<br>
<br>
</li>
<li>--het<br>
<br>
Calculates
a measure of heterozygosity on a per-individual basis. Specfically, the
inbreeding coefficient, F, is estimated&nbsp;for each individual
using
a method of moments. The resulting file has the suffix '.het'.<br>
<br>
</li>
<li>--hardy<br>
<br>
Reports a p-value for each&nbsp;site from a&nbsp;Hardy-Weinberg
Equilibrium&nbsp;test (as
defined by Wigginton, Cutler and Abecasis (2005)). The resulting file
(with suffix '.hwe') also contains the Observed numbers of Homozygotes
and Heterozygotes and the corresponding&nbsp;Expected numbers under
HWE.&nbsp;<br>
<br>
</li>
<li>--missing<br>
<br>
Generates
two files reporting the missingness on a per-individual and per-site
basis. The two files have suffixes '.imiss' and '.lmiss' respectively.<br>
<br>
</li>
<li>--hap-r2</li>
<li>--geno-r2</li>
<li>--geno-chisq</li>
<li>--ld-window &lt;int&gt;</li>
<li>--ld-window-bp &lt;int&gt;</li>
<!--<li>--ld-window &lt;int&gt;</li>
<li>--ld-window-bp &lt;int&gt;</li>-->
<li>--min-r2 &lt;float&gt;<br>
<br>
These options are used to report Linkage Disequilibrium (LD) statistics
as summarised by the r<sup>2</sup>, D, D', or Chi-squared statistics,. The --hap-r2
option informs vcftools to output a file reporting the r<sup>2</sup>, D, and D'
statistics using phased haplotypes. These are the traditional measures of LD often reported in the population genetics literature. 
If phased haplotypes are unavailable then the --geno-r2 option may be used, which calculates the squared correlation coefficient
between genotypes encoded as 0, 1 and 2 to represent the number of non-reference alleles in each individual. This is the same as the LD measure reported by PLINK. The D and D' statistics are only available for phased genotypes.
If your data contains sites with more than two alleles, then the --geno-chisq option can be used to test for genotype independence via the chi-squared statistic. 
The haplotype version outputs a file with the suffix '.hap.ld', whereas the genotype version outputs a file with the suffix '.geno.ld'. The haplotype version implies the option --phased.
<br><br>The --ld-window option defines the maximum number of SNPs between the SNPs being tested for LD.
Likewise, the --ld-window-bp option can be used to define the maximum 
physical separation (in base-pairs) of SNPs included in the LD calculation. Finally,
the --min-r2 sets a minimum value for&nbsp;r<sup>2</sup>
below which the LD statistic is not reported.<br>
<br>
</li>
<li>--SNPdensity &lt;int&gt;<br>
<br>
Calculates the number and density of SNPs in bins of size defined by
this option. The resulting output file has the suffix '.snpden'.<br>
<br>
</li>
<li>--TsTv &lt;int&gt;<br>
<br>
Calculates the Transition / Transversion ratio in bins of size defined by
this option. Only uses bi-allelic SNPs. The resulting output file has the suffix '.TsTv'. A summary is also supplied in a file with the suffix '.TsTv.summary'.<br>
<br>
</li>
<li>--TsTv-by-count<br>
<br>
Calculates the Transition / Transversion ratio as a function of alternative allele count. Only uses bi-allelic SNPs. 
The resulting output file has the suffix '.TsTv.count'. <br>
<br>
</li>
<li>--TsTv-by-qual<br>
<br>
Calculates the Transition / Transversion ratio as a function of SNP quality threshold. Only uses bi-allelic SNPs. 
The resulting output file has the suffix '.TsTv.qual'. <br>
<br>
</li>
<li>--FILTER-summary<br>
<br>
Generates a summary of the number of SNPs and Ts/Tv ratio for each FILTER category. The output file 
has the suffix '.FILTER.summary. <br>
<br>
</li>
<li>--filtered-sites<br>
<br>
Creates two files listing sites that have been kept or removed after filtering. The first file, with
suffix '.kept.sites', lists sites kept by vcftools after filters have been applied. The second file,
with the suffix '.removed.sites', list sites removed by the applied filters.<br>
<br>
</li>
<li>--singletons<br>
<br>
This option will generate a file detailing the location of singletons,
and the&nbsp;individual they occur in. The file reports both true
singletons, and private doubletons (i.e. SNPs where the minor allele
only occurs in a single individual and that individual is homozygotic
for that allele). The output file has the suffix '.singletons'.<br>
<br>
</li>
<li>--site-pi</li>
<li>--window-pi &lt;int&gt;&nbsp;(--window-pi-step &lt;int&gt;)<br>
<br>
These options are used to estimate levels of <a href="http://en.wikipedia.org/wiki/Nucleotide_diversity">nucleotide diversity</a>. The first
option does this on a per-site basis, and the output file has the suffix '.sites.pi'. The second option calculates the nucleotide diversity
in windows, with the window size defined in the option argument, and an optional window step. <br>
<br>
</li>
<li>--hist-indel-len<br>
<br>
This option will generate a histogram file of the length of all indels (including SNPs). It shows both the count and the percentage of all
indels for indel lengths that occur at least once in the input file. SNPs are considered indels with length zero. The output file has the 
suffix '.indel.hist'.<br>
<br>
</li>
<li>--TajimaD &lt;int&gt; <br>
<br>
Output Tajima's D statistic in bins of size &lt;int&gt;.
<br>
<br>
</li>
<!--
<li>--hapmap-fst-file &lt;filename&gt;<br>
<li>--hapmap-fst-gzfile &lt;filename&gt;<br>
<br>
Calculate F<sub>ST</sub> for a pair of VCF files, with the second file being specified by this option. 
F<sub>ST</sub> is currently calculated using the formula described in the supplementary material of the <a href="http://www.ncbi.nlm.nih.gov/pubmed/16255080">Phase I HapMap paper</a>.
Currently, only pairwise F<sub>ST</sub> calculations are supported, although this will likely change in the future. 
The --gzfst option can be used to read compressed VCF files. 
<u>This option is likely to be depreciated in future releases</u>, and it is recommended that you use the --hapmap-fst-pop or --weir-fst-pop options instead.
<br>
<br>
</li>
//-->
<li>--hapmap-fst-pop &lt;filename&gt;<br>
<li>--weir-fst-pop &lt;filename&gt;<br>
<li>--fst-window-size &lt;int&gt;<br>
<li>--fst-window-step &lt;int&gt;<br>
<br>
These options are used to calculate two estimates of F<sub>ST</sub>. The --hapmap-fst-pop option is used to calculate the estimator described in the
<a href="http://www.ncbi.nlm.nih.gov/pubmed/16255080">Phase I HapMap paper</a>. The second (preferred) estimator is from <a href="http://www.jstor.org/stable/2408641">Weir and Cockerham's 1984 paper</a>. These options are used to provide a file which lists the individuals in a given population (one individual per line), and is used multiple times to specific multiple populations. These options can be used to calculate F<sub>ST</sub> for more than two populations. <br><br>
The --fst-window-size and --fst-window-step options are optional and can be used to calculate F<sub>ST</sub> in windows. If they are not specified, F<sub>ST</sub> is calculated on a per-site basis.
<br>
<br>
</li>
</ul>
</div>
<div class="item">
<h1><a name="other_format" class="Q">Output in Other Formats</a></h1>
<ul>
<li>--012<br>
<br>
This option outputs the genotypes as a large matrix. Three files are
produced. The first, with suffix '.012', contains the genotypes of each
individual on a separate line. Genotypes are represented as 0, 1 and 2,
where the number represent that number of non-reference alleles. Missing 
genotypes are represented by -1. The second file, with suffix '.012.indv' 
details the individuals included in the main file. The third file, with 
suffix '.012.pos' details the site locations included in the main file.<br>
<br>
</li>
<li>--IMPUTE<br>
<br>
This option outputs phased haplotypes in <a href="http://mathgen.stats.ox.ac.uk/impute/impute.html">IMPUTE</a>
reference-panel format. As IMPUTE requires phased data, using this
option also implies --phased. Unphased individuals and genotypes are
therefore excluded. Only bi-allelic sites are included in the output. Using this option generates three files. The IMPUTE
haplotype file has the suffix '.impute.hap', and the IMPUTE legend file
has the suffix '.impute.hap.legend'. The third file, with suffix
'.impute.hap.indv', details the individuals included in the haplotype
file, although this file is not needed by IMPUTE.<br>
<br>
</li>
<li>--ldhat<br>
<li>--ldhat-geno<br>
<br>
These options output data in <a href="http://www.stats.ox.ac.uk/%7Emcvean/LDhat/">LDhat</a>
format. Use of these options &nbsp;also require the --chr option to
by used. The --ldhat option outputs phased data only, and therefore also implies --phased, leading to unphased
individuals and genotypes being excluded. Alternatively, the --ldhat-geno option treats all of the data as unphased, and therefore
outputs LDhat files in genotype/unphased format. In either case, two files are
generated with the suffixes '.ldhat.sites'
and '.ldhat.locs', which correspond to the
LDhat 'sites' and 'locs' input files respectively.<br>
<br>
</li>
<li>--BEAGLE-GL<br>
<li>--BEAGLE-PL<br>
<br>
These options output genotype likelihood information for input into the 
<a href="http://faculty.washington.edu/browning/beagle/beagle.html">BEAGLE</a> program. 
These options require the VCF file to contain FORMAT fields with GL or PL tags, which can generally be output by 
SNP callers such as <a href="http://www.broadinstitute.org/gsa/wiki/index.php/The_Genome_Analysis_Toolkit">the GATK</a>. 
Use of this option requires a chromosome to be specified via the --chr option. The resulting output 
file (with the suffix '.BEAGLE.GL' or '.BEAGLE.PL') contains genotype likelihoods for biallelic sites, and is 
suitable for input into BEAGLE via the 'like=' argument.<br>
<br>
</li>
<li>--plink<br>
<br>
This option outputs the genotype data in <a href="http://pngu.mgh.harvard.edu/%7Epurcell/plink/">PLINK</a>
PED format. Two files are generated, with suffixes '.ped' and '.map'. Note that only bi-allelic loci will be output.
Further details of these files can be found in the <a href="http://pngu.mgh.harvard.edu/%7Epurcell/plink/data.shtml#ped">PLINK
documentation</a>. <br>
<br>Note: This option can be very slow on large datasets. Using the --chr option to divide up the dataset is advised, or alternatively use the --plink-tped option below.<br>
<br>
</li>
<li>--plink-tped<br>
<br>
The --plink option above can be extremely slow on large datasets. An alternative that might be considerably quicker is 
to output in the <a href="http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr">PLINK transposed format</a>. 
This can be achieved using the --plink-tped option, which produces two files with suffixes '.tped' and '.tfam'.
<br>
<br>
</li>
<li>--recode<br>
<br>
The --recode option is used to generate a VCF file from the input VCF or BCF file
having applied the options specified by the user. The output file has
the suffix '.recode.vcf'. 
<br><br>
By default, the INFO fields are removed from the output file, as the INFO values may be invalidated 
by the recoding (e.g. the total depth may need to be recalculated if individuals are removed). 
This default functionality can be overridden by using the --recode-INFO &lt;string&gt; option, where &lt;string&gt;
defines the INFO key to keep in the output file. The --recode-INFO flag can be used multiple times. 
Alternatively, the option --recode-INFO-all can be used to retain all INFO fields.
<br><br>
</li>
<li>--recode-bcf<br>
<br>
The --recode-bcf option is used to generate a BCF file from the input VCF or BCF file
having applied the options specified by the user. The output file has
the suffix '.recode.bcf'. 
<br><br>
Similar to the --recode option, the INFO fields are removed from the output file but can be restored in the same manner.
<br><br>
</li>
<li>--recode-to-stream<br>
<br>
The --recode-to-stream option works in the same manner as the --recode option, but 
the resulting VCF file printed directly to the screen, rather than to file. Logging information normally
displayed on the screen is supressed (but still written to the log file). The advantage 
of this option is that the resulting VCF file can be piped to other programs for processing.
For example, to gzip the VCF file, one could use:<br><br>
<p class="codebox">vcftools --vcf file1.vcf --recode-to-stream | gzip -c - > output_file.vcf.gz</p>
<br><br>
</li>
<li>--recode-bcf-to-stream<br>
<br>
The --recode-bcf-to-stream option works in the same matter as the --recode-to-stream option,
printing the output BCF file to the screen. By default, all output BCF files are compressed in BGZF. 
This stream can be piped through gzip or bgzip to decompress.  
<br><br>
</li>
</ul>
</div>
<div class="item">
<h1><a name="misc" class="Q">Miscellaneous</a></h1>
<ul>
<li>--extract-FORMAT-info &lt;string&gt;<br>
<br>
Extract information from the genotype fields in the VCF file relating to a specfied FORMAT identifier. For example,
using the option '--extract-FORMAT-info GT' would extract the all of the GT (i.e. Genotype) entries. 
The resulting output file has the suffix '.&lt;FORMAT_ID&gt;.FORMAT'.
<br>
<br>
<li>--get-INFO &lt;string&gt;<br>
<br>
This option is used to extract information from the INFO field in the
VCF file. The &lt;string&gt; argument specifies the INFO tag to
be extracted, and the option can be used multiple times in order to
extract multiple INFO entries. The resulting file, with suffix '.INFO',
contains the required INFO information in a tab-separated table. For
example, to extract the NS and DB flags, one would use the command:
<br>
<br>
<p class="codebox">vcftools --vcf file1.vcf --get-INFO NS --get-INFO DB </p>
<br>
<li>--force-index-write<br>
<br>
When vcftools first reads a VCF file, it creates an index of the file for faster access next time. 
In some cases, the index file can get out of sync with the VCF file. This option can be used to force 
vcftools to regenerate the index file.
<br>
<br>
</li>
</ul>
</div>

<div class="item">
<h1><a name="diff" class="Q">VCF File Comparison Options</a></h1>
<ul>
<li>--diff &lt;filename&gt;<br>
<li>--gzdiff &lt;filename&gt;<br>
<li>--diff-bcf &lt;filename&gt;<br>
<br>
Select a VCF file for comparison with the file specified by the --vcf option. Outputs two files 
describing the sites and individuals common / unique to each file. These files have the suffixes 
'.diff.sites_in_files' and '.diff.indv_in_files' respectively. The --gzdiff version can be used to 
read compressed VCF files and the --diff-bcf version can be used to read BCF files.
<br><br>
</li>
<li>--diff-site-discordance<br>
<br>
Used in conjuction with the --diff option to calculate discordance on a site by site basis. The 
resulting output file has the suffix '.diff.sites'.
<br><br>
</li>
<li>--diff-indv-discordance<br>
<li>--diff-indv-map &lt;filename&gt;<br>
<br>
Used in conjuction with the --diff option to calculate discordance on a per-individual basis. The 
resulting output file has the suffix '.diff.indv'. The --diff-indv-map option can be used to specify
a mapping of individual IDs in the second file to those in the first file.
<br><br>
</li>
<li>--diff-discordance-matrix<br>
<br>
Used in conjuction with the --diff option to calculate a discordance matrix. This option only works 
with bi-allelic loci with matching alleles that are present in both files. The resulting output file 
has the suffix '.diff.discordance.matrix'.
<br><br>
</li>
<li>--diff-switch-error<br>
<br>
Used in conjuction with the --diff option to calculate phasing errors 
(specifically '<a href="http://www.cell.com/AJHG/retrieve/pii/S0002929707619788">switch errors</a>'). 
This option generates two output files describing switch errors found between sites, and the average switch 
error per individual. These two files have the suffixes '.diff.switch' and '.diff.indv.switch' respectively. 
<br><br>
</li>

</ul>
</div>

<div class="item">
<h1><a name="development" class="Q">Options still in development</a></h1>
<p>The following options are yet to be finalised, are likely to contain bugs, and are likely to change in the future.<br>
</p>
<ul>
<!--
<li>--pca <br>
<br>
Perform Principal Component Analysis.
<br>
<br>
</li>
//-->
<li>--LROH <br>
<br>
Identify Long Runs of Homozygosity. This option is currently <b>very</b> slow.
<br>
<br>
</li>
<li>--relatedness <br>
<br>
Output Individual Relatedness Statistics.
<br>
<br>
</li>
</ul>
</div>

