Merriman Selection Pipeline
=========================
Pipeline to take VCF through to Selection Analysis.

Software Requirements
---------------------

All executables should be added to you ${PATH} variable on linux systems.

- VCF tools required to be installed [vcftools](http://sourceforge.net/proj ects/vcftools/files/latest/download "Vcf Tools") 
- Shapeit [shapeit](http://www.shapeit.fr/ "Shapeit") executable should be called shapeit
- plink [plink](http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml "Plink")
- python [python](http://www.python.org/download/ "Python")
- R 3.0 [R](http://cran.at.r-project.org/ "R")

Once all tools are installed some package dependencies are required for python
and R. These are as follows.

- pyfasta [pyfasta](https://pypi.python.org/pypi/pyfasta/ pyfasta)
- PyVCF [PyVCF](https://github.com/jamescasbon/PyVCF PyVCF)
- rehh  [rehh] (http://cran.r-project.org/web/packages/rehh/index.html rehh)

(optional)

for multicore processing the R package multicore is required.

- multicore [multicore](http://cran.r-project.org/web/packages/multicore/ multicore)

Installation
------------

In the root directory run `./install.sh` will attempt to install as many dependencies as can be possible.



Follow-up Instructions
----------------------

If you installed the dependencies using install.sh some extra steps are required to use the selection pipeline.

- selection_pipeline executable should be added to your path automatically.




Usage
-----

Usage: selection_pipeline [options]

Options:


  -h, --help            show this help message and exit

  -v, --verbose         Print debug messages

  -q, --silent          Run Silently

  -i VCF_INPUT, --vcf=VCF_INPUT VCF input file

  -o OUTPUT_PREFIX, --out=OUTPUT_PREFIX Output file prefix

  -c CHROMOSOME, --chromosome=CHROMOSOME Chromosome

  -l LOG_FILE, --log-fire=LOG_FILE Log file for the pipeline process

  --maf=MAF             Minor allele-frequency filter

  --hwe=HWE             Hardy-Weinberg Equillibrium filter proportion

  --daf=DAF             Derived Allele Frequency filter proportion

  --remove-missing=REMOVE_MISSING Remove missing genotypes

  --config-file=CONFIG_FILE Config file

  --phased-vcf          Phased vcf file

  --population=POPULATION Population Code
  --imputation          Imputation

  --full-process        Run Entire Process

  --gzvcf               VCF input is in GZ file (optional)



Config File
-----------

The selection pipeline is configured and run using both the program arguments and a config file the specfies the more persistent features of a single program. 
