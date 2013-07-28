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



