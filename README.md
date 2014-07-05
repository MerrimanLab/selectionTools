selectionTools 1.0
=========================
Pipeline to take VCF through to Selection Analysis.

Software Prerequisites
---------------------

The selection pipeline was developed on a 64-bit Ubuntu 13.04 system and has been tested on 64-bit Centos and Ubuntu 13.10 installations. The pipeline should work on any 64-bit linux system and OSX (without Fay and Wu's H).

* Python >= 2.6
* Bourne-again Shell (Bash)
* Perl >=5
* R >= 3.0.0
* GNU Autotools
* GCC
* Git
* Java >= 1.7 (for beagle)

Python Dependencies

* python-setuptools
* python-numpy
* python-scipy

If you are using python < 2.7 the python package argparse will need to be installed. 
Installation
------------

After cloning or downloading and extracting running `./install.sh` in the root directory will attempt to install the pipeline and all required dependencies.

By default the pipeline executables will be added to $HOME/.local/bin, you should add this directory to your executable path.

Config File
-----------

Each run of the pipeline requires the specification of settings in a config file. A default config file is generated after installation in the selectionTools directory
named defaults.cfg. Detailed information on what the settings do and how to change them is avaliable in the pdf manual in the docs/ directory.

Single Population
-----------------

To run the selection pipeline on a single population

    selection_pipeline -c <chromosome number> -i <vcf input> --population <population name> \
    --config-file <config_file> --cores <cpu cores>

To view the other options run the help.
    
    selection_pipeline -h

Multiple Populations
--------------------

To run the selection pipeline on multiple populations.

    multipop_selection_pipeline -p <population file1> -p <population file2> \
    -i <merged input vcf> --config-file <config file> 

For more information on population files consult the PDF manual specifically section 3.3.

To view the other options run the help.

    multipop_selection_pipeline -h
