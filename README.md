# CRITCAL BUGFIX 18NOV15
https://github.com/smilefreak/selectionTools/issues/16
bug in aa_annotate.py has now been fixed **only in master (selectionTools1.1) and selectionTools1.1-dev branches or version 1.1.1+**


**selectionTools1.1 was merged onto the master branch on 12JUN15**


Original and minor updated versions of 1.0 can be found under releases

Citation
========

[![Join the chat at https://gitter.im/smilefreak/selectionTools](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/smilefreak/selectionTools?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)


Cadzow, Murray, et al. "A bioinformatics workflow for detecting signatures of selection in genomic data." Frontiers in genetics 5 (2014).

selectionTools 1.1
=========================
Pipeline to take VCF through to Selection Analysis.

the branch selectionTools1.1-dev is being used for further development 

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

__Docker__: in the docker directory there is a dockerfile and instructions to create a docker image. Also try this option if the normal installation method doesn't work

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
    -i <merged input vcf> --config-file <config file> -c <chromosome number>

For more information on population files consult the PDF manual specifically section 3.3.

To view the other options run the help.

    multipop_selection_pipeline -h
    
    
Notes
-----

With a Red Hat linux system you will need to extract the qctool scientific linux distribution located in the src directory.


