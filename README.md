selectionTools 1.0
=========================
Pipeline to take VCF through to Selection Analysis.

Software Prerequisites
---------------------

The selection pipeline was developed on a 64-bit Ubuntu 13.04 system and has been tested on 64-bit Centos and Ubuntu 13.10 installations. The pipeline should work on any 64-bit linux system and Mac OSX with some manual tweaking.

* Python >= 2.6
* Bourne-again Shell (Bash)
* Perl5
* R >= 3.0.0
* GNU Autotools
* GCC
* Git

Python Dependencies

* python-setuptools
* python-numpy
* python-scipy

If you are using python < 2.7 the python package argparse will need to be installed. 
Installation
------------

After cloning or downloading and extracting running `./install.sh` in the root directory will attempt to install the pipeline and all required dependencies.

By default the pipeline executables will be added to $HOME/.local/bin.

Config File
-----------

Each run of the pipeline requires the specification of config file

Single Population
-----------------


Multiple Populations

