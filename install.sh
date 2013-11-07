#!/bin/bash
#
# Install all the programs required for the install of the program.
#

PWD=`pwd`

if [ "$1" == "--standalone" ]; then 
	echo "Install Selection Pipeline"
	git submodule init
	git submodule update
 	(cd pyfasta && python3 setup.py install)
	(cd PyVCF && python3 setup.py install)
	python3 setup.py install
else
	mkdir -p bin
	mkdir -p lib/perl5
	echo "Installing Dependencies"
	echo "Install Zlib"
	tar xzf src/zlib-1.2.8.tar.gz
	#(cd zlib-1.2.8/ && ./configure --prefix ${PWD} && make install)	
	rm -Rf zlib-1.2.8
		
	echo "Installing VCF tools"
	tar xzf src/vcftools.tar.gz
	#(cd vcftools_0.1.11/ && make)
	cp vcftools_0.1.11/bin/* bin/
	cp vcftools_0.1.11/perl/*pm lib/perl5/
	rm -Rf vcftools_0.1.11
	echo "Installing QCTool"
	tar xzf src/qctool_v1.3-linux-x86_64.tgz
	mv qctool_v1.3-linux-x86_64/qctool bin/
	rm -Rf qctool_v1.3-linux-x86_64
	echo "Installing PLINK"
	unzip src/plink-1.07-x86_64.zip
	cp plink-1.07-x86_64/plink bin/
	rm -Rf plink-1.07-x86_64
	echo "Installing shapeit"
	tar xzf src/shapeit.v2.r727.linux.x64.tar.gz
	mv shapeit.v2.r727.linux.x64 bin/shapeit
	rm -Rf example
	echo "Installing Impute2"
	tar xzf src/impute_v2.3.0_x86_64_static.tgz
	mv impute_v2.3.0_x86_64_static/impute2 bin/
	rm -Rf impute_v2.3.0_x86_64_static/
	chmod 755 bin/impute2
	echo "Installing Tabix"
	tar -xjf src/tabix.tar.bz2
	#(cd tabix-0.2.6/ && make)
	cp tabix-0.2.6/bgzip bin/
	cp tabix-0.2.6/tabix bin/
	rm -Rf tabix-0.2.6
	echo "Installing Variscan"
 	tar -xzf src/variscan-2.0.3.tar.gz
	mv variscan-2.0.3/bin/Linux-i386/variscan bin/
	rm -Rf variscan-2.0.3
	echo "Install rehh"
	R CMD INSTALL src/rehh_1.11.tar.gz
	echo "Installing R Multicore"
	R CMD INSTALL src/multicore_0.1-7.tar.gz	
	echo "Install Selection Pipeline"
	git submodule init
	git submodule update
	echo "Generating Default Config File"
  # Because PWD contains slashes (/) need to use # as substitution
	sed 's#!SELECT_PIPELINE!#'"${PWD}"'#g' corescripts/defaults.cfg > defaults.cfg
	if [[ $EUID -eq 0 ]]; then
	(cd pyfasta && python3 setup.py install)
	(cd PyVCF && python3 setup.py install)
	python3 setup.py install
	else
	(cd pyfasta && python3 setup.py install --user)
	(cd PyVCF && python3 setup.py install --user)
	python3 setup.py install --user
	fi
fi
