#!/bin/bash
#
# Install all the programs required for the install of the program.
#
#



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
	which vcftools
	if [ "$?" -ne "0" ]; then 
		echo "Installing VCF tools"
		tar xzf src/vcftools.tar.gz
		(cd vcftools_0.1.11/ && make)
		cp vcftools_0.1.11/bin/* bin/
		cp vcftools_0.1.11/perl/*pm lib/perl5/
		rm -Rf vcftools_0.1.11
		
	fi
	which qctool
	if [ "$?" -ne "0" ]; then
		echo "Installing QCTool"
		tar xzf src/qctool_v1.3-linux-x86_64.tgz
		mv qctool_v1.3-linux-x86_64/qctool bin/
		rm -Rf qctool_v1.3-linux-x86_64
	fi

	which shapeit
	if [ $? -ne "0" ]; then
		echo "Installing shapeit"
		tar xzf src/shapeit.v2.r727.linux.x64.tar.gz
		mv shapeit.v2.r727.linux.x64 bin/shapeit
		rm -Rf example
	fi	
	which impute2
	if [ $? -ne 0 ]; then
		echo "Installing Impute2"
		tar xzf src/impute_v2.3.0_x86_64_static.tgz
		mv impute_v2.3.0_x86_64_static/impute2 bin/
		rm -Rf impute_v2.3.0_x86_64_static/
	fi
	which tabix
	if [ $? -ne 0]; then
		echo "Installing Tabix"
		tar -xjf src/tabix.tar.bz2
		(cd tabix-0.2.6/ && make)
		cp tabix-0.2.6/bgzip bin/
		cp tabix-0.2.6/tabix bin/
		rm -Rf tabix-0.2.6
	
	fi
	echo "Install rehh"
	R CMD INSTALL src/rehh_1.11.tar.gz	
	echo "Install Selection Pipeline"
	git submodule init
	git submodule update
	(cd pyfasta && python3 setup.py install)
	(cd PyVCF && python3 setup.py install)
	python3 setup.py install
fi
