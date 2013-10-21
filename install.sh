#!/bin/bash
#
# Install all the programs required for the install of the program.
#
#



if [ "$1" == "--standalone" ]; then 
	echo "Install Selection Pipeline"
	git submodule init
	git submodule update
	cd pyfasta && python setup.py install
	cd Pyfaste && python setup.py install
	python setup.py install
else
	mkdir -p bin
	echo "Installing Dependencies"
	which vcftools
	if [ "$?" -ne "0" ]; then 
		echo "Installing VCF tools"
		tar xzf src/vcftools.tar.gz
		(cd vcftools_0.1.11/ && make)
		cp vcftools_0.1.11/bin/* bin/
		rm -Rf vcftools_0.1.11
		
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
	echo "Install rehh"
	R CMD INSTALL src/rehh_1.11.tar.gz	
	
	
	
	echo "Install Selection Pipeline"
	git submodule init
	git submodule update
	cd pyfasta && python setup.py install
	cd Pyfaste && python setup.py install
	python setup.py install
fi
