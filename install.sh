#!/bin/bash
#
# Install all the programs required for the install of the program.
#

ORIG_DIR=`pwd`

# Argument to build each function
# $1 program name
# $2 folder
change_folder(){
    cd $1
}

check_success(){
    "$@"
    status=$?
     if [ $status -ne 0 ]; then
        echo "error with $1"
        exit 1
    fi
    return $status
}

orig_dir(){
    cd $ORIG_DIR
}

mkdir -p bin
mkdir -p lib/perl5
echo "Installing Dependencies"
echo "Install Zlib"
tar xzf src/zlib-1.2.8.tar.gz
prefix_zlib=${ORIG_DIR}
change_folder zlib-1.2.8 
echo $PWD
check_success ./configure --prefix=${ORIG_DIR}
check_success make install
orig_dir
rm -Rf zlib-1.2.8

echo "Installing VCF tools"
tar xzf src/vcftools.tar.gz
LIB_VAR="-lz -L${ORIG_DIR}/lib -I${ORIG_DIR}/include"
change_folder vcftools_0.1.11
check_success make LIB="${LIB_VAR}"
orig_dir
cp vcftools_0.1.11/bin/* bin/
cp vcftools_0.1.11/perl/*pm lib/perl5/
rm -Rf vcftools_0.1.11
echo "Installing QCTool"
if [ `uname` = "Darwin" ]; then 
    tar xzf src/qctool_v1.4-osx.tgz
    mv qctool_v1.4-osx/qctool bin/
    rm -Rf qctool_v1.4-osx
else
    tar xzf src/qctool_v1.4-linux-x86_64.tgz
    mv qctool_v1.4-linux-x86_64/qctool bin/
    rm -Rf qctool_v1.4-linux-x86_64
fi

echo "Installing Shapeit"
if [ `uname` = "Darwin" ]; then
    wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r778.MacOSX.tgz
    tar xzf shapeit.v2.r778.MacOSX.tgz
    mv shapeit bin/
    rm -Rf shapeit.v2.r778.MacOSX.tgz
else
    echo `uname`
    wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r778.Ubuntu_12.04.4.static.tar.gz
    tar xzf shapeit.v2.r778.Ubuntu_12.04.4.static.tar.gz
    mv shapeit bin/
    rm -Rf shapeit.v2.r778.Ubuntu_12.04.4.static.tar.gz
fi
rm -Rf example
rm -f LICENCE
echo "Installing PLINK"
if [ `uname` = "Darwin" ]; then
	unzip src/plink-1.07-mac-intel.zip
	cp plink-1.07-mac-intel/plink bin/
	rm -Rf plink-1.07-mac-intel
else
	unzip src/plink-1.07-x86_64.zip
	cp plink-1.07-x86_64/plink bin/
	rm -Rf plink-1.07-x86_64
fi
echo "Installing Impute2"
tar xzf src/impute_v2.3.0_x86_64_static.tgz
mv impute_v2.3.0_x86_64_static/impute2 bin/
rm -Rf impute_v2.3.0_x86_64_static/
chmod 755 bin/impute2
echo "Installing Tabix"
tar -xjf src/tabix.tar.bz2
change_folder tabix-0.2.6
check_success make
orig_dir
cp tabix-0.2.6/bgzip bin/
cp tabix-0.2.6/tabix bin/
rm -Rf tabix-0.2.6
echo "Installing Variscan"
if [ `uname` = "Darwin" ]; then
    echo "Cannot install on OSX"
else
    tar -xzf src/variscan-2.0.3.tar.gz
    (cd variscan-2.0.3/src/ && rm *o)
    change_folder  variscan-2.0.3
    check_success bash autogen.sh && make
    orig_dir
    mv variscan-2.0.3/src/variscan bin/
    rm -Rf variscan-2.0.3
fi
echo "Installing Beagle"
cp src/beagle.jar bin/
echo "Installing getopt"
check_success Rscript src/R_dependencies.R 'getopt'
echo "Installing R Multicore"
check_success Rscript src/R_dependencies.R 'multicore'
echo "Installing old rehh"
check_success Rscript src/R_dependencies.R 'rehh'
echo "Install rehh"
check_success R CMD INSTALL src/rehh_1.11.tar.gz
echo "Updating submodules"
git submodule init
git submodule update
echo "Generating Default Config File"
# Because PWD contains slashes (/) need to use # as substitution
sed 's#!SELECT_PIPELINE!#'"${PWD}"'#g' src/defaults.cfg > defaults.cfg


if [[ $EUID -eq 0 ]]; then
    echo "Installing PyFasta"
    change_folder pyfasta
    check_success python setup.py install
    orig_dir
    echo "Installing PyVCF"
    change_folder PyVCF
    check_success python setup.py install
    orig_dir
    echo "Installing selection_pipeline"
    check_success python setup.py install
    else
    echo "Installing PyFasta"
    change_folder pyfasta
    check_success python setup.py install --user
    orig_dir
    change_folder PyVCF
    echo "Install PyVCF"
    check_success python setup.py install --user
    orig_dir
    echo "Installing selection_pipeline"
    check_success python setup.py install --user
fi
