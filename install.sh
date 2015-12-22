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
wget -O vcftools-0.1.14.tar.gz 'https://github.com/vcftools/vcftools/releases/download/v0.1.14/vcftools-0.1.14.tar.gz'
tar xzf vcftools-0.1.14.tar.gz
#tar xzf src/vcftools.tar.gz
LIB_VAR="-lz -L${ORIG_DIR}/lib -I${ORIG_DIR}/include"
change_folder vcftools-0.1.14
./configure
check_success make LIB="${LIB_VAR}"
orig_dir
cp vcftools_0.1.14/bin/* bin/
cp vcftools_0.1.14/perl/*pm lib/perl5/
rm -rf vcftools-0.1.14
rm vcftools-0.1.14.tar.gz


echo "Installing VCFlib"
unzip src/vcflib.zip
cp vcflib-master/bin/vcfsnps bin/
chmod 755 bin/vcfsnps
rm -Rf vcflib-master

#echo "Installing selscan"
#unzip src/selscan_10_Aug_2015.zip
#cp src/selscan-master/bin/linux/selscan bin/
#chmod 755 bin/selscan
#rm -Rf selscan-master 


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
    wget -O shapeit.v2.r837.MacOSX.tgz 'https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.MacOSX.tgz'
    tar xzf shapeit.v2.r837.MacOSX.tgz
    rm shapeit.v2.r837.MacOSX.tgz
    mv shapeit bin/
    rm -Rf shapeit.v2.r837.MacOSX.tgz
else
    echo `uname`
    wget -O shapeit.v2.r837.GLIBCv2.20.Linux.static.tgz 'https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.GLIBCv2.20.Linux.static.tgz'
    #http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/old_versions/shapeit.v2.r790.Ubuntu_12.04.4.static.tar.gz
    #wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r790.Ubuntu_12.04.4.static.tar.gz
    tar xzf shapeit.v2.r837.GLIBCv2.20.Linux.static.tgz
    rm shapeit.v2.r837.GLIBCv2.20.Linux.static.tgz
    mv shapeit bin/
    rm -Rf shapeit.v2.r837.GLIBCv2.20.Linux.static
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
if [ `uname` = "Darwin" ]; then
    tar xzf src/impute_v2.3.1_MacOSX_Intel.tgz
    mv impute_v2.3.1_MacOSX_Intel/impute2 bin/
    rm -Rf impute_v2.3.1_MacOSX_Intel
else
    tar xzf src/impute_v2.3.1_x86_64_static.tgz
    mv impute_v2.3.1_x86_64_static/impute2 bin/
    rm -rf impute_v2.3.1_x86_64_static/
fi
chmod 755 bin/impute2
echo "Installing Tabix"
wget -O htslib-1.2.1.tar.bz2 'https://github.com/samtools/htslib/releases/download/1.2.1/htslib-1.2.1.tar.bz2'
tar -xjf htslib-1.2.1.tar.bz2
change_folder htslib-1.2.1
#tar -xjf src/tabix.tar.bz2
#change_folder tabix-0.2.6
./configure
check_success make LIBPATH="${LIB_VAR}"
orig_dir
cp htslib-1.2.1/bgzip bin/
cp htslib-1.2.1/tabix bin/
rm htslib-1.2.1.tar.bz2
rm -rf htslib-1.2.1
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
echo "Installing old rehh"
check_success Rscript src/R_dependencies.R 'rehh'
echo "Install rehh"
check_success R CMD INSTALL src/rehh_1.11.tar.gz
echo "Install selscan"
if [ `uname` = "Darwin" ]; then
    wget -O selscan-osx-1.1.0a.tar.gz 'https://github.com/szpiech/selscan/releases/download/1.1.0a/selscan-osx-1.1.0a.tar.gz'
    tar -xzf selscan-osx-1.1.0a.tar.gz
    cp selscan-osx-1.1.0a/norm bin/
    cp selscan-osx-1.1.0a/selscan bin/
    rm -rf selscan-osx-1.1.0a/
    rm -rf selscan-osx-1.1.0a.tar.gz
else
    wget -O selscan-linux-1.1.0a.tar.gz 'https://github.com/szpiech/selscan/releases/download/1.1.0a/selscan-linux-1.1.0a.tar.gz'
    tar -xzf selscan-linux-1.1.0a.tar.gz
    cp selscan-linux-1.1.0a/norm bin/
    cp selscan-linux-1.1.0a/selscan bin/
    rm -rf selscan-linux-1.1.0a/
    rm selscan-linux-1.1.0a.tar.gz
    rm -rf ._selscan-linux-1.1.0a

fi    
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
