#
# 1) Compiling
#   Type make in this directory
#
# 2) Installation
#   Edit the BINDIR and MODDIR below as necessary or pass the PREFIX variable
#   to the make command. When not set, the programs will be placed in "bin" 
#   and "lib" subdirectories in this directory.
#       PREFIX="/install/to/path/prefix" make install
#
#   Add the MODDIR to your PERL5LIB environment variable:
#       export PERL5LIB=${PREFIX}/lib:${PERL5LIB}
#


ifndef PREFIX
    export PREFIX = $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
endif
export BINDIR = ${PREFIX}/bin
export MODDIR = ${PREFIX}/lib/perl5/site_perl

DIRS = cpp perl
install:
	    @mkdir -p $(BINDIR); mkdir -p $(MODDIR); \
        for dir in $(DIRS); do cd $$dir && $(MAKE) $(MAKEFLAGS) && cd ..; done

clean:
		@for dir in $(DIRS); do cd $$dir && $(MAKE) clean && cd ..; done
