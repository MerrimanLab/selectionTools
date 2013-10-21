#!/bin/sh
# Installs the tabix program
#
bunzip2 -k src/tabix.tar.bz2
tar xf src/tabix.tar
rm -Rf src/tabix.tar
(cd tabix-0.2.6/ && make)
cp tabix-0.2.6/bgzip bin/
cp tabix-0.2.6/tabix bin/
rm -Rf tabix-0.2.6
