#!/usr/bin/env python
#
# Script for concatenating haps files. 
#
# @author James Boocock
# @date 09/08/2013
#

from optparse import OptionParser
import logging

logging.basicConfig(format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)
log.setLevel(level=logging.INFO)


def main():
    parser = OptionParser()
    parser.add_option("-i",'--file',dest="haps_files",help="Base names of haps / sample pairs")
    parser.add_option('-v','--verbose',action="store_true",help="Verbosity of the logging outputs")
    (options,args) = parser.parse_args()
    if(options.verbose == None):
        log.setLevel(logging.ERROR)
    else:
        log.setLevel(logging.DEBUG)
    haps_files = options.haps_files
    print(haps_files)
     


if __name__=="__main__":main()

