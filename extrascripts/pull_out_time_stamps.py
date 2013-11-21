#!/usr/bin/python

import os
from argparse import ArgumentParser


def main():
    parser=ArgumentParser(description='Generate summary information from multi_populatand selection_pipeline config_files.')
    parser.add_argument('input_configs',metavar='N',type=str,nargs='+', help='config_files_to_summarise')
    args = parser.parser_args()

if __name__=="__main__":main()
