"""
    Annotates a haps file by replacing all the SNPS with the genetic map positions
    Assumes the genetic map is in the shapeit format. 
"""

import argparse
import OrderedDict

def load_genetic_map(genetic_map,file_format):
    gmap_pos = OrderedDict()
    for i, line in enumerate(f):
        #shape it format is line seperated
        shapeit_line = line.split()
        gmap_pos[int(shapeit_line[0])]=float(shapeit_line[2])

    return gmap_list

def replace_positions(haps,output,genetic_map):
    with open(haps) as f:
        for line in f:
        

def main():
    parser = argparse.ArgumentParser(description="Annotate a haps file with Genetic map positions")
    parser.add_argument('--haps',dest='haps')
    parser.add_argument('--output',dest='output')
    parser.add_argument('--genetic-map',dest="gmap")
    args = parser.parse_args()
    assert args.haps is not None, \
        "Haps file needs to be specified"
    assert args.output is not None, \
        "Output file needs to be specified"
    assert args.gmap is not None, \
        "Genetic map needs to be specified"
    
    haps = arg.haps
    output = args.output
    genetic_map = args.gmap
    file_format = detect_file_format(genetic_map)
    gmap_dict = load_genetic_map(genetic_map,file_format)
    replace_positions(haps, output, genetic_map)
