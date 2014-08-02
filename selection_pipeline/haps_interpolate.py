"""
    Annotates a haps file by replacing all the SNPS with the genetic map positions
    Assumes the genetic map is in the shapeit format. 
"""

import argparse
from collections import OrderedDict
import tempfile
import decimal
from decimal import *

def get_genetic_map_format(genetic_map):
    """ Method to get the format of the genetic map files 
        used to work out whether we are dealing with shapeit or
        plink geneticmap
    """
    with open(genetic_map) as gmap:
        gmap_line = gmap.readline()
        # If it's 1 it must be a shapeit genetic map I think
        # Because a plink genetic map is tab seperated
        if(len(gmap_line.split('\t'))==1):
            return "shapeit"
        else:
            return "plink"

def plink_to_shapeit_gmap(genetic_map_output,new_genetic_map):
    """ Convert a plink genetic map to shapeit, as this normalises our
        genetic map for use further in the pipeline,

        e.g Enables shapeit to use a plink geneticmap
        User might prefer beagle and not use shapeit but
        we use a plink genetic map for rehh
    """   
    with open(genetic_map_output) as gmap:
        for i, line in enumerate(gmap):
            gmap_line = line.split('\t')
            if ( i == 0 ):
                continue
            else:
                position=int(gmap_line[1])
                recomb_rate=float(gmap_line[2])
                centi_morgans=float(gmap_line[3])
                new_genetic_map.write(str(position) + ' ' + str(recomb_rate) + ' ' + str(centi_morgans) + '\n')
    return new_genetic_map

def get_shapeit_genetic_map(genetic_map,temp_genetic_map):
    """ Returns either the original file name
        if the file is already in shapeit format 
    """
    file_format=get_genetic_map_format(genetic_map)
    if(file_format =='shapeit'):
        return(genetic_map)
    else:
        temp_genetic_map = open(temp_genetic_map,'w')
        plink_to_shapeit_gmap(genetic_map,temp_genetic_map)
        return(temp_genetic_map)

def load_genetic_map(genetic_map):
    gmap_pos = OrderedDict()
    genetic_map.seek(0)
    for i, line in enumerate(genetic_map):
        #shape it format is line seperated
        shapeit_line = line.split()
        gmap_pos[float(shapeit_line[0])]=Decimal(shapeit_line[2])
    return gmap_pos

def interpolate(start_position,end_position,x):
    start0 = Decimal(str(start_position[0]))
    end0 = Decimal(str(end_position[0]))
    slope = (end_position[1] - start_position[1])/(end0 - start0) 
    intercept=start_position[1]
    interp = intercept + ((x-start0) * slope)
    return interp

def replace_positions(haps,output,gmap_dict,physical_out):
    interpolate_list = []
    out = open(output,'w')
    phys_out = None
    if (physical_out is not None):
        phys_out = open(physical_out,'w')
    with open(haps) as f:
        gmap_dict = gmap_dict.items()
        haps_line = f.readline()
        dictionary_index = 1
        # Requires the gmap_dictionary is atleast  
        start_position=[0,Decimal("0.0")]
        end_position=gmap_dict[0]
        while(haps_line and dictionary_index < len(gmap_dict)):
            temp_line = haps_line.split()
            temp_pos = int(temp_line[2])
            if(temp_pos >= start_position[0] and temp_pos <= end_position[0]):
                t_inter = interpolate(start_position,end_position,temp_pos)
                temp_line[2] = str(t_inter)
                temp_inter = ' '.join(temp_line)
                out.write(temp_inter + '\n')
                if (phys_out != None):
                    phys_out.write(str(temp_pos)+ '\n')
                haps_line = f.readline()
            else:
                start_position = end_position
                end_position = gmap_dict[dictionary_index]
                dictionary_index += 1
    out.close()
    return interpolate_list

def main():
    parser = argparse.ArgumentParser(description="Annotate a haps file with Genetic map positions using linear interpolation")
    parser.add_argument('--haps',dest='haps')
    parser.add_argument('--output',dest='output')
    parser.add_argument('--genetic-map',dest="gmap")
    parser.add_argument('--physical-position-output',dest='physical_positions')
    args = parser.parse_args()
    assert args.haps is not None, \
        "Haps file needs to be specified"
    assert args.output is not None, \
        "Output file needs to be specified"
    assert args.gmap is not None, \
        "Genetic map needs to be specified"
    haps = args.haps
    output = args.output
    genetic_map = args.gmap
    physical_out = args.physical_positions
    file_format = get_genetic_map_format(genetic_map)
    if(file_format == 'plink'):
        temp_file = tempfile.TemporaryFile()
        genetic_map=plink_to_shapeit_gmap(genetic_map,temp_file)    
    else:
        genetic_map = open(genetic_map)
    gmap_dict = load_genetic_map(genetic_map)
    replace_positions(haps, output, gmap_dict, physical_out=physical_out)

if __name__ =="__main__":
    main()
