#
# Murray Cadzow
# December 2015
# University of Otago
#
#
# A script to convert the output of haps_to_selscan
# into files with matched markers for running selscan --xpehh

from optparse import OptionParser

def selscan_xpehh_convert(options):
    pop1_map = open(options.pref1 + ".selscanmap", 'r')
    pop2_map = open(options.pref2 + ".selscanmap",'r')

    # read in positons from map file
    pos1 = {}
    i = 0
    for line in pop1_map:
        pos1[int(line.split()[3])] = i #store position with line number (0 indexed)
        i += 1

    pos2 = {}
    i = 0
    for line in pop2_map:
        pos2[int(line.split()[3])] = i  #store position with line number (0 indexed)
        i += 1

    # create intersection of map file positions
    dict1 ={}
    dict2 = {}
    for ele in pos1.keys():
        dict1[ele] = 1

    for ele in pos2.keys():
        if ele in dict1:
            dict2[ele] = 1

    intersect = dict2.keys()
    # create two lists of indexes which will correspond to line numbers
    # of what lines to keep
    extract1 = []
    for ele in intersect:
        extract1.append(pos1[ele])
    extract1=sorted(extract1)

    extract2 = []
    for ele in intersect:
        extract2.append(pos2[ele])
    extract2 = sorted(extract2)

    # write out common map file
    # positions should be the same if using pop1_map or pop2_map but marker names may be different
    out1_map = open(options.out +options.pop1 +'_' + options.pop2 +'_'+ options.chr+ '.xpehh_selscanmap','w')
    pop1_map.seek(0)
    i = 0
    for ind in extract1:
        while i < ind:
            pop1_map.readline()
            i += 1
        if i == ind:
            i+=1
            out1_map.write(pop1_map.readline())
    out1_map.close()
    pop1_map.close()

    # sanity check - write second map file from pop2_map, positions should match exactly between the 2 files
    #out1_map = open(options.out +options.pop1 +'_' + options.pop2 +'_'+ options.chr+ '.xpehh_selscanmap','w')
    #pop2_map.seek(0)
    #i = 0
    #for ind in extract2:
    #    while i < ind:
    #        pop2_map.readline()
    #        i += 1
    #    if i == ind:
    #        i+=1
    #        out2_map.write(pop2_map.readline())
    #out2_map.close()
    #pop2_map.close()

    # write out selected hap1 columns
    out1_hap = open(options.out + options.pop1 + "_" + options.chr + '.matches_' + options.pop2+'.xpehh_selscanhaps','w')
    pop1_hap = open(options.pref1 + ".selscanhaps", 'r')
    for line in pop1_hap:
        hap = line.split()
        new_hap = [hap[i] for i in extract1]
        out1_hap.write(' '.join(new_hap) + '\n')
    out1_hap.close()
    pop1_hap.close()

    #write out selected hap2 columns
    out2_hap = open(options.out + options.pop2 + "_" + options.chr + '.matches_' + options.pop1 +'.xp_ehh_selscanhaps','w')
    pop2_hap = open(options.pref2 + ".selscanhaps",'r')
    for line in pop2_hap:
        hap = line.split()
        new_hap = [hap[i] for i in extract2]
        out2_hap.write(' '.join(new_hap)+ '\n')
    out2_hap.close()
    pop2_hap.close()




def main():
    parser = OptionParser()
    parser.add_option('--pop1-prefix', dest='pref1', help = 'prefix for .selscanhaps and selscanmap files for population 1')
    parser.add_option('--pop2-prefix', dest='pref2', help = 'prefix for .selscanhaps and selscanmap files for population 2 ')
    parser.add_option('--out', dest = 'outpref', help = 'prefix for output files')
    parser.add_option('-c', dest = 'chr', help = 'Chromosome')
    parser.add_option('--pop1-name', dest = 'pop1', help = 'name of population 1')
    parser.add_option('--pop2-name', dest = 'pop2', help = 'name of population 2')
    (options, args) = parser.parse_args()

    assert options.pref1 is not None and options.pref2 is not None,\
        "prefixes for both population files is required"
    assert options.pop1 is not None and options.pop2 is not None,\
        "Population names are required"
    assert options.chr is not None,\
        "chromosome is required"
    assert options.outpref is not None,\
        "output prefix is required"

    selscan_xpehh_convert(options)

if __name__ == "__main__":
    main()
