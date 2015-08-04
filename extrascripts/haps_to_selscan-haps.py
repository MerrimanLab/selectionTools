import argparse

def readFiles(haps, pos, output, chr):
    hapsFile = open(haps, 'r')
    posFile = open(pos, 'r')
    output_haps = open(output + ".selscanhaps", 'w')
    output_map = open(output + ".selscanmap", 'w')
    
    hapsArray = []
    #posArray = []
    positions = []
    for hapsLine in hapsFile:
        #output_haps.write( changeLine(hapsLine.strip()) )
        posLine=posFile.readline().strip().split(' ')
        hapsLine = hapsLine.strip().split(' ')
        if(len(posLine) == 1):
            i=0
        else:
            i=1
        if(posLine[i] not in positions):
            positions.append(posLine[i])
            hapsArray.append(hapsLine[5:])
            output_map.write(str(chr) +' '+ hapsLine[0] +' '+ hapsLine[2] +' '+ posLine[i]+ "\n")
            #posArray.append(posLine)
    #transpose haps array and write out
    #each row is now a sample (haplotype), each coloum is now a loci
    for i in xrange(0,len(hapsArray[0])): # number of loci
        line = []
        for j in xrange(0,len(hapsArray)): # number of samples
            line.append(hapsArray[j][i])    
        output_haps.write(" ".join(line) + "\n")     
    #close files    
    output_map.close()
    output_haps.close()
    hapsFile.close()
    posFile.close()




def changeLine(line):
    l = line.split(' ')
    hap1 = []
    hap2 = []
    #haps file has marker marker pos a1 a2 hap0 hap0 hap1 hap1...
    for i in xrange(5,len(l)):
        if(i % 2 == 0):
            hap1.append(l[i])
        elif (i % 2 == 1):
            hap2.append(l[i])
    newHaps = " ".join(hap1) + "\n" + " ".join(hap2) + "\n"
    return newHaps

def main():
    parser = argparse.ArgumentParser(description="Change haps file to selscan haps file\nInput files are ouput of haps_interpolate")
    parser.add_argument('--haps',dest='haps',help="haps file from haps_interpolate")
    parser.add_argument('--pos',dest='pos', help="pos file from haps_interpolate")
    parser.add_argument('--output',dest='output',help='prefix for output files')
    parser.add_argument('--chr', dest='chr', help='chomosome')
    args = parser.parse_args()
    assert args.haps is not None, \
            "Haps file needs to be specified"
    assert args.pos is not None, \
            "Pos file needs to be specified"
    assert args.output is not None, \
            "Output file needs to be specified"
    assert args.chr is not None, \
            "chromosome needs to be specified"
    readFiles(args.haps, args.pos, args.output, args.chr)
if __name__ =="__main__":
    main()
