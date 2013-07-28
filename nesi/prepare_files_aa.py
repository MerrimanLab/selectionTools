

import sys
import os

#$1 haps file to split
#$2 Size of split files <window size>
#$3 overlap = 100000
#$4 population
# Remove final file

def main():
    window = int(sys.argv[2])
    overlap = int(sys.argv[3])
    population=(sys.argv[4])
    i = 1	
    overlapping_file = open(population+str(i+1)+'.phaps','w')
    first_file = open(population+str(i) + '.phaps','w')
    temp=overlap
    #print(window * i)
    with open(sys.argv[1],'r') as f:
        for line in f:
            #if(i == 1):
            #    temp=overlap
            #    overlap = 0
            #else:
            #    overlap=temp
    # Ignore offset because it doesnt matter
            if(int(line.split()[2]) < ((window -overlap) * i + overlap) and int(line.split()[2]) >= ((window - overlap ) * (i - 1))):
                first_file.write(line)
            elif(int(line.split()[2]) >= ((window -overlap)* (i) + overlap)):
                i = i + 1
                first_file.close()
                first_file = overlapping_file 
                overlapping_file = open(population+str(i+1)+'.phaps','w')
            if(int(line.split()[2]) >= ((window - overlap)*i)):
                overlapping_file.write(line)
    overlapping_file.close()
    first_file.close()
    print("Done Success")


if __name__=="__main__":main()

