__author__ = 'ksil91'

#!/usr/bin/python

"""To convert pyRAD .loci output to individual nexus files for each locus.
"""

__author__ = 'Katherine Silliman'
__version__ = '0.0.3'

#imports
import sys
import os
from os import listdir
from os.path import isfile, join

#constants

#functions

def Loci2gPhocs(infileName,outfileSuf):
    numLoci = 1
    numTaxa = 0
    seqList = []
    outfile = open("temp.txt", "w")
    with open(infileName, "r") as INFILE:
        for line in INFILE:
            if "//" in line:
                outfile.write("\n\nLocus_"+str(numLoci)+" "+str(numTaxa) +" "+ str(numChar)+"\n")
                for seq in seqList:
                    outfile.write(seq)
                numTaxa = 0
                seqList = []
                numLoci +=1
            else:
                newline = line[1:]
                name = newline.split()
                nnn = newline.replace("-","N")
                seqList.append(nnn)
                numChar = len(name[1])
                numTaxa +=1
    INFILE.close()
    outfile.close()
    temp = open("temp.txt", "r")
    finalOut = open(outfileSuf+".txt","w")
    finalOut.write(str(numLoci))
    finalOut.write(temp.read())
    finalOut.close()
    temp.close()


def main(argv):

    #get arguments from command line

    infile_name = argv[1]
    outfile_Suf = argv[2]
    Loci2gPhocs(infile_name,outfile_Suf)


if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)