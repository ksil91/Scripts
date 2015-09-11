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

def Loci2IndNex(infileName,outfileSuf,GivennumTaxa):
    numLoci = 1
    taxaList = []
    outfile = open("temp.nex", "w")
    with open(infileName, "r") as INFILE:
        for line in INFILE:
            if "//" in line:
                outfile.close()
                finalOut = open("{0}{1}".format(outfileSuf+str(numLoci),".nex"), "w")
                finalOut.write("#NEXUS\nBEGIN DATA;\n\tDIMENSIONS NTAX="+str(GivennumTaxa)+ " NCHAR="+str(numChar)+";\n\tFORMAT " \
                "DATATYPE=DNA MISSING=N GAP=- INTERLEAVE=NO;\n\tMATRIX\n")
                temp = open("temp.nex","r")
                finalOut.write(temp.read())
                finalOut.close()
                temp.close()
                os.remove("temp.nex")
                outfile = open("temp.nex","w")
                numLoci+=1
            else:
                newline = line[1:]
                name = newline.split()
                if name[0] not in taxaList:
                    taxaList.append(name[0])
                outfile.write(newline)
                numChar = len(name[1])

    os.remove("temp.nex")
    INFILE.close()
    onlyfiles = [ f for f in listdir('.') if isfile(join('.',f)) ]
    for fil in onlyfiles:
        if "nex" in fil:
            editNames = []
            with open(fil, 'r+') as edit:
                for l in xrange(5):
                    next(edit)
                for lin in edit:
                    names = lin.split()
                    if names[0] not in editNames:
                        editNames.append(names[0])
                    ncharE = len(names[1])
                staxaList = set(taxaList)
                seditNames = set(editNames)
                missing = staxaList.difference(seditNames)
                if len(missing) != 0:
                    Nstring = "N" * ncharE
                    for n in missing:
                        edit.write(n+ '      ' +Nstring+'\n')
                edit.write(';\nEnd;')
            edit.close()


def main(argv):

    #get arguments from command line
    for dirpath, dirnames, files in os.walk('.'):
    	infile_name = argv[1]
    	outfile_Suf = argv[2]
    	num_Taxa = argv[3]
    	Loci2IndNex(infile_name,outfile_Suf,num_Taxa)
        break

if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)
