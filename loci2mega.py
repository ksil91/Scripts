#!/usr/bin/python

"""To convert pyRAD .loci output to MEGA input
"""

__author__ = 'Katherine Silliman'
__version__ = '0.0.1'

#imports
import sys
import os
from os import listdir
from os.path import isfile, join

#constants

#functions

def Loci2Mega(infileName,outfileName,genes,genesFil,taxaFil):
    numLoci = 1
    thisName = []
    pidsList = []
    taxaList = open(taxaFil, "r")
    tList = []
    finalOut = open("{0}{1}".format(outfileName,".meg"), "w")
    locitemp = open("locitemp.meg", "w")
    if genes == "both":
        prot_ids = open(genesFil, "r")
        for p in prot_ids:
            pidsList.append(p.strip())
        prot_ids.close()

    for t in taxaList:
        tList.append(t.strip())
    taxaList.close()
    tset = set(tList)

    with open(infileName, "r") as f:
        text = f.read().strip()
        totLoci = text.count("|")
    f.close()

    finalOut.write("#mega\n!Title "+outfileName+";\n!Format DataType=DNA NTaxa="+str(len(tset))+" NSeqs="+str(totLoci)\
    +" Indel=- Missing=?;\n")
    with open(infileName, "r") as INFILE:
        for line in INFILE:
            if "//" in line:
                locitemp.close()
                locitemp = open("locitemp.meg", "r")
                locid = line.strip().split('|')[1]
                if genes == "both":
                    if locid in pidsList:
                        prop = "Coding"
                    else:
                        prop = "Noncoding"
                elif genes == "coding":
                    prop = "Coding"
                else:
                    prop = "Noncoding"
                finalOut.write("\n!Gene=loci"+str(numLoci)+" Property="+prop+";\n")
                finalOut.write(locitemp.read())
                locitemp.close()
                os.remove("locitemp.meg")
                nameset = set(thisName)
                missing = tset.difference(nameset)
                if len(missing) != 0:
                    Nstring = "?"*numChar
                    for n in missing:
                        finalOut.write("#"+n+"_{species."+n.split("_")[2]+"}    "+Nstring+"\n")
                locitemp = open("locitemp.meg","w")
                numLoci+=1
                thisName = []
            else:
                newline = line[1:]
                name = newline.split()
                seqed = name[1].replace("N","?")
                thisName.append(name[0])
                locitemp.write("#"+name[0]+"_{species."+name[0].split("_")[2]+"}    "+seqed+"\n")
                numChar = len(seqed)

    os.remove("locitemp.meg")
    INFILE.close()


def main(argv):
    #get arguments from command line
    infile_name = argv[1]
    outfile_name = argv[2]
    genes = argv[3]
    genes_file=argv[4]
    taxa_file = argv[5]
    Loci2Mega(infile_name,outfile_name,genes,genes_file,taxa_file)

if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)
