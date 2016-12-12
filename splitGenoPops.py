##Creates a new .geno and .ind file for each pair combination of populations. This is to be used for Fst and Pi scripts.

import sys
from itertools import combinations


def getPopList(popfile):
    pop = open(popfile, "r")
    poplist = []
    samplelist = []
    for line in pop:
        sampleID, U, popID = line.strip().split()
        poplist.append(popID)
        samplelist.append(sampleID)
    pop.close()
    return(poplist,samplelist)

def makeNewGenoPair(suf,pop12,poplist,samplelist):
    IN = open(suf+".geno","r")
    OUT = open(suf+"_"+pop12[0]+"_"+pop12[1]+".geno","w")
    indfile = open(suf+"_"+pop12[0]+"_"+pop12[1]+".ind","w")
    for line in IN:
        genos = list(line.strip())
        for i in range(0,len(genos)):
            if poplist[i] in pop12:
                OUT.write(genos[i])
        OUT.write("\n")
    for i in range(0,len(samplelist)):
        if poplist[i] == pop12[0]:
            indfile.write(samplelist[i]+"    0\n")
        elif poplist[i] == pop12[1]:
            indfile.write(samplelist[i]+"    1\n")
    IN.close()
    OUT.close()
    indfile.close()

def makeNewGenoSingle(suf,pop,poplist,samplelist):
    IN = open(suf+".geno","r")
    OUT = open(suf+"_"+pop+".geno","w")
    indfile = open(suf+"_"+pop+".ind","w")
    for line in IN:
        genos = list(line.strip())
        for i in range(0,len(genos)):
            if poplist[i] == pop:
                OUT.write(genos[i])
        OUT.write("\n")
    for i in range(0,len(samplelist)):
        if poplist[i] == pop:
            indfile.write(samplelist[i]+"    0\n")
    IN.close()
    OUT.close()
    indfile.close()


def main(argv):
    #get arguments from command line
    suf = argv[1]
    popfile = argv[2]
    poplist,samplelist = getPopList(popfile)
    popset = set(poplist)
    #for c in combinations(popset,2):
        #makeNewGeno(suf,c,poplist,samplelist)
    for s in popset:
        makeNewGenoSingle(suf,s,poplist,samplelist)
    

if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)

