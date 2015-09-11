__author__ = 'ksil91'

import copy
from copy import deepcopy
from string import maketrans
import re

alleleDict = {}
indDict = {}
indADict = {}

intab = "RYSWKM"
outtab1 = "ACGAGA"
ot2 = "GTCTTC"
trantab1 = maketrans(intab,outtab1)
trantab2 = maketrans(intab,ot2)


## read in your data file
infile = open("/home/ksil91/Projects/Alpheus/alpheus_git/pyRAD_2_12/alpheus_combo/combo_6/outfiles/alpheus_c2.snps", "r")

## create an empty output file
outfile = open("test_combo4.str", "w")

h1 = r'[RKWYSM]+'
h = r'[RKWYSM]{2,}'
for line in infile:
    if "##" not in line:
        locis = line.split()
        numLoci = len(locis)
        indDict[locis[0]] = deepcopy(locis)


for x in range(1,numLoci):
    acount = 0
    for inds in indDict.keys():
        a = deepcopy(indDict[inds][x])
        if a == "_":
            break
        hetero1 = re.search(h1,a)
        heterom = re.search(h,a)
        if inds not in indADict:
            indADict[inds] = [" "," "," "]
        if "N" in a:
            indADict[inds][0].append(-9)
            indADict[inds][1].append(-9)
            indADict[inds][2].append(x)
        elif hetero1:
            if heterom:
                indADict[inds][0].append(-9)
                indADict[inds][1].append(-9)
                indADict[inds][2].append(x)
            else:
                a1 = a.translate(trantab1)
                a2 = a.translate(trantab2)
                if a1 not in alleleDict:
                    acount += 1
                    alleleDict[a1] = acount
                if a2 not in alleleDict:
                    acount += 1
                    alleleDict[a2] = acount
                indADict[inds][0].append(deepcopy(alleleDict[a1]))
                indADict[inds][1].append(deepcopy(alleleDict[a2]))
                indADict[inds][2].append(x)
        else:
            if a not in alleleDict:
                acount += 1
                alleleDict[a] = acount
            indADict[inds][0].append(deepcopy(alleleDict[a]))
            indADict[inds][1].append(deepcopy(alleleDict[a]))
            indADict[inds][2].append(x)
    alleleDict = {}

infile.close()

SF = list(indADict.keys())
SF.sort()
print >> outfile, " ".join(str(e1) for e1 in indADict[SF[0]][2])
for i in SF:
    print >> outfile, i + " ".join(str(e) for e in indADict[i][0])
    print >> outfile, i +" ".join(str(ef) for ef in indADict[i][1])
outfile.close()



