#!/usr/bin/python

__author__ = 'Katherine Silliman'
__version__ = '0.0.1'

#imports
import re
import copy
from copy import deepcopy

#constants
#def findLoci(coverage):
f = open("Projects/alpheus-gbs/pyRAD_2_12/alpheus_combo/combo_6/outfiles/alpheus_c2.loci","r")
file_loci = f.read()
f.close()
    
istList = []
covDict = {} 
#all 3
#reg1 = r'98_52_ist'
#reg2 = r'98_53_ist'
#reg3 = r'98_54_ist'
#mincov = 8
#istList =[reg1,reg2,reg3]
#2
#reg1 = r'98_52_ist'
#reg2 = r'98_53_ist'
#istList = [reg1,reg2]
#mincov = 7
#1
reg1 = r'98_52_ist'
istList.append(reg1)
mincov = 6 
    
all_loci = re.split('//', file_loci)
nameReg = r'>\w*'
nameList = []
mistList = []
for loci in all_loci:
    cov = loci.count('>')
    if cov >= mincov:
        for reg in istList:
            m = re.search(reg, loci)
            mistList.append(m)
        if None not in mistList:
            nameList = re.findall(nameReg,loci)
            names = tuple(deepcopy(nameList))
            if names not in covDict:
                covDict[names] = 1
            else:
                covDict[names] +=1
        mistList = []
        nameList = []
 
#print covDict
dkeys = covDict.keys()
print str(len(dkeys))
maxCov = 0
maxCombo = []
for combo in covDict.keys():
    if covDict[combo] > maxCov:
        maxCov = covDict[combo]
        maxCombo = deepcopy(combo)
print "maxCov: " + str(maxCov)
print "maxCombo: " + str(maxCombo)
                
                
   
#        least loci:
#            06_271_ist (4 total)
#            98_163_form (5)
#            D_235_pan (5)
#            C330_F177_col (4)
#            98_308_est (4)
#            98_216_mal (499) (2)
#    