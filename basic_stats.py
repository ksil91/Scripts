#!/usr/bin/python

"""To get basic popgen stats from MEGA output as produced by loci2mega.py.
"""

__author__ = 'Katherine Silliman'
__version__ = '0.0.1'

#imports
import re
import numpy
import sys
import getopt

#constants

#functions

def buildNameDict(infile):
    namedict = {}
    with open(infile, "r") as infil:
        loc_r = re.compile(r'(?:Gene=)(\w*)')
        for l in xrange(4):
            next(infil)
        for l in infil:
            if "!" in l:
                locid = re.findall(loc_r,l)
            elif "#" in l:
                (name,seq) = l[1:].split()
                if name not in namedict.keys():
                    namedict[name] = {}
                if "NNNNN" not in seq:
                    namedict[name][locid[0]] = seq
                else:
                    namedict[name][locid[0]] = None
    infil.close()
    return namedict

def nucratio(s):
    agct = [s.count("A"), s.count("G"), s.count("C"), s.count("T")]
    return agct

def findOcc(s,ch):
    return [i for i, letter in enumerate(s) if letter == ch]

def hamdist(str1, str2):
    diffs = 0
    naset = set().union(set(findOcc(str1,"N")), set(findOcc(str1,"-")),set(findOcc(str2,"N")),set(findOcc(str2,"-")))
    news1 = "".join([char for idx, char in enumerate(str1) if idx not in naset ])
    news2 = "".join([char for idx, char in enumerate(str2) if idx not in naset ])

    for ch1, ch2 in zip(news1, news2):
        if ch1 != ch2:
            diffs += 1
    return [diffs,len(news1)]


def getGroups(gfile):
    gf = open(gfile,"r")
    gdict = {}
    for l in gf:
        (name,pop) = l.split()
        if pop not in gdict:
            gdict[pop] = [name]
        else:
            gdict[pop].append(name)
    return gdict


def withindist(ndict,gdict):
    widict = {}
    for pop in gdict.keys():
        compset = set()
        plist = []
        nuclen = []
        indlist = gdict[pop]
        for ind in indlist:
            totseqlen = 0
            totdist = 0
            locids = ndict[ind].keys()
            for comp in indlist:
                if (ind != comp) and (comp not in compset):
                    for lid in locids:
                        indseq = ndict[ind][lid]
                        compseq = ndict[comp][lid]
                        if (indseq != None) and (compseq != None):
                            dlist = hamdist(indseq,compseq)
                            seqlen = dlist[1]
                            totdist += dlist[0]
                            totseqlen += seqlen
                    plist.append(float(totdist)/float(totseqlen))
                    nuclen.append(totseqlen)
            compset.add(ind)
        widict[pop] = [numpy.mean(plist,dtype=numpy.float64),numpy.mean(nuclen),numpy.std(nuclen)]
    return widict


def main(argv):
    inputfile = ''
    outputfile = ''
    popfile = ''
    method = ''

    try:
        opts, args = getopt.getopt(argv,"hi:o:p:m:",["ifile=","ofile=","pfile=","m="])
    except getopt.GetoptError:
        print 'basic_stats.py -i <inputfile> -o <outputfile> -p <popfile> -m <within, nucratio>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'basic_stats.py -i <inputfile> -o <outputfile> -p <popfile> -m <within, nucratio>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt == "-p":
            popfile = arg
        elif opt == "-m":
            method = arg
    if method == "within":
        namedict = buildNameDict(inputfile)
        groupdict = getGroups(popfile)
        widict = withindist(namedict,groupdict)
        for p in widict.keys():
            print p+'\nAvg proportion of differences: '+ str(widict[p][0])+'\nMean # of nucleotides per comparision: '\
            + str(widict[p][1]) + '\nStd in # of nucleotides: '+ str(widict[p][2])
    if method 

if __name__ == "__main__":
   main(sys.argv[1:])