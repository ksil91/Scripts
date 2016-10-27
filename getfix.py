__author__ = 'ksil91'


import sys
import numpy

def getLength(gphocs):
    gFile = open(gphocs,"r")
    locinum = 0
    specs = set()
    specD = {"colest": 0,"formpan": 0,"istmal": 0,"colmal":0,"colform":0,"colpan": 0,"colist":0,"estmal":0,"estform": 0,\
    "estpan": 0, "estist":0,"formist":0,"formmal": 0,"panist":0,"panmal":0} 
    gFile.next()
    sq = gFile.next()
    sqlen = int(sq.strip().split()[2])
    totseqlen = 0
    for line in gFile:
        if "locus" in line:
            if "col" in specs and "est" in specs:
                specD["colest"] += sqlen
            if "form" in specs and "pan" in specs:
                specD["formpan"] += sqlen
            if "ist" in specs and "mal" in specs:
                specD["istmal"] += sqlen
            if "col" in specs and "mal" in specs:
                specD["colmal"] += sqlen
            if "col" in specs and "form" in specs:
                specD["colform"] += sqlen
            if "col" in specs and "pan" in specs:
                specD["colpan"] += sqlen
            if "col" in specs and "ist" in specs:
                specD["colist"] += sqlen
            if "est" in specs and "mal" in specs:
                specD["estmal"] += sqlen
            if "est" in specs and "form" in specs:
                specD["estform"] += sqlen
            if "est" in specs and "pan" in specs:
                specD["estpan"] += sqlen
            if "est" in specs and "ist" in specs:
                specD["estist"] += sqlen
            if "form" in specs and "ist" in specs:
                specD["formist"] += sqlen
            if "form" in specs and "mal" in specs:
                specD["formmal"] += sqlen
            if "pan" in specs and "mal" in specs:
                specD["panmal"] += sqlen
            if "pan" in specs and "ist" in specs:
                specD["panist"] += sqlen                      
            totseqlen += sqlen            
            sqlen = int(line.strip().split()[2])
            specs = set()
        elif len(line) > 2:
            name = line.split()[0].split("_")[2]
            specs.add(name)
    gFile.close()    
    print specD
    print "Total nucs: "+str(totseqlen)



def getFix(snps,samples1,samples2):
    sFile = open(snps, "r")
    samps1 = open(samples1,"r")
    samps2 = open(samples2,"r")
    sampList1 = []
    sampList2 = []
    for line in samps1:
        sampList1.append(line.strip())
    for line in samps2:
        sampList2.append(line.strip())
    samps1.close()
    samps2.close()
    
    sDict1 = {}
    sDict2 = {}
    fixed = []
    fixcount = 0
    s1fixL = []
    s2fixL = []
    for line in sFile:
        if "##" not in line:
            s = line.split()
            snps = ''.join(s[1:])
            numsnps = len(snps)                                      
            if s[0] in sampList1:
                sDict1[s[0]] = snps
            elif s[0] in sampList2:
                sDict2[s[0]] = snps
    sFile.close()
    for i in range(0,numsnps):
        fix1 = False
        fix2 = False        
        s1 = []
        s2 = []
        for s in sDict1.keys():
            if sDict1[s][i] in list("ATGC"):
                s1.append(sDict1[s][i])
            elif "_" in sDict1[s][i]:
                s1 = []                
                break
            elif sDict1[s][i] in list("RKYSWM"):
                s1 = []
                break
        if len(s1) > 0:
            s1set = set(s1)
            if len(s1set) == 1:
                s1fix = s1set.pop()
                fix1 = True            
        if fix1:
            for s in sDict2.keys():
                if sDict2[s][i] in list("ATGC"):
                    s2.append(sDict2[s][i])
                elif sDict2[s][i] in list("RKYSWM"):
                    s2 = []
                    break
            if len(s2) > 0:            
                s2set = set(s2)
                if len(s2set) == 1:
                    s2fix = s2set.pop()
                    fix2 = True   
        if fix1 and fix2:
            if s2fix != s1fix:
                fixcount += 1
                s1fixL.append(s1fix)
                s2fixL.append(s2fix)
                fixed.append(i)
    print fixed[0:10]
    print "Total fixed: "+ str(fixcount)
    print "s1 fix: "+ str(s1fixL[0:10])
    print "s2 fix: "+ str(s2fixL[0:10])      

def main(argv):
    #get arguments from command line
    snpfile = argv[1]
    gfile = argv[2]
    samplelist1 = argv[3]
    samplelist2 = argv[4]
    getFix(snpfile,samplelist1,samplelist2)
    getLength(gfile)


if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)
