__author__ = 'ksil91'


import sys
import numpy

def unstruct(amb):
    amb = amb.upper()
    " returns bases from ambiguity code"
    D = {"R":["G","A"],
         "K":["G","T"],
         "S":["G","C"],
         "Y":["T","C"],
         "W":["T","A"],
         "M":["C","A"],
         "A":["A","A"],
         "T":["T","T"],
         "G":["G","G"],
         "C":["C","C"],
         "N":["N","N"],
         "-":["-","-"]}
    return D.get(amb)


def getLength(gphocs):
    gFile = open(gphocs,"r")
    lenList = []
    for line in gFile:
        if "locus" in line:
            sqlen = line.strip().split()[2]
            lenList.append(int(sqlen))
    gFile.close()
    return lenList


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
    return [float(diffs),len(news1)]


def getSNPs(snps,samples):
    sFile = open(snps, "r")
    samps = open(samples,"r")
    sampList = []
    for line in samps:
        sampList.append(line.strip())
    samps.close()
    sDictA = {}
    sDictB = {}
    sDlist = []
    for line in sFile:
        if "##" not in line:
            s = line.split()
            if s[0] in sampList:
                sDictA[s[0]] = []
                sDictB[s[0]] = []
                for snps in s[1:]:
                    if any((c in list("ATGCRKYSWM")) for c in snps):
                        snpA = ""
                        snpB = ""
                        for snp in snps:
                            if snp in list("RKYSWM"):
                                het = unstruct(snp)
                                snpA += het[0]
                                snpB += het[1]
                            else:
                                snpA += snp
                                snpB += snp
                        sDictA[s[0]].append(snpA)
                        sDictB[s[0]].append(snpB)
                    elif "_" not in snps:
                        sDictA[s[0]].append("none")
                        sDictB[s[0]].append("none")
                    else:
                        sDictA[s[0]].append("_")
                        sDictB[s[0]].append("_")

    sFile.close()
    sDlist.append(sDictA)
    sDlist.append(sDictB)
    return sDlist

def getPi(lD,sDlist,outfile):
    sDA = sDlist[0]
    sDB = sDlist[1]
    out = open(outfile, "w")
    names = sDA.keys()
    pi_alllist = []
    pi_nozeros = []
    for snp in range(0,len(sDA[names[0]])):
        if sDA[names[0]][snp] != "_":
            snplen = len(sDA[names[0]][snp])
            snpcomp = 0.0
            missingsamp = 0
            for x in range(0,len(names)):
                if sDA[names[x]][snp] != "none":
                    snpdist = hamdist(sDA[names[x]][snp],sDB[names[x]][snp])
                    missing = snplen - snpdist[1]
                    new_snpcomp = snpdist[0]/(lD[snp]-missing)
                    snpcomp += new_snpcomp
                    for y in range(x+1,len(names)):
                        if sDA[names[y]][snp] != "none":
                            snpdist = hamdist(sDA[names[x]][snp],sDA[names[y]][snp])
                            missing = snplen - snpdist[1]
                            new_snpcomp = snpdist[0]/(lD[snp]-missing)
                            snpcomp += new_snpcomp
                            snpdist = hamdist(sDA[names[x]][snp],sDB[names[y]][snp])
                            missing = snplen - snpdist[1]
                            new_snpcomp = snpdist[0]/(lD[snp]-missing)
                            snpcomp += new_snpcomp
                            snpdist = hamdist(sDB[names[x]][snp],sDA[names[y]][snp])
                            missing = snplen - snpdist[1]
                            new_snpcomp = snpdist[0]/(lD[snp]-missing)
                            snpcomp += new_snpcomp
                            snpdist = hamdist(sDB[names[x]][snp],sDB[names[y]][snp])
                            missing = snplen - snpdist[1]
                            new_snpcomp = snpdist[0]/(lD[snp]-missing)
                            snpcomp += new_snpcomp
                else:
                    missingsamp += 1
            if missingsamp == len(names):
                out.write("NaN\n")
            else:
                numseq = float(2*(len(names)-missingsamp))
                #print(numseq)
                #print(snpcomp)
                pi = (1.0/(numseq**2))*snpcomp
                pi_alllist.append(float(pi))
                out.write(str(pi)+"\n")
                if pi != 0.0:
                    pi_nozeros.append(float(pi))
        else:
            out.write("_\n")
            pi_alllist.append(0.0)

    out.close()
    pi_all_array = numpy.array(pi_alllist)
    pi_nz_array = numpy.array(pi_nozeros)
    pi_all_mean = numpy.mean(pi_all_array)
    pi_nz_mean = numpy.mean(pi_nz_array)
    print("Pi all: "+str(len(pi_alllist))+", "+str(pi_all_mean))
    print("Pi no zeros: "+str(len(pi_nozeros))+", "+str(pi_nz_mean))
    return pi_all_array


def main(argv):
    #get arguments from command line
    gphocs = argv[1]
    snpfile = argv[2]
    output = argv[3]
    samplelist = argv[4]
    lD = getLength(gphocs)
    sDL = getSNPs(snpfile,samplelist)
    pis = getPi(lD,sDL,output)

if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)
