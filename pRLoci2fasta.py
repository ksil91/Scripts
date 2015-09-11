__author__ = 'ksil91'

#for older versions of pyRAD .loci

def Loci2fasta(infileName,outfileName):

    lociDict = {}
    ## read in your data file
    infile = open(infileName, "r")


    ## create an empty output file
    outfile = open(outfileName, "w")
    # outform = open("alph_form4.fasta", "w")
    # outest = open("alph_est4.fasta", "w")
    # outist = open("alph_ist4.fasta", "w")
    # outmal = open("alph_mal4.fasta", "w")
    # outcol = open("alph_col4.fasta", "w")

    ## parse the loci in your data file
    loci = infile.read().split("|\n")[:-1]

## write the first read from each locus to the output file in fasta format
locinum = 1en(loci)
for loc in loci:
    lociDict[locinum] = []
    reads = loc.split("\n")
    for r in range(0,len(reads)-1):
        name, seq = reads[r].split()
        sps = name.split("_")
        sp = sps[2]
        if sp not in lociDict[locinum]:
            newseq = seq.replace("-","")
            newname = name.replace(">","")
            print >>outfile, ">"+str(locinum)+"|"+newname+"|"+sp+"\n"+newseq
            lociDict[locinum].append(sp)
    locinum += 1


## close the output file
outfile.close()
infile.close()

locinum = 1
for loc in loci:
    lociDict[locinum] = []
    reads = loc.split("\n")
    for r in range(0,len(reads)-1):
        name, seq = reads[r].split()
        sps = name.split("_")
        sp = sps[2]
        if sp == "pan":
            newseq = seq.replace("-","")
            newname = name.replace(">","")
            print >>outfile, ">"+str(locinum)+"|"+newname+"|"+sp+"\n"+newseq
            lociDict[locinum].append(sp)
    locinum += 1