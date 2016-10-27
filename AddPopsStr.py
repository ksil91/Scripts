import sys

def addPops(str_infile, popfile, outfile, pop2int):
    IN = open(str_infile, "r")
    OUT = open(outfile, "w")
    pops = open(popfile, "r")
    pop2int = open(pop2int, "r")
    popdict = {}
    pop2intdict = {}
    for line in pop2int:
        popID = line.split()[0]
        intID = line.strip().split()[1]
        pop2intdict[popID] = intID
    for line in pops:
        sampleID = line.split()[0]
        popID = line.strip().split()[1]
        popdict[sampleID] = popID
    for line in IN:
        linelist = line.split()
        intID = pop2intdict[popdict[linelist[0]]]
        linelist.insert(1, intID)
        print >> OUT, "\t".join(str(e) for e in linelist)
    IN.close()
    OUT.close()
    pops.close()
    pop2int.close()

def main(argv):
    #get arguments from command line
    inf = argv[1]
    outf = argv[2]
    popf = argv[3]
    pop2int = argv[4]
    addPops(inf,popf, outf, pop2int)

if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)

