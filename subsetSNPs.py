__author__ = 'Katherine Silliman'

## Code to subset one SNP per GBS locus from a VCF file. Chooses the SNP
## with the highest sample coverage. If there is a tie, chooses the 1st SNP in the loci. (may change to random)
## May be specific to VCF format output from ipyrad.

import sys
import linecache
def subsetSNPs(inputfile,outputfile):
    locidict = {}
    lineNum = []
    IN = open(inputfile, "r")
    OUT = open(outputfile, "w")

    n = 1
    for line in IN:
        if "#" not in line:
            linelist = line.split()
            if "loc" in linelist[0]:
                loci = linelist[0]
            else:
                loci = int(linelist[0])
            #Column 8 is INFO column of VCF file
            NS = float(linelist[7].split(";")[0].split("=")[1])
            if loci not in locidict.keys():
                locidict[loci] = [NS,n]
            else:
                if locidict[loci][0] < NS:
                    locidict[loci] = [NS,n]
        else:
            OUT.write(line)
        n += 1
    IN.close()
    print("Total SNPS: "+str(n)+"\nUnlinked SNPs: "+str(len(locidict.keys())))

    for locus in sorted(locidict.keys()):
        line = linecache.getline(inputfile, locidict[locus][1])
        OUT.write(line)
    OUT.close()

def main(argv):
    #get arguments from command line, first name of the infile .vcf then name of the outfile .vcf
    invcf = argv[1]
    outfile = argv[2]
    subsetSNPs(invcf, outfile)
    

if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)

