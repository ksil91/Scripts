__author__ = 'ksil91'

INFILE = open("Alpheus2c2.alleles","r")
new_loci = open("A2c2_sp3.gphocs", "w")
current_chunk = ""
num_species = 0
species = set()
tot_loci = 0
locus_tax = 0
for line in INFILE:
    if "//" in line:
        num_species = len(species)
        if num_species >= 3:
            new_loci.write("locus"+str(tot_loci)+" "+str(locus_tax)+" "+str(seq_len)+"\n")
            new_loci.write(current_chunk+"\n")
            tot_loci += 1
        current_chunk = ""
        locus_tax = 0
        species = set()
    else:
        locus_tax += 1
        newline = line[1:].replace("-","N")
        name = newline.split()
        seq_len = len(name[1])
        species.add(name[0].split("_")[2])
        current_chunk += newline

print str(tot_loci)
INFILE.close()
new_loci.close()
