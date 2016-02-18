__author__ = 'ksil91'

INFILE = open("Alpheus2c2.loci","r")
new_loci = open("A2c2_3spec.loci", "w")
current_chunk = ""
num_species = 0
species = set()
tot_loci = 0
for line in INFILE:
    if "//" in line:
        num_species = len(species)
        if num_species >= 3:
            current_chunk+= line
            new_loci.write(current_chunk)
            tot_loci += 1
        current_chunk = ""
        species = set()
    else:
        name = line.split()
        species.add(name[0].split("_")[2])
        current_chunk += line

print str(tot_loci)
INFILE.close()
new_loci.close()
