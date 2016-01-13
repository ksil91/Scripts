__author__ = 'ksil91'

blast_out = open("23_6_blast.out","r")
rp_prot = open("rna_prot.fasta","r")
gene_loc = open("23_6_genes.out","w")

ref_dict = {}

for line in blast_out:
    if "blastn" not in line:
        names = line.split(",")
        ref = names[1].split("_")
        ref_id = ref[3]
        query = names[0].split("_")
        q_id = query[3]
        ref_dict[ref_id] = q_id


for line in rp_prot:
    if ">" in line:
        names = line.split("|")
        full_name = names[0].split("_")
        rp_id = full_name[3]
        rp_id_23_6 = ref_dict.get(rp_id, "poop")
        if rp_id_23_6 != "poop":
            print >> gene_loc, str(rp_id_23_6)


gene_loc.close()
rp_prot.close()
blast_out.close()

f = open("raxml_23_6.loci","r")
f_stuff = f.read()
f.close()

import re

new_loci = re.split(r'\|\', f_stuff)



