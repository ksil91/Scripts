loci = open("raxml_23_6.loci","r")
pro_ids = open("23_6_genes.out","r")
pro_loci = open("pro_23_6.loci", "w")
nonpro_loci = open("non_23_6.loci", "w")
current_chunk = ""
id_set = set()

for l in pro_ids:
    id_set.add(int(l))
pro_ids.close()

for line in loci:
    if line[:2] == "//":
        loc_id = int(line.split("|")[1])
        if loc_id in id_set:
            current_chunk+= line                        
            pro_loci.write(current_chunk)
        else:
            current_chunk += line            
            nonpro_loci.write(current_chunk)
        current_chunk = ""
    else:
        current_chunk += line

loci.close()
pro_loci.close()
nonpro_loci.close()
            
		


