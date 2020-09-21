import bx.align.maf
import bx.align.tools
from Bio import AlignIO

########################################################
#RESTRICT THE ALIGNMENT TO REGIONS WITH ALL FIVE SPECIES
input_handle = open("chrY.ancestor.pruned.maf", "r")
output_handle = open("chrY.ancestor.pruned.5species.maf", "w")

alignments = AlignIO.parse(input_handle, "maf")

species_list=["Anc0", "panPan_Y", "ponAbe_Y", "panTro_Y", "gorGor_Y", "hg_Y"]
species_list.sort()
print(species_list)

for record in alignments :
    species_present=[]
    for row in record :
        #print(row.seq)
        #print(row.id)
        species=(row.id).split(".",1)[0]
        #if species not in species_present:
        species_present.append(species)
    #print("%s " % (record.id))
    species_present.sort()
    if (species_present == species_list) :
        #print("Alignment of length %i" % record.get_alignment_length())
        #print(species_present)
        AlignIO.write(record, output_handle, "maf")
        
output_handle.close()
input_handle.close()

print("All alignments containing all five species were written to the new file 5species.maf")









