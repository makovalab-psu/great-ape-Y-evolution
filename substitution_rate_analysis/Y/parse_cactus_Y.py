import bx.align.maf
import bx.align.tools
from Bio import AlignIO

########################################################
#RESTRICT THE ALIGNMENT TO REGIONS WITH ALL FIVE SPECIES
input_handle = open("alignment.Anc0_centric.20191126.pruned.maf", "r")
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

####################################################
#RESTRICT THE ALIGNMENT TO X-DEGENERATE REGIONS ONLY
XDG_hg38_Y_coordinates = [
   ('2781766', '3049682'),
   ('6748711', '7604183'),
   ('11749732', '13983906'),
   ('14058734', '15874593'),
   ('15905215', '16159393'),
   ('16425966', '17455476'),
   ('18870335', '20054272'),
   ('20351234', '21335775'),
   ('26311846', '26660611')
]

input_handle = open("chrY.ancestor.pruned.5species.maf", "rU") #open("alignment.hg_Y_centric.20191126.maf", "rU")
output_handle= open("XDEG.maf","w")
align=bx.align.maf.Reader(input_handle)
parsed_MAF=bx.align.maf.Writer(output_handle)

print("Finished reading the .maf file.")

for msa in align :
    #print("###")
    msa_start=msa.get_component_by_src('hg_Y.chrY').get_forward_strand_start()
    msa_end=msa.get_component_by_src('hg_Y.chrY').get_forward_strand_end()

    #print(msa_start)
    #print(msa_end)
    #for each alignment block, check if it's in the XDEG region
    for region in XDG_hg38_Y_coordinates:
        region_start=region[0]
        region_end=region[1]

        #print(region_start)
        #print(region_end)

        if ((int(msa_start)) >= int(region_start)) and ((int(msa_end)) <= int(region_end)): 
            parsed_MAF.write(msa)
parsed_MAF.close()

print("All alignments where human coordinates intersect XDEG region were written to the new file XDEG.maf")









