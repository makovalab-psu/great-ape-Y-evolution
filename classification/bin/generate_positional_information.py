import argparse
parser = argparse.ArgumentParser()
parser.add_argument("file", help="Provide chrY_MAPPABILITY.bed file")
parser.add_argument("assembly_name", help="Provide assembly name")
parser.add_argument("reference_length", help="Provide assembly length")

args = parser.parse_args()
print(args.file)
print(args.assembly_name)
print(args.reference_length)

mappability=args.file #'/nfs/brubeck.bx.psu.edu/scratch5/rahul/GTEx/REF/GEM-MAPPABILITY/chrY_MAPPABILITY.bed'
mapFH=open(mappability,"rU")
MAP=mapFH.read()
#def mappability_to_vector(bedfile):
map_val = []
for line in MAP.splitlines():
	col=line.split("\t")
	map_val+=[float(col[4])]*(int(col[2])-int(col[1]))

mapFH.close()

output=args.assembly_name+"/"+"chrY_MAPPABILITY_"+args.assembly_name+"_bypos.txt"
file = open(output, "w")
file.write("Mappability\n")
for index in range(len(map_val)):
    file.write(str(map_val[index])+"\n")

#because mappability uses k-mers, last n positions remain unfilled, thus should be filled with zeroes
gap_to_be_padded=int(args.reference_length)-len(map_val)
print("gap_to_be_padded: " + str(gap_to_be_padded))

for index in range(gap_to_be_padded):
    file.write(str("0\n"))

file.close()