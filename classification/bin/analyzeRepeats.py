import argparse

parser = argparse.ArgumentParser()
parser.add_argument("repeatmasker_file", help="Provide reference out file for Y chromosome")
parser.add_argument("assembly_name", help="Provide assembly name")
parser.add_argument("reference_length", help="Provide assembly length")
args = parser.parse_args()
print(args.repeatmasker_file)
print(args.assembly_name)
print(args.reference_length)

repeatmasker_output=args.assembly_name+"/chrY_"+args.assembly_name+"_REPEAT_MASKED.txt" 

repeatmasker=args.repeatmasker_file #'/nfs/brubeck.bx.psu.edu/scratch5/rahul/GTEx/REF/RM/hg38.sorted_chrY.bed'
rmasFH=open(repeatmasker,"rU")
MAS=rmasFH.read()
#def mappability_to_vector(bedfile):
mas_val = [1]*int(args.reference_length)
for line in MAS.splitlines():
	col=line.split("\t")
	mas_val[int(col[1]):int(col[2])]=[0]*(int(col[2])-int(col[1]))

rmasFH.close()


file = open(repeatmasker_output, "w")
file.write("RepeatMasked\n")
for index in range(len(mas_val)):
	file.write(str(mas_val[index])+ "\n")

file.close()

