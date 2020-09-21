import argparse
from time import sleep
import os

parser = argparse.ArgumentParser()
parser.add_argument("summary", help="Provide Summary as .tab file")
parser.add_argument("offset", help="Provide the offset for switching assembly to a continous reference")
args = parser.parse_args()
print(args.summary)
print(args.offset)

bump=50 #how many nucleotides around the contig/scaffold boundary should be affected?

with open(args.offset) as f: 
	f.readline() #skip the first line of the offset
	for line in f:
		array=line.rstrip().split('\t')
		#print(array)
		boundary=int(array[2])
		start=boundary-bump+1 #increment is accounting for the header
		end=boundary+bump+1 #increment accounting for the header
		offset_lookup = "awk 'FNR==" + str(start) + ",FNR==" + str(end) + "{$2=0};1' " +  args.summary #continous reference and scaffold assembly have different offsets
		#print(offset_lookup)
		os.popen(offset_lookup).read() #set the mappability for the boundary
f.close()
