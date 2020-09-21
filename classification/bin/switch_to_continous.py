import argparse
from time import sleep
import os

parser = argparse.ArgumentParser()
parser.add_argument("windows_file", help="Provide windows as .tab file")
parser.add_argument("offset", help="Provide the offset for switching assembly to a continous reference")
args = parser.parse_args()
#print(args.windows_file)
#print(args.offset)

with open(args.windows_file) as f: 
	for line in f:
		array=line.rstrip().split('\t')
		chr_name=array[0]
		start=int(array[1])
		end=int(array[2])
		annotation=array[3]
		window_name=chr_name + "_" + str(start) + "_" + str(end)

		offset_lookup = "egrep -B 1 " + chr_name + " " + args.offset + " | head -n 1 | cut -f3" #continous reference and scaffold assembly have different offsets
		offset=int(os.popen(offset_lookup).read()) #store the offset
		
		if (annotation=="Control"):
			#This is a control window
			print(window_name + "\t" + str(start+offset) + "\t" + str(end+offset) + "\t"+ "Control" + "\t"+ str("\t".join(array[4:])))
		else:
			#This is a normal window
			print(window_name + "\t" + str(start+offset) + "\t" + str(end+offset) + "\t"+ window_name + "\t"+ str("\t".join(array[4:])))
f.close()
