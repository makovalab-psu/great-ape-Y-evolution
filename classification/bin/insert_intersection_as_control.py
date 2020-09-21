import argparse
from time import sleep

parser = argparse.ArgumentParser()
parser.add_argument("intersection_file", help="Provide intersection bed file")
parser.add_argument("windows_file", help="Provide bed file with all windows")
parser.add_argument("assembly_name", help="Provide reference name")
args = parser.parse_args()
print(args.intersection_file)
print(args.windows_file)
print(args.assembly_name)

with open(args.intersection_file) as f:
    controls = map(str.strip, list(f))
    controls = [c.replace('\t', '_') for c in controls]

output=args.assembly_name+"/"+"window_definition."+args.assembly_name+".tab"

file = open(output, "w")

with open(args.windows_file) as fp:  
	line = fp.readline()
	while line:
		key="_".join(line.split('\t')[0:3]) #only tab separated chr start end are used as a key
		key=str(key.rstrip())
		#print(key)
		#Should this window be a Control window?
		if key in controls:
			#this is a control window
			array=line.split('\t')
			file.write(array[3]) #unique window name first
			file.write("\t")
			file.write(str(int(array[1])+1)) #bed is 0-based, but we want 1-based coordinates
			file.write("\t")
			file.write("\t".join(array[2:3]))
			file.write("\tControl\t")
			file.write("\t".join(array[4:]))
		else:
			#this is a regular window
			array=line.split('\t')
			file.write(array[3]) #unique winow name first
			file.write("\t")
			file.write(str(int(array[1])+1)) #bed is 0-based, but we want 1-based coordinates
			file.write("\t")
			file.write("\t".join(array[2:3]))
			file.write("\t")
			file.write(array[3])
			file.write("\t")
			file.write("\t".join(array[4:]))
		line = fp.readline()


file.close()
