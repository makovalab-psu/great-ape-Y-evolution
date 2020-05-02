import argparse
from Bio import SeqIO
import pybedtools
import numpy

parser = argparse.ArgumentParser()
parser.add_argument("file", help="Provide continous reference fasta file")
parser.add_argument("assembly_name", help="Provide assembly name")

args = parser.parse_args()
print(args.file)
print(args.assembly_name)

output_file=args.assembly_name + "/chrY_GC_Content" + args.assembly_name + ".txt"

def chunks(seq, window):
	seqlen = len(seq)
	for i in range(0,seqlen,window):
		if (i+window>seqlen):
			j = seqlen 
		else:
			j = i+window
		if i-250 < 0:
			l=0
		else:
			l=i-250
		if j+250 > seqlen:
			r=seqlen
		else:
			r=j+250
		yield seq[l:r],i+1,j
		if j==seqlen:
			break


output = open(output_file, "w") #open for writing
output.write("GC\n")

refFH=open(args.file,"rU") #'/galaxy/home/rxv923/refs/Hg/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
for chr_record in SeqIO.parse(refFH, "fasta"): #should be a single sequence
		chr = chr_record.seq
		GC=[]
		for win in chunks(chr, 1):
			N = win[0].count("N")+win[0].count("n")
			A = win[0].count("A")+win[0].count("a")
			T = win[0].count("T")+win[0].count("t")
			G = win[0].count("G")+win[0].count("g")
			C = win[0].count("C")+win[0].count("c")
			#id=str(win[1])+"-"+str(win[2])
			GC_value=((G+C)/float(A+T+G+C+0.0000001))*100
			GC+=[GC_value]
			output.write(str(GC_value)+ "\n")
		

refFH.close()
output.close()