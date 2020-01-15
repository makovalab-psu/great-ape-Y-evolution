#!/usr/bin/env python
"""
Read a fasta file and report the name and nucleotide content (counts) of all
the sequences.
"""

import sys

def main():

	assert (len(sys.argv) == 1), "give me no arguments"

	# read the fasta file

	print "#%s" % "\t".join(["name","length","nonN","A","C","G","T","a","c","g","t","N","n","other"])
	for (name,sequence) in fasta_sequences(sys.stdin):
		nucToCounts = {nuc:0 for nuc in "ACGTacgtNn" }
		nucToCounts[None] = 0
		for nuc in sequence:
			if (nuc not in "ACGTacgtNn"): nuc = None
			nucToCounts[nuc] += 1

		nonN = sum([nucToCounts[nuc] for nuc in "ACGTacgt"])
		print "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d" \
		    % (name,len(sequence),nonN,
		       nucToCounts["A"],nucToCounts["C"],nucToCounts["G"],nucToCounts["T"],
		       nucToCounts["a"],nucToCounts["c"],nucToCounts["g"],nucToCounts["t"],
		       nucToCounts["N"],nucToCounts["n"],nucToCounts[None])


# fasta_sequences--
#	Read the fasta sequences from a file

def fasta_sequences(f):
	seqName = None
	seqNucs = None

	for line in f:
		line = line.strip()

		if (line.startswith(">")):
			if (seqName != None):
				yield (seqName,"".join(seqNucs))
			seqName = sequence_name(line)
			seqNucs = []
		elif (seqName == None):
			assert (False), "first sequence has no header"
		else:
			seqNucs += [line]

	if (seqName != None):
		yield (seqName,"".join(seqNucs))


# sequence_name--
#	Extract the sequence name from a fasta header.
#	$$$ this may need to be improved $$$

def sequence_name(s):
	s = s[1:].strip()
	if (s == ""): return ""
	else:         return s.split()[0]


if __name__ == "__main__": main()

