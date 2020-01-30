#!/usr/bin/env python
"""
Read a fasta file and report the name and lengths of all the sequences.

This version doesn't use the bx library.
"""

import sys

def main():

	lengthFirst = False
	assert (len(sys.argv) in [1,2]), "give me no arguments"
	if (len(sys.argv) == 2):
		assert (sys.argv[1] == "--lengthfirst"), "unknown argument: %s" % sys.argv[1]
		lengthFirst = True

	# read the fasta file

	for (name,sequence) in fasta_sequences(sys.stdin):
		if (lengthFirst): print "%d %s" % (len(sequence),name)
		else:             print "%s %d" % (name,len(sequence))


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

