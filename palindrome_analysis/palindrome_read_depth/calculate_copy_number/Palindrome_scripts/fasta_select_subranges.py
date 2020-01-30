#!/usr/bin/env python
"""
Pick a specfied set of subintervals from fasta sequences.

	cat fasta_file | fasta_select_subranges interval_specs_file > fasta_file

The interval_specs_file looks like this (below).  Columns 2 and 3 form an
interval, origin-zero half-open (this can be overridden on the command line).
A fourth column can contain strand ("+" or "-").  Positions are always counted
along the forward strand.

	FX4JPGW02F85TR 0  31
	FX4JPGW02F6RB5 48 115
	FX4JPGW02HX4PB 43 101
	FX4JPGW02HLCJG 38 197
	 ...

"""

from sys    import argv,stdin,stderr,exit
from string import maketrans

def main():

	# parse the command line

	origin             = "zero"
	wrap               = None
	intervalsFile      = None
	intervalsHaveNames = False
	maskChar           = None
	maskNot            = False   # true => output masking NOT specified
	appendToName       = True
	copyOthers         = False

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--origin=")):
			origin = argVal
			if (origin == "0"): origin = "zero"
			if (origin == "1"): origin = "one"
			assert (origin in ["zero","one"]), \
			      "unknown argument: %s=%s" % (arg,val)
		elif (arg == "--intervalshavenames") or (arg == "--withnames"):
			intervalsHaveNames = True
		elif (arg.startswith("--wrap=")):
			wrap = int(argVal)
			assert (wrap > 0)
		elif (arg == "--mask"):
			(maskChar,maskNot) = ("X",False)
		elif (arg == "--masknot"):
			(maskChar,maskNot) = ("X",True)
		elif (arg.startswith("--mask=")):
			(maskChar,maskNot) = (argVal,False)
			assert (len(maskChar) == 1)
		elif (arg.startswith("--masknot=")):
			(maskChar,maskNot) = (argVal,True)
			assert (len(maskChar) == 1)
		elif (arg == "--names=original"):
			appendToName = False
		elif (arg == "--copy=others"):
			copyOthers = True
		elif (arg.startswith("--")):
			assert (False), "unknown argument: %s" % arg
		elif (intervalsFile == None):
			intervalsFile = arg
		else:
			assert (False), "unknown argument: %s" % arg

	assert (intervalsFile != None), \
	       "you have to tell me the intervals you're interested in"

	# read the intervals

	f = file(intervalsFile,"rt")

	nameToIntervals = {}
	haveStrands     = False

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == "") or (line.startswith("#")): continue

		fields = line.split()
		if (intervalsHaveNames):
			assert (len(fields) >= 4), \
			      "not enough fields (line %s): %s" % (lineNumber,line)
		else:
			assert (len(fields) >= 3), \
			      "not enough fields (line %s): %s" % (lineNumber,line)

		try:
			name  = fields[0]
			start = int(fields[1])
			end   = int(fields[2])
			if (origin == "one"): start -= 1
			if (start < 0):    raise ValueError
			if (start >= end): raise ValueError
			if (intervalsHaveNames):
				intervalName = fields[3]
		except ValueError:
			assert (False), \
			      "bad line (line %s): %s" % (lineNumber,line)

		strand = "+"
		if (len(fields) >= 4) and (fields[3] in "+-"):
			strand      = fields[3]
			haveStrands = True

		if (name not in nameToIntervals): nameToIntervals[name] = []
		if (intervalsHaveNames):
			nameToIntervals[name] += [(start,end,strand,intervalName)]
		else:
			nameToIntervals[name] += [(start,end,strand)]

	f.close()

	# read the fasta file and output intervals of interest

	for (name,sequence) in fasta_sequences(stdin):
		if (name not in nameToIntervals):
			if (copyOthers):
				print ">%s" %  name
				if (wrap != None): write_fasta_wrapped(sequence,wrap)
				else:              print sequence
			continue

		seqLen = len(sequence)

		for info in nameToIntervals[name]:
			if (intervalsHaveNames):
				(start,end,strand,intervalName) = info
			else:
				(start,end,strand) = info

			assert (end <= seqLen), \
			      "interval beyond end of sequence, %s %d %d (length is %d)" \
			    % (name,start,end,seqLen)

			if (maskChar == None):
				subseq = sequence[start:end]
			elif (maskNot):
				subseq = (maskChar * start) \
				       + sequence[start:end] \
				       + (maskChar * (seqLen-end))
			else:
				subseq = sequence[:start] \
				       + (maskChar * (end-start)) \
				       + sequence[end:]

			if (strand == "-"): subseq = reverse_complement(subseq)

			extra = ""
			if (haveStrands): extra = strand

			if (not intervalsHaveNames):
				if (appendToName): intervalName = "%s_%d_%d%s" % (name,start,end,extra)
				else:              intervalName = name

			if (origin == "one"): start += 1
			print ">%s" % intervalName
			if (wrap != None): write_fasta_wrapped(subseq,wrap)
			else:              print subseq


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


# reverse_complement

complementMap = maketrans("ACGTacgt","TGCAtgca")

def reverse_complement(nukes):
	return nukes[::-1].translate(complementMap)


# write_fasta_wrapped

def write_fasta_wrapped(seq,perLine=100):
	for i in range(0,len(seq),perLine):
		print seq[i:i+perLine]


if __name__ == "__main__": main()

