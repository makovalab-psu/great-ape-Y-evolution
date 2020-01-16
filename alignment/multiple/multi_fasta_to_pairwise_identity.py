#!/usr/bin/env python
"""
Read an alignment in multi-fasta format and output pairwise identity stats.
"""

from sys         import argv,stdin,stdout,stderr,exit
from math        import ceil
from collections import OrderedDict 


def usage(s=None):
	message = """

usage: cat maf | multi_fasta_to_pairwise_identity <species> [options]
  <species>  (cumulative) species name; if any species is specified, at least
             two must be; if no species is specified, all species in the input
             are used; stats will be reported for all pairs of species

The output is a tab-delimited table, like the one shown below. Each row is
stats for one species ("ref") aligned to another ("other"). nonN is the count
of non-N (and non-gap) bases in the seqeunce. "m", "mm", "i" and "d" are the
number of positions in the reference that have a match, mismatch (substitu-
tion), insertion, or deletion alignment event, respectively. Insertions are
when the sequence has a nucleotide and the other sequence has a gap. Deletions
are when the sequence has a gap and the other sequence has a nucleotide.

"Identity" is the percentage of the aligned bases that are matches, m/(m+mm+i).
Note that deletions aren't involved in the calculation. However, when the roles
of the two sequences are reversed, these become insertions and do contribute to
the calculation.

"m/(m+mm)" is the percentage of aligned bases that are matches.

   ...
  #ref   other  nonN    m       mm     i      d      identity m/(m+mm)
  apple  orange 4195502 3787652 261205 139924 468180 90.42%   93.55%
  orange apple  4522619 3787652 266772 468180 141063 83.75%   93.42%
   ...

The input is a fasta file representing a multiple alignment. Thus it may
contain gap characters. Counting gap characters, each sequence *must* have
exactly the same length.

Note: This is designed to work at the scale of alignments of a small chromosome
(e.g human and chimp Y). For whole mammalian genomes the memory requirements
will probably be prohibitive."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global pairToData,nameToLookup
	global debug

	# parse the command line

	speciesNames = []
	truncateTo   = None
	debug        = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--truncate=")):  # (unadvertised)
			truncateTo = int(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			speciesNames += [arg]

	if (speciesNames == []):
		speciesNames = None
	else:
		speciesNames = list(OrderedDict.fromkeys(speciesNames)) # remove duplicates
		if (len(speciesNames) < 2):
			usage("you have to give me the names of at least two species")

	# read the sequences

	if (speciesNames == None): discoveredNames = []
	nameToSequence = {}
	sequenceLength = None

	for (name,sequence) in fasta_sequences(stdin):
		if (speciesNames == None):
			discoveredNames += [name]
		elif (name not in speciesNames):
			continue

		assert (name not in nameToSequence), \
		       "input contains more than one sequence named \"%s\"" % name

		if (sequenceLength == None):
			sequenceLength = len(sequence)
		else:
			assert (len(sequence) == sequenceLength), \
			       "\"%s\" has a different length that the others (%d, expected %d)" \
			     % (name,len(sequence),sequenceLength)

		if (truncateTo != None): sequence = sequence[:truncateTo]
		nameToSequence[name] = sequence

	if (speciesNames == None):
		if (len(discoveredNames) == 0):
			usage("input contains no sequences, I need at least two")
		elif (len(discoveredNames) < 2):
			usage("input contains only one sequence, I need at least two")
		speciesNames = list(OrderedDict.fromkeys(discoveredNames))

	for name in speciesNames:
		assert (name in nameToSequence), \
			       "\"%s\" wasn't found in the input" % name

	# compute and report stats

	print "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" \
	   % ("ref","other","nonN","m","mm","i","d","identity","m/(m+mm)")

	for (refName,otherName) in all_pairs(speciesNames,ordered=True):
		refNucs   = nameToSequence[refName].upper()
		otherNucs = nameToSequence[otherName].upper()

		refLen = 0
		eventToCount = {"m":0,"mm":0,"i":0,"d":0,"n":0}
		for (refCh,otherCh) in zip(refNucs,otherNucs):
			refIsGap   = (refCh == "-")
			otherIsGap = (otherCh == "-")
			if (refIsGap) and (otherIsGap): continue
			if (not refIsGap):              refLen += 1
			if (refCh == "N"):              event = "n"
			elif (refCh == otherCh):        event = "m"
			elif (refIsGap == otherIsGap):  event = "mm"
			elif (otherIsGap):              event = "i"
			else:                           event = "d"
			eventToCount[event] += 1

		m      = eventToCount["m"]
		mm     = eventToCount["mm"]
		i      = eventToCount["i"]
		d      = eventToCount["d"]
		nCount = eventToCount["n"]

		column8 = "%.2f%%" % (100.0*m/(m+mm+i)) if (m+mm+1 != 0) else "NA"
		column9 = "%.2f%%" % (100.0*m/(m+mm))   if (m+mm   != 0) else "NA"

		print "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s" \
		   % (refName,otherName,refLen,m,mm,i,d,column8,column9)


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


# all pairs--
#	yield all pairs from a set
#
# if ordered is false, pairs (a,b) and (b,a) are considered the same, and only
# one is generated

def all_pairs(items,ordered=False):
	for ix in xrange(len(items)-1):
		for iy in xrange(ix+1,len(items)):
			yield (items[ix],items[iy])
			if (ordered): yield (items[iy],items[ix])


if __name__ == "__main__": main()
