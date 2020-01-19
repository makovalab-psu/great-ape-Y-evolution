#!/usr/bin/env python
"""
Read alignments in maf format and output pairwise identity stats.
"""

from sys         import argv,stdin,stdout,stderr,exit
from math        import ceil
from string      import maketrans
from collections import OrderedDict 
from gzip        import open as gzip_open
import copy
from maf_reader  import maf_alignments

try:                from fasta_file import FastaFile
except ImportError: FastaFile = None


def usage(s=None):
	message = """

usage: cat maf | maf_to_pairwise_identity <species> [options]
  <species>                (required,cumulative) species name; at least two
                           species are required; stats will be reported for all
                           pairs of species
  --fasta:<species>=<file> (required) specify the sequence for each species
  --discard:weeds          discard any components that go beyond the end of the
                           sequence ('off in the weeds')
                           (by default we treat these as errors and halt)
  --events=<filename>      write the base-by-base alignment events vectors to
                           files; <filename> must contain "{ref}" and
                           "{other}", which will be replaced by species names
  --skip=<number>          skip the first so many maf blocks
  --head=<number>          limit the number of maf blocks read
  --progress=<number>      periodically report how many maf blocks we've
                           processed

The output is a tab-delimited table, like the one shown below. Each row is
stats for one species ("ref") aligned to another ("other"). refNonN is the
count of non-N bases in the reference fasta file. "m", "mm", "i" and "d" are
the number of positions in the reference that have a match, mismatch
(substitution), insertion, or deletion alignment event, respectively. Matches
take precedence over the other events, so if a position has any alignment block
for which that position is a match to the other species, it counts as a match,
and only as a match. Similarly, substituions have precedence over indels, and
insertions have precedence over deletions. Insertions are when the reference
has a nucleotide and the other sequence has a gap. Deletions are undercounted
since a delete-run of any length in the reference will only be counted in one
position.

"Aligned" is the percentage of non-N reference bases that have an alignment,
(m+mm+i)/nonNCount. "Identity" is the percentage of the aligned bases that are
matches, m/(m+mm+i). Note that deletions aren't involved in either calculation.
However, when the roles of "ref" and "other" are reversed, these become
insertions and contribute to those calculations.

The reported "identity" may be inflated due to the precedence given to matches.
For example, if one reference segment aligns to two segments of the other
species, and both have undergone 10% random substitutions, the reported identity
will be roughly 99% (not 90%).

  #ref   other  refNonN  m        mm     i      d   aligned identity
  apple  orange 23636355 15056518 401444 892123 120 69.17%  92.09%
  orange apple  25284285 18447869 543260 73134  552 75.40%  96.77%
   ...

The input is a maf file, as per the spec at
    http://genome.ucsc.edu/FAQ/FAQformat.html#format5

Note: This is designed to work at the scale of alignments of a small chromosome
(e.g human and chimp Y). For whole mammalian genomes the memory requirements
will probably be prohibitive."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global pairToData,nameToLookup
	global debug

	# parse the command line

	speciesNames   = []
	nameToFilename = {}
	discardWeeds   = False
	eventsFilename = None
	skipCount      = None
	headLimit      = None
	reportProgress = None
	debug          = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--fasta:")):
			(name,filename) = arg.split(":",1)[1].split("=",1)
			assert (name not in nameToFilename), \
			       "(in \"%s\") \"%s\" appears for the second time" % (arg,name)
			nameToFilename[name] = filename
		elif (arg in ["--discard:weeds","--discard=weeds"]):
			discardWeeds = True
		elif (arg.startswith("--events=")):
			eventsFilename = argVal
			assert ("{ref}" in eventsFilename) and ("{other}" in eventsFilename), \
			       "(in \"%s\") either or both {ref} and {other} are missing" % (eventsFilename)
		elif (arg.startswith("--skip=")):
			skipCount = int_with_unit(argVal)
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			speciesNames += [arg]

	speciesNames = list(OrderedDict.fromkeys(speciesNames)) # remove duplicates
	if (len(speciesNames) < 2):
		usage("you have to give me the names of at least two species")

	for name in speciesNames:
		assert (name in nameToFilename), "you didn't give me a fasta file for \"%s\"" % name
	for name in nameToFilename:
		assert (name in speciesNames), "I don't need a fasta file for \"%s\"" % name

	# open fasta files

	print >>stderr, "species names: %s" % ",".join(speciesNames)

	assert (FastaFile != None), \
		   "the FastaFile module could not be imported"

	nameToLookup = {}
	for name in nameToFilename:
		filename = nameToFilename[name]
		if (filename.endswith(".gz")) or (filename.endswith(".gzip")):
			nameToLookup[name] = FastaFile(gzip_open(filename,"rt"))
		else:
			nameToLookup[name] = FastaFile(filename)

	# process the alignments

	pairToData = {}		# base-by-base alignment events
	nameToLengths = {}	# species.contig lengths as reported in maf file

	skipsRemaining = skipCount
	mafBlockNumber = mafBlocksSkipped = 0
	for a in maf_alignments(stdin,yieldHeader=True,discardWeeds=discardWeeds):
		if (type(a) == str):  # this is a maf header
			continue

		if (skipsRemaining != None):
			mafBlocksSkipped += 1
			if (mafBlocksSkipped <= skipsRemaining): continue
			print >>stderr, "first %s maf blocks skipped" % commatize(skipCount)
			skipsRemaining = None

		mafBlockNumber += 1
		if (headLimit != None) and (mafBlockNumber > headLimit):
			print >>stderr, "limit of %s input maf blocks reached" % commatize(headLimit)
			break
		if (reportProgress != None):
			if (mafBlockNumber == 1) or (mafBlockNumber % reportProgress == 0):
				print >>stderr, "progress: maf block %s (line %s)" \
				              % (commatize(mafBlockNumber),commatize(a.lineNum))

		for (c1,c2) in pairwise_alignments(a.block,speciesNames):
			if ("maf" in debug):
				print_as_maf(stderr,c1,c2)

			# $$$ could sanity check that the maf file contains consistent contig lengths
			if (c1.ref not in nameToLengths): nameToLengths[c1.ref] = {}
			if (c1.contig not in nameToLengths[c1.ref]): nameToLengths[c1.ref][c1.contig] = c1.srcSize
			if (c2.ref not in nameToLengths): nameToLengths[c2.ref] = {}
			if (c2.contig not in nameToLengths[c2.ref]): nameToLengths[c2.ref][c2.contig] = c2.srcSize

			accumulate_base_by_base(c1,c2)

	if (reportProgress != None) and (headLimit == None):
		print >>stderr, "progress: %s maf blocks total" \
		              % (commatize(mafBlockNumber))

	# report stats

	refToNonNCount = {}  # counts non-Ns in a species

	print "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" \
	   % ("ref","other","refNonN","m","mm","i","d","aligned","identity")

	for (refName,otherName) in all_pairs(speciesNames,ordered=True):
		speciesPair = (refName,otherName)
		if (speciesPair not in pairToData):
			print "%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA" \
			    % (refName,otherName)
			continue

		if (refName not in refToNonNCount):
			refLookup = nameToLookup[refName]
			refContigSeen = set()
			refToNonNCount[refName] = 0
			for refContig in refLookup.keys():
				refContigSeen.add(refContig)
				refSeq = refLookup[refContig].upper()
				refToNonNCount[refName] += len([nuc for nuc in refSeq if (nuc != "N")])
				if    (refContig in nameToLengths[refName]) \
				  and (len(refSeq) != nameToLengths[refName][refContig]):
					assert (False), \
					       "maf input says len(%s.%s) is %s, but in fasta file it is %s" \
					     % (refName,refContig,
					        commatize(nameToLengths[refName][refContig]),
					        commatize(len(refSeq)))
			for refContig in nameToLengths[refName]:
				assert (refContig in refContigSeen), \
				       "maf input refers to %s.%s, but fasta file has no such sequence" \
				     % (refName,refContig)

		nonNCount = refToNonNCount[refName]

		refData = pairToData[speciesPair]
		eventToCount = {matchEvent:0,subEvent:0,insEvent:0,delEvent:0,}
		for refContig in refData:
			contigEventVector = refData[refContig]
			refPositions = contigEventVector.keys()
			refPositions.sort()
			for refPos in refPositions:
				eventToCount[contigEventVector[refPos]] += 1

		m  = eventToCount[matchEvent]
		mm = eventToCount[subEvent]
		i  = eventToCount[insEvent]
		d  = eventToCount[delEvent]

		print "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.2f%%\t%.2f%%" \
		   % (refName,otherName,nonNCount,m,mm,i,d,
		      100.0*(m+mm+i)/nonNCount,
		      100.0*m/(m+mm+i))

		# report events vector

		if (eventsFilename != None):
			filename = eventsFilename.replace("{ref}",refName)
			filename = filename.replace("{other}",otherName)
			print >>stderr, "writing events vector to \"%s\"" % filename
			if (filename.endswith(".gz")) or (filename.endswith(".gzip")):
				eventsF = gzip_open(filename,"wt")
			else:
				eventsF = file(filename,"wt")

			for refContig in refData:
				contigEventVector = refData[refContig]
				refPositions = contigEventVector.keys()
				refPositions.sort()
				for refPos in refPositions:
					print >>eventsF, "%s\t%d\t%s"\
					               % (refContig,refPos,eventToCh[contigEventVector[refPos]])
			eventsF.close()


# pairwise_alignments--
#	yield reference-centric pairwise blocks from a maf alignment block
#
# The first component in each reported pairwise block is always on the "+"
# strand. Note that each block is reported twice, once with one species as the
# first component, once with the other as first.

def pairwise_alignments(block,speciesNames):
	if (len(block) < 2): return

	# identify all the relevant components

	nameToComponents = {}
	speciesPresent = []

	for (componentIx,c) in enumerate(block):
		if (c.ref not in speciesNames): continue
		if (c.ref not in nameToComponents):
			nameToComponents[c.ref] = [componentIx]
			speciesPresent += [c.ref]
		else:
			nameToComponents[c.ref] += [componentIx]

	if (len(speciesPresent) < 2): return

	# yield the pairwise alignments

	componentToReverse = {}

	for (refName,otherName) in all_pairs(speciesPresent,ordered=True):
		for componentIx1 in nameToComponents[refName]:
			c1 = block[componentIx1]
			flipForStrand = False
			if (c1.strand == "-"):
				if (componentIx1 not in componentToReverse):
					componentToReverse[componentIx1] = reverse_of_component(c1)
				c1 = componentToReverse[componentIx1]
				flipForStrand = True

			for componentIx2 in nameToComponents[otherName]:
				c2 = block[componentIx2]
				if (flipForStrand):
					if (componentIx2 not in componentToReverse):
						componentToReverse[componentIx2] = reverse_of_component(c2)
					c2 = componentToReverse[componentIx2]

				yield (c1,c2)


# reverse_of_component--
#	reverse complement a component

def reverse_of_component(c):
	c = copy.copy(c)
	c.start  = c.srcSize - (c.start + c.length)
	c.strand = "+" if (c.strand == "-") else "-"
	c.nucs   = reverse_complement(c.nucs)
	return c


# accumulate_base_by_base--
#	accumulate each reference base into a base-by-base vector of alignment
#	events

matchEvent = 4	# match
subEvent   = 3	# substitution
insEvent   = 2	# inserted in ref
delEvent   = 1	# deleted from ref

eventToCh = {matchEvent:"M",subEvent:"S",insEvent:"I",delEvent:"D"}

def accumulate_base_by_base(cRef,cOther):
	global pairToData

	assert (cRef.strand == "+")

	refName     = cRef.ref
	refContig   = cRef.contig		# note that this can be None
	otherName   = cOther.ref
	otherContig = cOther.contig		# note that this can be None

	refNucs   = cRef.nucs.upper()
	otherNucs = cOther.nucs.upper()

	speciesPair = (refName,otherName)
	if (speciesPair not in pairToData): pairToData[speciesPair] = {}
	refData = pairToData[speciesPair]

	if (refContig not in refData): refData[refContig] = {}
	contigEventVector = refData[refContig]

	refPos = cRef.start
	for (refCh,otherCh) in zip(refNucs,otherNucs):
		refIsGap   = (refCh == "-")
		otherIsGap = (otherCh == "-")
		if (refIsGap) and (otherIsGap): continue
		if (refCh == "N"):              continue
		if (refCh == otherCh):          event = matchEvent
		elif (refIsGap == otherIsGap):  event = subEvent
		elif (otherIsGap):              event = insEvent
		else:                           event = delEvent

		if (refPos not in contigEventVector):
			contigEventVector[refPos] = event
		elif (event > contigEventVector[refPos]):
			contigEventVector[refPos] = event
		
		if (refCh != "-"): refPos += 1


# print_as_maf--

def print_as_maf(f,c1,c2):
	fields1 = component_to_fields(c1)
	fields2 = component_to_fields(c2)

	fmt = "s"
	for (ix,field1) in enumerate(fields1):
		field2 = fields2[ix]
		if (ix == 0):
			fmt += " %%-%ds" % max(len(field1),len(field2))
		else:
			fmt += " %%%ds"  % max(len(field1),len(field2))
	print >>f
	print >>f, "a"
	print >>f, fmt % fields1
	print >>f, fmt % fields2


# component_to_fields--

def component_to_fields(c):
	fields = []
	fields += [c.ref+"."+c.contig if (c.contig != None) else c.ref]
	fields += [str(c.start)]
	fields += [str(c.length)]
	fields += [str(c.strand)]
	fields += [str(c.srcSize)]
	fields += [str(c.nucs)]
	return tuple(fields)


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


# reverse_complement--

complementMap = maketrans("ACGTSWRYMKBDHVNacgtswrymkbdhvn",
                          "TGCASWYRKMVHDBNtgcaswyrkmvhdbn")

def reverse_complement(nukes):
	return nukes[::-1].translate(complementMap)


# int_with_unit--
#	Parse a string as an integer, allowing unit suffixes

def int_with_unit(s):
	if (s.endswith("K")):
		multiplier = 1000
		s = s[:-1]
	elif (s.endswith("M")):
		multiplier = 1000 * 1000
		s = s[:-1]
	elif (s.endswith("G")):
		multiplier = 1000 * 1000 * 1000
		s = s[:-1]
	else:
		multiplier = 1

	try:               return          int(s)   * multiplier
	except ValueError: return int(ceil(float(s) * multiplier))


# commatize--
#	Convert a numeric string into one with commas.

def commatize(s):
	if (type(s) != str): s = str(s)
	(prefix,val,suffix) = ("",s,"")
	if (val.startswith("-")): (prefix,val) = ("-",val[1:])
	if ("." in val):
		(val,suffix) = val.split(".",1)
		suffix = "." + suffix

	try:    int(val)
	except: return s

	digits = len(val)
	if (digits > 3):
		leader = digits % 3
		chunks = []
		if (leader != 0):
			chunks += [val[:leader]]
		chunks += [val[ix:ix+3] for ix in xrange(leader,digits,3)]
		val = ",".join(chunks)

	return prefix + val + suffix


if __name__ == "__main__": main()
