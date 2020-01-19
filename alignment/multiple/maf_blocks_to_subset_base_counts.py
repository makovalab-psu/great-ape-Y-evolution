#!/usr/bin/env python
"""
Read alignments in maf format and, conceptually, partition blocks by the
set of species present. Collect and report species-specific stats (defined
below) within each subset.
"""

from sys        import argv,stdin,stdout,stderr,exit
from math       import ceil
from maf_reader import maf_alignments


def usage(s=None):
	message = """

usage: cat maf | maf_blocks_to_subset_base_counts <reference> [options]
  --species=<list>       (cumulative) comma-separated list of species of
                         interest; required only for some operations
  --discard:weeds        discard any components that go beyond the end of the
                         sequence ('off in the weeds')
                         (by default we treat these as errors and halt)
  --skip=<number>        skip the first so many maf blocks
  --head=<number>        limit the number of maf blocks read
  --progress=<number>    periodically report how many maf blocks we've
                         processed

Read alignments in maf format and, conceptually, partition blocks by the set
of species present. Collect and report species-specific stats (defined below)
within each subset.

Three counts are colected for each subset:
  - "One-to-one" counts bases in segments that are in blocks that have no
    duplicate in any species.
  - "Unique" counts bases in segments that are in blocks for which the species
    has no other segment in the same block. This *includes* the bases counted
    in "one-to-one".
  - "Duplicate" counts bases in segments for which the species does have other
    segments in the same block (i.e. it has an aligning segment in the same
    species).

Examples:

  This block ...                 contributes these counts
    species segment_length         species one-to-one unique duplicate
    human         27               human        0        0     27+30
    human         30               chimp        0       24       0
    chimp         24               bonobo       0        0    30+22+26
    bonobo        30               gorilla      0        0     20+28
    bonobo        22               orang        0       23       0
    bonobo        26
    bonobo        28
    gorilla       20
    gorilla       28
    orang         23

  This block ...                 contributes these counts
    species segment_length         species one-to-one unique duplicate
    human         23               human       23       23       0
    chimp         30               chimp       30       30       0
    bonobo        26               bonobo      26       26       0
    gorilla       21               gorilla     21       21       0
    orang         25               orang       25       25       0

The input is a maf file, as per the spec at
    http://genome.ucsc.edu/FAQ/FAQformat.html#format5"""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global debug

	# parse the command line

	speciesOfInterest = None
	discardWeeds      = False
	skipCount         = None
	headLimit         = None
	reportProgress    = None
	debug             = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--species=")):
			if (argVal.strip() == ""): usage("unrecognized option: %s" % arg)
			if (speciesOfInterest == None): speciesOfInterest = []
			for species in argVal.split(","):
				species = species.strip()
				if (species not in speciesOfInterest):
					speciesOfInterest += [species]
		elif (arg in ["--discard:weeds","--discard=weeds"]):
			discardWeeds = True
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
			usage("unrecognized option: %s" % arg)

	# process the alignments

	refSubsetSeen = set()
	refSubsetToOneToOne  = {}
	refSubsetToUnique    = {}
	refSubsetToDuplicate = {}

	skipsRemaining = skipCount
	mafBlockNumber = mafBlocksSkipped = 0
	for a in maf_alignments(stdin,discardWeeds=discardWeeds):
		if (skipsRemaining != None):
			mafBlocksSkipped += 1
			if (mafBlocksSkipped <= skipsRemaining): continue
			print >>stderr, "first %s maf blocks skipped" % commatize(skipCount)
			skipsRemaining = None

		mafBlockNumber += 1
		if (headLimit != None) and (mafBlockNumber > headLimit):
			print >>stderr, "limit of %s input maf blocks reached" % commatize(headLimit)
			break
		if (reportProgress != None) and (mafBlockNumber % reportProgress == 0):
			print >>stderr, "progress: maf block %s (line %s)" \
			              % (commatize(mafBlockNumber),commatize(a.lineNum))

		# determine the ref subset in this block

		refSubset = set([c.ref for c in a.block])
		if (speciesOfInterest == None):
			refSubset = list(refSubset)
			refSubset.sort()
		else:
			refSubset = [ref for ref in speciesOfInterest if (ref in refSubset)]
			if (refSubset == []): continue
		refSubset = tuple(refSubset)

		# determine whether this block is one-to-one or not

		refToCount = {ref:0 for ref in refSubset}
		for c in a.block:
			if (c.ref not in refSubset): continue
			refToCount[c.ref] += 1

		isOneToOne = (max([refToCount[ref] for ref in refSubset]) == 1)

		# incorporate this block into the counts

		if (refSubset not in refSubsetSeen):
			refSubsetSeen.add(refSubset)
			refSubsetToOneToOne [refSubset] = {species:0 for species in refSubset}
			refSubsetToUnique   [refSubset] = {species:0 for species in refSubset}
			refSubsetToDuplicate[refSubset] = {species:0 for species in refSubset}

		for c in a.block:
			if (c.ref not in refSubset): continue
			if (isOneToOne):
				refSubsetToOneToOne [refSubset][c.ref] += c.length
				refSubsetToUnique   [refSubset][c.ref] += c.length
			elif (refToCount[c.ref] == 1):
				refSubsetToUnique   [refSubset][c.ref] += c.length
			else: # if (refToCount[c.ref] >=2 1):
				refSubsetToDuplicate[refSubset][c.ref] += c.length

	if (reportProgress != None) and (headLimit == None):
		print >>stderr, "progress: %s maf blocks total" \
		              % (commatize(mafBlockNumber))

	# report stats

	if (speciesOfInterest != None):
		speciesSeen = speciesOfInterest
	else:
		speciesSeen = set()
		for subset in refSubsetToOneToOne:
			speciesSeen = speciesSeen.union(set(subset))
		speciesSeen = list(speciesSeen)

	speciesToIndex = {}
	for (speciesIx,species) in enumerate(speciesSeen):
		speciesToIndex[species] = speciesIx

	refSubsets = []
	for refSubset in refSubsetToOneToOne:
		subsetIndices = []
		for (speciesIx,species) in enumerate(speciesSeen):
			if (species not in refSubset): continue
			subsetIndices += [speciesIx]
		refSubsets += [(-len(refSubset),tuple(subsetIndices),refSubset)]
	refSubsets.sort()

	print "#class\t%s" % ("\t".join(speciesSeen))

	for (_,_,refSubset) in refSubsets:
		counts = ["NA"] * len(speciesSeen)
		for species in refSubset:
			counts[speciesToIndex[species]] = refSubsetToOneToOne[refSubset][species]
		print "onetoone\t%s" % "\t".join(map(str,counts))

		counts = ["NA"] * len(speciesSeen)
		for species in refSubset:
			counts[speciesToIndex[species]] = refSubsetToUnique[refSubset][species]
		print "unique\t%s" % "\t".join(map(str,counts))

		counts = ["NA"] * len(speciesSeen)
		for species in refSubset:
			counts[speciesToIndex[species]] = refSubsetToDuplicate[refSubset][species]
		print "duplicate\t%s" % "\t".join(map(str,counts))


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
