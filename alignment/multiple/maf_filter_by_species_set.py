#!/usr/bin/env python
"""
Read alignments in maf format and output those blocks that have a specified
set of species, and no other species.
"""

from sys              import argv,stdin,stdout,stderr,exit
from math             import ceil
from maf_to_lzgeneral import maf_alignments


def usage(s=None):
	message = """

usage: cat maf | maf_filter_by_species_set <reference> [options]
  --species=<list>       (cumulative) comma-separated list of species of
                         interest
  --discard:duplicates   discard any *blocks* that contain more than one
                         component (a duplicate) for any of the species of
                         interest
  --discard:weeds        discard any components that go beyond the end of the
                         sequence ('off in the weeds')
                         (by default we treat these as errors and halt)
  --skip=<number>        skip the first so many maf blocks
  --head=<number>        limit the number of maf blocks read
  --progress=<number>    periodically report how many maf blocks we've
                         processed

Read alignments in maf format and output those blocks that have a specified
set of species, and no other species.  The block that are output are not
modified.

The input is a maf file, as per the spec at
    http://genome.ucsc.edu/FAQ/FAQformat.html#format5"""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global debug

	# parse the command line

	speciesOfInterest = None
	discardDuplicates = False
	discardWeeds      = False
	skipCount         = None
	headLimit         = None
	reportProgress    = None
	debug             = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--species=")):
			if (speciesOfInterest == None): speciesOfInterest = set()
			for species in argVal.split(","):
				species = species.strip()
				speciesOfInterest.add(species)
		elif (arg in ["--discard:duplicates","--discard=duplicates"]):
			discardDuplicates = True
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

	if (speciesOfInterest == None):
		usage("a species list is required")

	# process the alignments

	nothingHasPrinted = True

	skipsRemaining = skipCount
	mafBlockNumber = mafBlocksSkipped = 0
	mafBlocksKept = mafBlocksDiscarded = 0
	for a in maf_alignments(stdin,saveLines=True,discardWeeds=discardWeeds):
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
			print >>stderr, "progress: maf block %s (line %s), %s kept, %s discarded" \
			              % (commatize(mafBlockNumber),commatize(a.lineNum),
			                 commatize(mafBlocksKept),commatize(mafBlocksDiscarded))

		refs = set([c.ref for c in a.block])
		if (refs != speciesOfInterest):
			mafBlocksDiscarded += 1
			continue

		if (discardDuplicates):
			refs = set()
			hasDuplicate = False
			for c in a.block:
				if (c.ref in refs):
					hasDuplicate = True
					break
				refs.add(c.ref)
			if (hasDuplicate):
				mafBlocksDiscarded += 1
				continue

		mafBlocksKept += 1
		if (nothingHasPrinted): nothingHasPrinted = False
		else:                   print
		print "a"
		for c in a.block:
			print c.line

	if (reportProgress != None) and (headLimit == None):
		print >>stderr, "progress: %s maf blocks total" \
		              % (commatize(mafBlockNumber))

	if (mafBlocksKept + mafBlocksDiscarded == 0):
		print >>stderr, "no maf blocks kept, none discarded"
	else:
		print >>stderr, "%s maf blocks kept (%.2f%%), %s discarded" \
		              % (commatize(mafBlocksKept),
		                 100.0*mafBlocksKept/(mafBlocksKept+mafBlocksDiscarded),
		                 commatize(mafBlocksDiscarded))


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
