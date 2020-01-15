#!/usr/bin/env python
"""
fasta_N_intervals-- report the runs of Ns in fasta sequences.
"""

from sys  import argv,stdin,stderr,exit
from math import ceil


def usage(s=None):
	message = """
usage: cat <fasta_file> | fasta_N_intervals [options]
  --skip=<length>      don't report N-intervals shorter than <length>
  --unnamed=<name>     name for any unnamed sequences (of course, this makes
                       sense only when the file contains a single unnamed
                       sequence)
  --complement         report the mathematical complement of N-intervals
  --progress=<number>  periodically report how many bp we've processed

Intervals are reported as origin zero half open."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse args

	unnamed        = None
	skipLength     = 0
	doComplement   = False
	reportProgress = None

	for arg in argv[1:]:
		if (arg.startswith("--unnamed=")):
			unnamed = arg.split("=",1)[1]
		elif (arg.startswith("--skip=")) or (arg.startswith("--boundary=")):
			skipLength = int(arg.split("=",1)[1])
		elif (arg in ["--complement","--opposite","--invert"]):
			doComplement = True
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(arg.split("=",1)[1])
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# hunt for Ns

	for (name,start,end) in runs_of_n(stdin,doComplement=doComplement,
	                                  reportProgress=reportProgress):
		if (end-start <= skipLength):
			continue
		if (name == None):
			if (unnamed != None): name = unnamed
			else:                 name = ""
		print "%s\t%s\t%s" % (name,start,end)


# runs_of_n--
#	read a fasta file and yield all runs of N

def runs_of_n(f,doComplement=False,reportProgress=None):
	name  = None
	pos   = 0
	start = None

	lineNum = 0
	for line in f:
		lineNum += 1
		line = line.strip()
		if (line.startswith(">")):
			if (reportProgress != None) and (name != None):
				print >>stderr, "%s:%d" % (name,pos)
			if (start != None):
				yield (name,start,pos)
			name  = line[1:].strip()
			pos   = 0
			start = None
			if (reportProgress != None):
				print >>stderr, "%s:%d" % (name,1)
				prevProgress = 0
			continue

		if (reportProgress != None) and (pos / reportProgress > prevProgress):
			print >>stderr, "%s:%d" % (name,pos+1)
			prevProgress = pos / reportProgress

		for ch in line:
			if (ch == " "): continue
			assert (ch in "ACGTacgtNn"), \
			       "non-ACGTN character in line %d" % lineNum

			inRun = (ch in "Nn")
			if (doComplement): inRun = not inRun

			if (not inRun):
				if (start != None):
					yield (name,start,pos)
					start = None
			elif (start == None):
				start = pos

			pos += 1

	if (start != None):
		yield (name,start,pos)


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

	try:               return            int(s)   * multiplier
	except ValueError: return int(ceil(float(s) * multiplier))


if __name__ == "__main__": main()
