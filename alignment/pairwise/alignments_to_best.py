#!/usr/bin/env python
"""
From an alignment, collect intervals that are "best" as measured by some
column.  An interval may be split, if only parts of that interval are best.
"""

from sys    import argv,stdin,stdout,stderr,exit
from math   import ceil


def usage(s=None):
	message = """

usage: cat alignments | alignments_to_best [options]
  --interval=<column names> the columns defining the interval, as a comma-
                            separated list
                            (default is "name1,zstart1,end1")
  --interval=sequence1      the interval is "name1,zstart1,end1"
  --interval=sequence2      the interval is "name2,zstart2+,end2+"
  --key=<column name>       (required) the column to order by
  --format=<list>           provide comma-separated list of the names of the
                            columns, in order;  these must include the field
                            names for the interval
  --format=auto             read column names from the first line of the input,
                            which must begin with a "#"
  --output=intervals        output only the best intervals
  --head=<number>           limit the number of input lines
  --progress=<number>       periodically report how many lines we've read

The input is alignments in lastz's general format.  Required columns are as
indicated above in the description of the --format=<list> option.  This is
typical input but other columns may be included, and there is no pre-sorting
requirement.

	#name1 zstart1  end1     id
	chrX   14852396 14852633 97.7
	chrX   14852704 14852809 95.9
	chrX   15135564 15135684 98.6
	chrX   15135670 15135811 95.9
	chrX   15135672 15135862 99.0
	chrX   15809077 15809128 96.9
	chrX   15848181 15848280 99.6
	chrX   15864022 15864112 96.8
	chrX   15978732 15978840 98.8
	 ...

Output is in the same format, but with three additional columns added on the
left: chrom, start, end.

Caveat: we expect the intervals to be from relatively short chromosomes/contigs
or overall to cover a relatively short interval on any chromosome/contig.  We
allocate a list of positions on the chromosome/contig being worked on, so if
the chromosome/contig is large so will be the memory used."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global columnNames,intervalColumns,keyColumn
	global debug,reportProgress

	# parse the command line

	intervalColumns = ["name1","zstart1","end1"]
	keyColumn       = None
	columnNames     = None
	outputWhat      = "full"
	headLimit       = None
	reportProgress  = None
	debug           = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg == "--interval=sequence1"):
			intervalColumns = ["name1","zstart1","end1"]
		elif (arg == "--interval=sequence2"):
			intervalColumns = ["name2","zstart2+","end2+"]
		elif (arg.startswith("--interval=")):
			intervalColumns = argVal.split(",")
			if (len(intervalColumns) != 3):
				usage("\"%s\" has the wrong number of names (3)" % arg)
		elif (arg.startswith("--key=")) or (arg.startswith("--order=")):
			keyColumn = argVal
		elif (arg == "--format=auto") or (arg == "--format=automatic"):
			columnNames = "automatic"
		elif (arg.startswith("--format=general:")):
			argVal = argVal.split(":",1)[1]
			columnNames = column_names(argVal.split(","))
		elif (arg.startswith("--format=")):
			columnNames = column_names(argVal.split(","))
		elif (arg == "--output=intervals"):
			outputWhat = "intervals"
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

	if (keyColumn == None):
		usage("you must provide a column to define the order")
	elif (keyColumn in intervalColumns):
		usage("\"%s\" is not a valid ordering" % keyColumn)
	elif (type(columnNames) == list) and (keyColumn not in columnNames):
		usage("ordering \"%s\" is not one of the columns" % keyColumn)

	if (columnNames == None):
		usage("you must tell me the input column names")
	elif (columnNames == "automatic"):
		columnNames = None

	# collect the alignments

	chroms            = []
	chromToAlignments = {}

	alignmentNumber = 0
	for a in read_alignments(stdin):
		alignmentNumber += 1
		if (reportProgress != None) and (alignmentNumber % reportProgress == 0):
			print >>stderr, "progress: %s alignments read" % commatize(alignmentNumber)

		if (headLimit != None) and (alignmentNumber > headLimit):
			print >>stderr, "limit of %d alignments reached" % headLimit
			break

		if (a.name not in chromToAlignments):
			chroms += [a.name]
			chromToAlignments[a.name] = []
		length = a.end - a.start
		chromToAlignments[a.name] += [(a.key,length,-a.lineNumber,a.start,a.end,a.line)]

	if (reportProgress != None):
		print >>stderr, "progress: %s alignments read, total" % commatize(alignmentNumber)

	# process chromosome-by-chromosome

	if (outputWhat == "intervals"):
		print "#%s" % " ".join(intervalColumns)
		for chrom in chroms:
			if (reportProgress != None):
				print >>stderr, "progress: processing alignments on %s" % chrom
			alignments = pick_best_covered(chromToAlignments[chrom])
			for (start,end,_) in alignments:
				print chrom,start,end
	else: # if (outputWhat == "full"):
		print "#chrom start end %s" % (header)
		for chrom in chroms:
			if (reportProgress != None):
				print >>stderr, "progress: processing %s alignments on %s" \
				              % (commatize(len(chromToAlignments[chrom])),chrom)
			alignments = pick_best_covered(chromToAlignments[chrom])
			for (start,end,line) in alignments:
				print chrom,start,end,line


# read_alignments--

class Alignment: pass

def read_alignments(f):
	global header,columnNames

	(nameColumn,startColumn,endColumn) = intervalColumns

	if (columnNames != None):
		columnsNeeded = 1 + max([columnNames[name] for name in columnNames])

	header = None
	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line.startswith("#")):
			if (columnNames != None): continue
			header = line[1:]
			fields = line.split()
			fields[0] = fields[0][1:]
			if (type(columnNames) == list) and (keyColumn not in fields):
				assert False, "ordering \"%s\" is not one of the file's columns" % keyColumn
			columnNames = column_names(fields)
			columnsNeeded = 1 + max([columnNames[name] for name in columnNames])
			continue

		assert (columnNames != None), \
		       "input column names are not provided within the file"

		fields = line.split()
		assert (len(fields) >= columnsNeeded), \
		      "not enough columns at line %d (%d, expected %d)" \
		    % (lineNumber,len(fields),columnsNeeded)

		a = Alignment()

		a.lineNumber = lineNumber
		a.line       = line
		a.name       = fields[columnNames[nameColumn ]]
		a.start      = fields[columnNames[startColumn]]
		a.end        = fields[columnNames[endColumn  ]]
		a.key        = fields[columnNames[keyColumn  ]]

		try:
			a.start = int(a.start)
			a.end   = int(a.end)
			if (a.start >= a.end): raise ValueError
		except ValueError:
			assert (False), "bad alignment (at line %d), first start/end\n%s" % (lineNumber,line)

		try:
			a.key = int(a.key)
		except ValueError:
			try:
				a.key = float(a.key)
			except ValueError:
				pass

		yield a


# column_names--

def column_names(names):
	columnNames = {}
	for (ix,name) in enumerate(names):
		if (name in columnNames):
			usage("\"%s\" appears more than once in --format" % name)
		columnNames[name] = ix
	for name in intervalColumns:
		if (name not in columnNames):
			usage("--format lacks required name \"%s\"" % name)
	return columnNames


# pick_best_covered--
#	input alignments  are (key,key2,key3,start,end,line)
#	output alignments are (start,end,line)
#
# $$$ there may be better algorithms for this

def pick_best_covered(alignments):

	# make our own copy of the alignments and sort 'em

	if (reportProgress != None):
		print >>stderr, "progress: sorting alignments"

	alignments = list(alignments)
	alignments.sort()
	alignments.reverse()

	if (reportProgress != None):
		print >>stderr, "progress: alignments sorted"

	# scan the alignments from best to worst, occupying positions in the
	# chromosome as we go

	minPos = min([candStart for (_,_,_,candStart,candEnd,_) in alignments])
	maxPos = max([candEnd   for (_,_,_,candStart,candEnd,_) in alignments])
	occupied = [None] * (maxPos-minPos)

	if (reportProgress != None):
		print >>stderr, "progress: positional range is %d..%d" % (minPos,maxPos)

	for (candIx,(_,_,_,candStart,candEnd,_)) in enumerate(alignments):
		for pos in xrange(candStart,candEnd):
			if (occupied[pos-minPos] == None): occupied[pos-minPos] = candIx
		if (reportProgress != None) and ((1+candIx) % reportProgress == 0):
			print >>stderr, "progress: %s alignments processed" % commatize(1+candIx)

	if (reportProgress != None):
		print >>stderr, "progress: %s alignments processed, total" % commatize(len(alignments))

	# scan the occupied vector and output the corresponding intervals and
	# alignments

	bestIntervals = []

	prevEnd = None
	prevOcc = None

	for pos in xrange(maxPos-minPos):
		occ = occupied[pos]
		if (occ != prevOcc):
			if (prevOcc != None):
				(_,_,_,_,_,candLine) = alignments[prevOcc]
				bestIntervals += [(prevStart,minPos+pos,candLine)]
			prevOcc   = occ
			prevStart = minPos+pos

	if (prevOcc != None):
		(_,_,_,_,_,candLine) = alignments[prevOcc]
		bestIntervals += [(prevStart,maxPos,candLine)]

	return bestIntervals


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
