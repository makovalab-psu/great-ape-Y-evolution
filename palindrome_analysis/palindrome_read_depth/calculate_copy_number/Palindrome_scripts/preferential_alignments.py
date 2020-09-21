#!/usr/bin/env python
"""
Given pairwise alignments from lastz's general format, sorted by some feature,
select alignments preferentially (description in the help text).
"""

from sys    import argv,stdin,stdout,stderr,exit
from math   import ceil


def usage(s=None):
	message = """

usage: cat alignments | preferential_alignments [options]
  --format=<list>      provide comma-separated list of the names of the
                       columns, in order; these must include the field names
                       listed in detail, below
  --format=auto        read column names from the first line of the input,
                       which must begin with a "#"
  --head=<number>      limit the number of input lines
  --progress=<number>  periodically report how many lines we've read

The input is alignments in lastz's general format.  Required columns are name2,
zstart2+, and end2+. The example below shows typical input but other columns
may be included.

#name1           zstart1+ end1+  name2     strand2 zstart2  end2     align1 align2
donkey.contig12  142287   142351 mule.chr1 -       9833952  9834016  CTG... CTG...
donkey.contig12  30291    30359  mule.chr1 +       8336423  8336491  CAG... CAG...
donkey.contig12  160824   160886 mule.chr1 -       25360782 25360844 TTT... TTT...
donkey.contig163 51135    51192  mule.chr1 -       11648921 11648978 ATT... ATT...
donkey.contig163 57665    57733  mule.chr1 -       11283997 11284065 GAA... GAA...
donkey.contig163 56549    56617  mule.chr1 +       18960023 18960091 TGA... TGA...
     ...

Output is in the same format.

The input alignments are assumed to be sorted in order of decreasing
preference.  We keep the first alignment, discard any alignments that intersect
that one (in sequence 2), then keep the first remaining alignment, and so on."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))

requiredColumns = ["name2","zstart2+","end2+"]


def main():
	global columnNames
	global debug

	# parse the command line

	columnNames     = None
	headLimit       = None
	reportProgress  = None
	debug           = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg == "--format=auto") or (arg == "--format=automatic"):
			columnNames = "automatic"
		elif (arg.startswith("--format=general:")):
			argVal = argVal.split(":",1)[1]
			columnNames = argVal.split(",")
		elif (arg.startswith("--format=")):
			columnNames = argVal.split(",")
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

	if (columnNames == None):
		usage("you must tell me the input column names")
	elif (columnNames == "automatic"):
		columnNames = None
	else:
		columnNames = column_names(columnNames)

	# process the alignments

	secondSequenceName = None
	coveredIntervals = []

	alignmentNumber = 0
	for a in read_alignments(stdin):
		alignmentNumber += 1
		if (reportProgress != None) and (alignmentNumber % reportProgress == 0):
			progressCount = commatize(alignmentNumber)
			print >>stderr, "progress: %s alignments read" % commatize(alignmentNumber)

		if (headLimit != None) and (alignmentNumber > headLimit):
			print >>stderr, "limit of %d alignments reached" % headLimit
			break

		if (secondSequenceName == None):
			secondSequenceName = a.name2
		else:
			assert (a.name2 == secondSequenceName), \
				   "inconsistent second sequence at line %d, %s and %s" \
				 % (a.lineNumber,secondSequenceName,a.name2)

		if (is_overlapped(coveredIntervals,a.start2,a.end2)):
			continue

		coveredIntervals += [(a.start2,a.end2)]
		print a.line


# read_alignments--

class Alignment: pass

def read_alignments(f):
	global columnNames
	global name2Column,start2Column,end2Column

	columnsNeeded = None

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
			columnNames = column_names(fields)
			print "#" + "\t".join(fields)
			continue

		assert (columnNames != None), \
		       "input column names are not provided within the file"

		if (columnsNeeded == None):
			columnsNeeded = 1 + max([columnNames[name] for name in columnNames])

			(name2Column,start2Column,end2Column) = requiredColumns

			name2Column  = columnNames[name2Column ]
			start2Column = columnNames[start2Column]
			end2Column   = columnNames[end2Column  ]

		fields = line.split()
		assert (len(fields) >= columnsNeeded), \
		      "not enough columns at line %d (%d, expected %d)" \
		    % (lineNumber,len(fields),columnsNeeded)

		a = Alignment()
		a.lineNumber = lineNumber
		a.line       = line
		a.name2      = fields[name2Column ]
		a.start2     = fields[start2Column]
		a.end2       = fields[end2Column  ]

		try:
			a.start2 = int(a.start2)
			a.end2   = int(a.end2)
			if (a.start2 >= a.end2): raise ValueError
		except ValueError:
			assert (False), "bad alignment (at line %d), second species start/end\n%s" \
			              % (lineNumber,line)

		yield a


# is_overlapped--

def is_overlapped(intervals,start,end):
	for (s,e) in intervals:
		if (e > start) and (s < end): return True
	return False


# column_names--

def column_names(names):
	columnNames = {}
	for (ix,name) in enumerate(names):
		actualName = name
		if (name not in requiredColumns): continue
		if (name in columnNames):
			usage("\"%s\" (or alias) appears more than once in --format" % actualName)
		columnNames[name] = ix
	for name in requiredColumns:
		if (name not in columnNames):
			usage("--format lacks required name \"%s\"" % name)
	return columnNames


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
