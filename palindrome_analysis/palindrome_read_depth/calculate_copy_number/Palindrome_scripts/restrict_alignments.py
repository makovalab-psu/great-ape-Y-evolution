#!/usr/bin/env python
"""
Given pairwise alignments from lastz's general format, restrict them to
specified intervals.
"""

from sys    import argv,stdin,stdout,stderr,exit
from math   import ceil


def usage(s=None):
	message = """

usage: cat alignments | restrict_alignments [options]
  --restriction1=<filename>  file containing restriction intervals for the
                             first species (see note below)
  --restriction2=<filename>  file containing restriction intervals for the
                             second species (see note below)
  --onebyone                 restriction intervals are applied one-by-one to
                             the corresponding line of the input alignments
                             (by default, restriction intervals apply
                             collectively to all input alignments)
  --format=<list>            provide comma-separated list of the names of the
                             columns, in order; these must include the field
                             names listed in detail, below
  --format=auto              read column names from the first line of the
                             input, which must begin with a "#"
  --minincomingid=<value>    discard incoming alignments with identity% lower
                             than specified; value should be in 0..100 range
  --minincomingcon=<value>   discard incoming alignments with continuity% lower
                             than specified; value should be in 0..100 range
  --minoutgoingid=<value>    discard any resulting alignments with identity%
                             lower than specified; value should be in 0..100
                             range;  this automatically enables --recompute:id
  --minoutgoingcon=<value>   discard any resulting alignments with continuity%
                             lower than specified; value should be in 0..100
                             range;  this automatically enables --recompute:con
  --minoutalignment=<length> discard any resulting alignments that have too
                             few alignment columns
  --recompute:id             compute the identity% of the resulting alignments,
                             replacing that value the column named "id"
  --recompute:con            compute the continuity% of the resulting
                             alignments, replacing that value the column named
                             "con"
  --recompute:cov1           compute the cov1 of the resulting alignments;
                             the incoming alignments must have cov1 in the form
                             <length1>/<size1>
  --recompute:cov2           compute the cov2 of the resulting alignments;
                             the incoming alignments must have cov2 in the form
                             <length2>/<size2>
  --compute:matchstats       compute the number of matches, mismatches, and gap
                             columns in the resulting alignments, inserting
                             the result in a new column just after the column
                             named "id"; the format is a colon-delimited list;
  --compute:matchstats+      same as compute:matchstats, except that gap counts
                             are split per-species, and the stats are inserted
                             as four new columns instead of one
  --compute:diff             compute the diff field, showing the differences
                             between align1 and align2, inserting the result in
                             a new column just after the column named "align1"
  --head=<number>            limit the number of input lines
  --progress=<number>        periodically report how many lines we've read

The input is alignments in lastz's general format.  Required columns are name1,
zstart1, end1, name2, strand2, zstart2+, end2+, align1 (or text1), and align2
(or text2).  The example below shows typical input but other columns may be
included, and there is no pre-sorting requirement.

    #name1    zstart1  end1     name2            strand2 zstart2+ end2+  align1 align2
    mule.chr1 9833952  9834016  donkey.contig12  -       142287   142351 CTG... CTG...
    mule.chr1 8336423  8336491  donkey.contig12  +       30291    30359  CAG... CAG...
    mule.chr1 25360782 25360844 donkey.contig12  -       160824   160886 TTT... TTT...
    mule.chr1 11648921 11648978 donkey.contig163 -       51135    51192  ATT... ATT...
    mule.chr1 11283997 11284065 donkey.contig163 -       57665    57733  GAA... GAA...
    mule.chr1 18960023 18960091 donkey.contig163 +       56549    56617  TGA... TGA...
     ...

Output is in the same format as input.

The restriction interval files indicate the intervals that are allowed;
alignments or parts of alignments outside these intervals are discarded.  The
files contain one interval per line, in the form chrom start end.  The interval
coordinates are origin-zero, half open.  There is no pre-sorting requirement."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))

requiredColumns = ["name1","zstart1","end1",
                   "name2","strand2","zstart2+","end2+",
                   "align1","align2"]
columnAliases   = {"s"     : "strand2",
                   "s2"    : "strand2",
                   "text1" : "align1",
                   "text2" : "align2",
                   "id%"   : "id",
                   "con%"  : "con",
                   "cov1"  : "cov1",
                   "cov2"  : "cov2"}


def main():
	global recomputedColumns,computedColumns,specialColumns,columnNames
	global debug

	# parse the command line

	restrictionFilename1 = None
	restrictionFilename2 = None
	oneByOne             = False
	columnNames          = None
	minIncomingId        = None
	minIncomingCon       = None
	minOutgoingId        = None
	minOutgoingCon       = None
	minOutgoingLength    = 0
	recomputedColumns    = []
	computedColumns      = []
	specialColumns       = []
	headLimit            = None
	reportProgress       = None
	debug                = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--restriction1=")) or (arg.startswith("--restrict1=")):
			restrictionFilename1 = argVal
		elif (arg.startswith("--restriction2=")) or (arg.startswith("--restrict2=")):
			restrictionFilename2 = argVal
		elif (arg == "--onebyone") or (arg == "--one-by-one"):
			oneByOne = True
		elif (arg == "--format=auto") or (arg == "--format=automatic"):
			columnNames = "automatic"
		elif (arg.startswith("--format=general:")):
			argVal = argVal.split(":",1)[1]
			columnNames = argVal.split(",")
		elif (arg.startswith("--format=")):
			columnNames = argVal.split(",")
		elif (arg.startswith("--minincomingid=")):
			if (argVal.endswith("%")): argVal = argVal[:-1]
			minIncomingId = float(argVal)
			assert 50 <= minIncomingId < 100
			minIncomingId /= 100.0
		elif (arg.startswith("--minincomingcon=")):
			if (argVal.endswith("%")): argVal = argVal[:-1]
			minIncomingCon = float(argVal)
			assert 50 <= minIncomingCon < 100
			minIncomingCon /= 100.0
		elif (arg.startswith("--minoutgoingid=")):
			if (argVal.endswith("%")): argVal = argVal[:-1]
			minOutgoingId = float(argVal)
			assert 50 <= minOutgoingId < 100
			minOutgoingId /= 100.0
			if ("id" not in recomputedColumns): recomputedColumns += ["id"]
			if ("id" not in specialColumns):    specialColumns    += ["id"]
		elif (arg.startswith("--minoutgoingcon=")):
			if (argVal.endswith("%")): argVal = argVal[:-1]
			minOutgoingCon = float(argVal)
			assert 50 <= minOutgoingCon < 100
			minOutgoingCon /= 100.0
			if ("con"  not in recomputedColumns): recomputedColumns += ["con"]
			if ("con"  not in specialColumns):    specialColumns    += ["con"]
		elif (arg.startswith("--minoutalignment=")) or (arg.startswith("--minoutgoingalignment=")):
			minOutgoingLength = int(argVal)
		elif (arg.startswith("--recompute:")):
			for name in arg.split(":",1)[1].split(","):
				name = name.strip()
				if (name == "id"):
					if ("id"   not in recomputedColumns): recomputedColumns += ["id"]
					if ("id"   not in specialColumns):    specialColumns    += ["id"]
				elif (name == "con"):
					if ("con"  not in recomputedColumns): recomputedColumns += ["con"]
					if ("con"  not in specialColumns):    specialColumns    += ["con"]
				elif (name == "cov1"):
					if ("cov1" not in recomputedColumns): recomputedColumns += ["cov1"]
					if ("cov1" not in specialColumns):    specialColumns    += ["cov1"]
				elif (name == "cov2"):
					if ("cov2" not in recomputedColumns): recomputedColumns += ["cov2"]
					if ("cov2" not in specialColumns):    specialColumns    += ["cov2"]
				else:
					usage("I don't know how to recompute \"%s\" (in \"%s\"" % (name,arg))
		elif (arg.startswith("--compute:")):
			for name in arg.split(":",1)[1].split(","):
				name = name.strip()
				if (name == "matchstats"):
					if ("m:mm:g" not in computedColumns): computedColumns += ["m:mm:g"]
					if ("id"   not in specialColumns):    specialColumns  += ["id"]
				elif (name == "matchstats+"):
					if ("m"    not in computedColumns):   computedColumns += ["m"]
					if ("mm"   not in computedColumns):   computedColumns += ["mm"]
					if ("g1"   not in computedColumns):   computedColumns += ["g1"]
					if ("g2"   not in computedColumns):   computedColumns += ["g2"]
					if ("id"   not in specialColumns):    specialColumns  += ["id"]
				elif (name == "diff"):
					if ("diff" not in computedColumns):   computedColumns += ["diff"]
				else:
					usage("I don't know how to compute \"%s\" (in \"%s\"" % (name,arg))
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

	if (restrictionFilename1 == None) and (restrictionFilename2 == None):
		usage("you must provide at least one restriction file")

	if (columnNames == None):
		usage("you must tell me the input column names")
	elif (columnNames == "automatic"):
		columnNames = None
	else:
		columnNames = column_names(columnNames)

	if (oneByOne):
		if (minIncomingId != None):
			usage("--minincomingid cannot be used with --onebyone")
		if (minIncomingCon != None):
			usage("--minincomingcon cannot be used with --onebyone")

	# read the restriction files

	if (oneByOne):
		chrom1ToIntervals = None
		chrom2ToIntervals = None

		restrictionList1 = None
		if (restrictionFilename1 != None):
			restrictionList1 = read_intervals(restrictionFilename1,asList=True)

		restrictionList2 = None
		if (restrictionFilename2 != None):
			restrictionList2 = read_intervals(restrictionFilename2,asList=True)

		if (restrictionList1 != None) and (restrictionList2 != None):
			assert (len(restrictionList1) == len(restrictionList2)), \
			       ("for --onebyone, restriction lists must have the same length" \
			     +  ", but yours have %d and %d intervals") \
			     % (len(restrictionList1),len(restrictionList2))
	else:
		chrom1ToIntervals = None
		if (restrictionFilename1 != None):
			chrom1ToIntervals = read_intervals(restrictionFilename1)

		chrom2ToIntervals = None
		if (restrictionFilename2 != None):
			chrom2ToIntervals = read_intervals(restrictionFilename2)

	# process the alignments

	idColumn   = None
	conColumn  = None
	cov1Column = None
	cov2Column = None

	alignmentNumber = 0
	for a in read_alignments(stdin):
		alignmentNumber += 1
		if (reportProgress != None) and (alignmentNumber % reportProgress == 0):
			progressCount = commatize(alignmentNumber)
			print >>stderr, "progress: %s alignments read" % commatize(alignmentNumber)

		if (headLimit != None) and (alignmentNumber > headLimit):
			print >>stderr, "limit of %d alignments reached" % headLimit
			break

		# if we're in one-by-one mode, pop the next restriction interval off
		# the list

		if (oneByOne):
			if (restrictionList1 != None):
				assert (restrictionList1 != []), \
				       "we've run out of restriction intervals but there are more input alignments"
				(chrom,start,end) = restrictionList1.pop(0)
				chrom1ToIntervals = { chrom: [(start,end)] }
				if ("onebyone" in debug):
					print >>stderr, "%d: %s" % (alignmentNumber,chrom1ToIntervals)

			if (restrictionList2 != None):
				assert (restrictionList2 != []), \
				       "we've run out of restriction intervals but there are more input alignments"
				(chrom,start,end) = restrictionList2.pop(0)
				chrom2ToIntervals = { chrom: [(start,end)] }
				if ("onebyone" in debug):
					print >>stderr, "%d: %s" % (alignmentNumber,chrom2ToIntervals)

		# if we're enforcing an identity lower bound, make sure this alignment
		# is satisfactory;  note that we don't trust the incoming identity
		# field because (a) there might not be one, (b) we might not know what
		# column it is, and (c) it might be wrong

		if (minIncomingId != None) or (minIncomingCon != None):
			(m,mm,g1,g2) = match_mismatch_gap(a.align1,a.align2)

		if (minIncomingId != None):
			denom = m + mm
			if (denom == 0): continue
			if (m < minIncomingId * denom): continue

		if (minIncomingCon != None):
			denom = m + mm + g1 + g2
			if (denom == 0): continue
			if (m+mm < minIncomingCon * denom): continue

		# if either chromosome is completely disallowed, just discard this
		# alignment

		if (chrom1ToIntervals != None) and (a.name1 not in chrom1ToIntervals):
			continue                       # none of this first species chrom is allowed
		if (chrom2ToIntervals != None) and (a.name2 not in chrom2ToIntervals):
			continue                       # none of this second species chrom is allowed

		# split into sub-alignments along the first species

		if ("unrejected" in debug):
			print >>stderr, "%s %d %d vs %s%s %d %d" \
			                % (a.name1,a.start1,a.end1,
			                   a.name2,a.strand2,a.start2,a.end2)

		if (chrom1ToIntervals == None):
			subAlignments = [a]            # no restrictions for first species
		else:
			intervals = chrom1ToIntervals[a.name1]
			allowed = restriction_of(intervals,a.start1,a.end1,chrom=a.name1)
			if (allowed == []): continue   # none of start1,end1 is allowed
			subAlignments = [restrict_species1_to(a,s1,e1) for (s1,e1) in allowed]
			subAlignments = [b for b in subAlignments
			                   if  (b != None)
			                   and (len(b.align1) >= minOutgoingLength)]

		if (subAlignments == []): continue # no sub-alignments survived splitting

		# split into sub-alignments along the second species

		if (chrom2ToIntervals == None):
			alignments = subAlignments    # no restrictions for second species
		else:
			intervals = chrom2ToIntervals[a.name2]
			alignments = []
			for b in subAlignments:
				allowed = restriction_of(intervals,b.start2,b.end2,chrom=b.name2)
				if (allowed == []): continue  # none of this start2,end2 is allowed
				subSubAlignments = [restrict_species2_to(b,s2,e2) for (s2,e2) in allowed]
				alignments += [b for b in subSubAlignments
				                 if  (b != None)
				                 and (len(b.align1) >= minOutgoingLength)]

		if (alignments == []): continue   # no alignments left

		if ("id" in specialColumns) and (idColumn == None):
			idColumn = columnNames["id"]
		if ("con" in specialColumns) and (conColumn == None):
			conColumn = columnNames["con"]
		if ("cov1" in specialColumns) and (cov1Column == None):
			cov1Column = columnNames["cov1"]
		if ("cov2" in specialColumns) and (cov2Column == None):
			cov2Column = columnNames["cov2"]

		for b in alignments:
			if (b.align1 == "") or (b.align2 == ""):
				continue  # (don't print empty alignments)

			fields = b.line.split()

			fields[name1Column  ] =     b.name1
			fields[start1Column ] = str(b.start1)
			fields[end1Column   ] = str(b.end1)
			fields[name2Column  ] =     b.name2
			fields[strand2Column] =     b.strand2
			fields[start2Column ] = str(b.start2)
			fields[end2Column   ] = str(b.end2)
			fields[align1Column ] =     b.align1
			fields[align2Column ] =     b.align2

			if ("id"      in recomputedColumns) \
			or ("con"     in recomputedColumns) \
			or ("m:mm:g"  in computedColumns) \
			or ("m"       in computedColumns) \
			or ("mm"      in computedColumns) \
			or ("g"       in computedColumns) \
			or ("g1"      in computedColumns) \
			or ("g2"      in computedColumns):
				(m,mm,g1,g2) = match_mismatch_gap(b.align1,b.align2)

			if ("id" in recomputedColumns):
				denom = m + mm
				if (denom == 0):
					idValue = 0.0
					fields[idColumn] = "0.0"
				else:
					idValue = m / float(denom)
					fields[idColumn] = "%.1f" % (100 * idValue)

			if ("con" in recomputedColumns):
				denom = m + mm + g1 + g2
				if (denom == 0):
					conValue = 0.0
					fields[conColumn] = "0.0"
				else:
					conValue = (m+mm) / float(denom)
					fields[conColumn] = "%.1f" % (100.0 * conValue)

			if ("cov1" in recomputedColumns):
				len1 = b.end1 - b.start1
				fields[cov1Column] = "%d/%d" % (len1,a.size1)

			if ("cov2" in recomputedColumns):
				len2 = b.end2 - b.start2
				fields[cov2Column] = "%d/%d" % (len2,a.size2)

			if ("m:mm:g" in computedColumns):
				fields.insert(idColumn+1,"%d:%d:%d" % (m,mm,g1+g2))
			if ("g2" in computedColumns):
				fields.insert(idColumn+1,"%d" % g2)
			if ("g1" in computedColumns):
				fields.insert(idColumn+1,"%d" % g1)
			if ("g" in computedColumns):
				fields.insert(idColumn+1,"%d" % g1+g2)
			if ("mm" in computedColumns):
				fields.insert(idColumn+1,"%d" % mm)
			if ("m" in computedColumns):
				fields.insert(idColumn+1,"%d" % m)

			if (minOutgoingId != None) and (idValue < minOutgoingId):
				continue
			if (minOutgoingCon != None) and (conValue < minOutgoingCon):
				continue

			if ("diff" in computedColumns):
				insertColumn = align1Column
				if (idColumn != None) and (insertColumn > idColumn):
					insertColumn += 1
				diff = diff_text(b.align1,b.align2)
				fields.insert(insertColumn+1,diff)

			print "\t".join(fields)

	# if we're in one-by-one mode, make sure we've used up all the restriction
	# intervals

	if (oneByOne):
		if (restrictionList1 != None):
			assert (restrictionList1 == []), \
				   "we've run out of input alignments but there are restriction intervals left over"
		if (restrictionList2 != None):
			assert (restrictionList2 == []), \
				   "we've run out of input alignments but there are restriction intervals left over"


# read_alignments--

class Alignment: pass

def read_alignments(f):
	global columnNames
	global name1Column,start1Column,end1Column
	global name2Column,strand2Column,start2Column,end2Column
	global align1Column,align2Column

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
			idColumn = None
			if ("m:mm:g" in computedColumns):
				idColumn = columnNames["id"]
				fields.insert(idColumn+1,"m:mm:g")
			if ("g2" in computedColumns):
				idColumn = columnNames["id"]
				fields.insert(idColumn+1,"g2")
			if ("g1" in computedColumns):
				idColumn = columnNames["id"]
				fields.insert(idColumn+1,"g1")
			if ("g" in computedColumns):
				idColumn = columnNames["id"]
				fields.insert(idColumn+1,"g")
			if ("mm" in computedColumns):
				idColumn = columnNames["id"]
				fields.insert(idColumn+1,"mm")
			if ("m" in computedColumns):
				idColumn = columnNames["id"]
				fields.insert(idColumn+1,"m")
			if ("diff" in computedColumns):
				align1Column = requiredColumns[7]
				insertColumn = columnNames[align1Column]
				if (idColumn != None) and (insertColumn > idColumn):
					insertColumn += 1
				fields.insert(insertColumn+1,"diff")
			print "#" + "\t".join(fields)
			continue

		assert (columnNames != None), \
		       "input column names are not provided within the file"

		if (columnsNeeded == None):
			columnsNeeded = 1 + max([columnNames[name] for name in columnNames])

			(name1Column,start1Column,end1Column,
             name2Column,strand2Column,start2Column,end2Column,
             align1Column,align2Column) = requiredColumns

			name1Column   = columnNames[name1Column  ]
			start1Column  = columnNames[start1Column ]
			end1Column    = columnNames[end1Column   ]
			name2Column   = columnNames[name2Column  ]
			strand2Column = columnNames[strand2Column]
			start2Column  = columnNames[start2Column ]
			end2Column    = columnNames[end2Column   ]
			align1Column  = columnNames[align1Column ]
			align2Column  = columnNames[align2Column ]

		fields = line.split()
		assert (len(fields) >= columnsNeeded), \
		      "not enough columns at line %d (%d, expected %d)" \
		    % (lineNumber,len(fields),columnsNeeded)

		a = Alignment()
		a.lineNumber = lineNumber
		a.line       = line
		a.name1      = fields[name1Column  ]
		a.start1     = fields[start1Column ]
		a.end1       = fields[end1Column   ]
		a.name2      = fields[name2Column  ]
		a.strand2    = fields[strand2Column]
		a.start2     = fields[start2Column ]
		a.end2       = fields[end2Column   ]
		a.align1     = fields[align1Column ]
		a.align2     = fields[align2Column ]

		try:
			a.start1 = int(a.start1)
			a.end1   = int(a.end1)
			if (a.start1 >= a.end1): raise ValueError
		except ValueError:
			assert (False), "bad alignment (at line %d), first species start/end\n%s" \
			              % (lineNumber,line)

		try:
			a.start2 = int(a.start2)
			a.end2   = int(a.end2)
			if (a.start2 >= a.end2): raise ValueError
		except ValueError:
			assert (False), "bad alignment (at line %d), second species start/end\n%s" \
			              % (lineNumber,line)

		assert (a.strand2 in ["+","-"]), "bad alignment (at line %d), second species strand\n%s" \
		                               % (lineNumber,line)

		if ("cov1" in specialColumns):
			cov1 = fields[columnNames["cov1"]]
			(a.len1,a.size1) = cov1.split("/",1)
			(a.len1,a.size1) = (int(a.len1),int(a.size1))

		if ("cov2" in specialColumns):
			cov2 = fields[columnNames["cov2"]]
			(a.len2,a.size2) = cov2.split("/",1)
			(a.len2,a.size2) = (int(a.len2),int(a.size2))

		yield a


# diff_text--
#	Create the same field as lastz's format=general:diff field
#		matches        . (period)
#		transitions    : (colon)
#		transversions  X
#		gaps           - (hyphen).

nucPairToDiff = {}
for nuc1 in "ACGTN":
	nucPairToDiff[("-",nuc1)] = "-"
	nucPairToDiff[(nuc1,"-")] = "-"
	for nuc2 in "ACGTN":
		nucPairToDiff[(nuc1,nuc2)] = "X"
	nucPairToDiff[(nuc1,nuc1)] = "."
nucPairToDiff[("A","G")] = ":"
nucPairToDiff[("G","A")] = ":"
nucPairToDiff[("C","T")] = ":"
nucPairToDiff[("T","C")] = ":"

def diff_text(align1,align2):
	return "".join([nucPairToDiff[pair] for pair in zip(align1.upper(),align2.upper())])


# column_names--

def column_names(names):
	columnNames = {}
	for (ix,name) in enumerate(names):
		actualName = name
		if (name in columnAliases): name = columnAliases[name]
		if (name not in requiredColumns + specialColumns): continue
		if (name in columnNames):
			usage("\"%s\" (or alias) appears more than once in --format" % actualName)
		columnNames[name] = ix
	for name in requiredColumns:
		if (name not in columnNames):
			usage("--format lacks required name \"%s\"" % name)
	for name in specialColumns:
		if (name not in columnNames):
			usage("--format lacks to-be-computed name \"%s\"" % name)
	return columnNames


# read_intervals--

def read_intervals(filename,asList=False):

	# read the intervals

	if (asList): intervalsList    = []
	else:        chromToIntervals = {}

	f = file(filename,"rt")

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line.startswith("#")): continue

		fields = line.split()
		assert (len(fields) == 3), \
			  "wrong number of fields at line %d (%d, expected %d)" \
			% (lineNumber,len(fields),3)

		try:
			chrom  =     fields[0]
			start  = int(fields[1])
			end    = int(fields[2])
			if (start >= end): raise ValueError
		except ValueError:
			assert (False), "bad interval at line %d in %s\n%s" \
						  % (lineNumber,filename,line)
		if (asList):
			intervalsList += [(chrom,start,end)]
		else: 
			if (chrom not in chromToIntervals): chromToIntervals[chrom] = []
			chromToIntervals[chrom] += [(start,end)]

	f.close()

	# if we're to return a list, we're done

	if (asList):
		return intervalsList

	# merge overlaps

	for chrom in chromToIntervals:
		intervals = chromToIntervals[chrom]
		intervals.sort()

		mergedIntervals = []
		start = None
		for (s,e) in intervals:
			if (start == None):
				(start,end) = (s,e)
			elif (s > end):
				mergedIntervals += [(start,end)]
				(start,end) = (s,e)
			elif (e > end):
				end = e

		if (start != None):
			mergedIntervals += [(start,end)]

		chromToIntervals[chrom] = mergedIntervals

	return chromToIntervals


# restriction_of--
#	Determine what portions (if any) of an interval are allowed

def restriction_of(allowedIntervals,start,end,chrom=None):
	# we assume intervals have been sorted and do not overlap nor abut

	if ("restriction_of" in debug):
		if (chrom != None): print >>stderr, "restricting %s %d %d" % (chrom,start,end)
		else:               print >>stderr, "restricting %d %d"    % (start,end)

	# find the first interval that overlaps this one;  if there is no such
	# interval, we return an empty list
	# $$$ binary search would probably be quicker

	ixFirst = None
	for (ix,(s,e)) in enumerate(allowedIntervals):
		if (e <= start): continue
		if (s < end):    ixFirst = ix
		break
	if (ixFirst == None): return []   # no part of start,end overlaps any interval
	assert (s < end) and (e > start)

	restriction = []

	for (s,e) in allowedIntervals[ixFirst:]:
		if (s >= end): break
		restriction += [(max(start,s),min(end,e))]

	if ("restriction_of" in debug):
		for (s,e) in restriction:
			print >>stderr, "  allowing %d %d" % (s,e)

	return restriction


# restrict_species1_to--
#	Restrict an alignment to a subinterval in the first species

def restrict_species1_to(a,start1,end1):
	b = Alignment()
	b.lineNumber = a.lineNumber
	b.line       = a.line
	b.name1      = a.name1
	b.name2      = a.name2
	b.strand2    = a.strand2

	# if there is no actual restriction just return a copy of the alignment

	start1 = max(a.start1,start1)
	end1   = min(a.end1,  end1)

	if (start1 == a.start1) and (end1 == a.end1):
		b.start1 = a.start1
		b.end1   = a.end1
		b.align1 = a.align1
		b.start2 = a.start2
		b.end2   = a.end2
		b.align2 = a.align2
		return b

	# locate the alignment columns corresponding to this subinterval
	# nota bene: we intend to keep columns leftIx..rightIx, inclusive

	leftIx  = pos_to_column(a.start1,a.end1,"+",a.align1,start1)
	rightIx = pos_to_column(a.start1,a.end1,"+",a.align1,end1-1)

	if ("pos_to_column" in debug):
		print >>stderr, "%d %d to start %d" % (a.start1,a.end1,start1)
		print >>stderr, " %s" % a.align1
		print >>stderr, "%*s^ %d" % (1+leftIx, " ",leftIx)
		print >>stderr, "%*s^ %d" % (1+rightIx," ",rightIx)

	# if the sub-alignment ends in gaps in the other species, we need to
	# trim the sub-alignment to remove them;  this is an interative process,
	# trimming for one species than the other, until we no longer have any end
	# gaps or until the sub-alignment is empty

	while (True):
		(leftGap,rightGap) = count_end_gaps(a.align2[leftIx:rightIx+1])
		if (leftGap + rightGap == 0): break
		leftIx  += leftGap
		rightIx -= rightGap
		if (leftIx == rightIx): return None

		(leftGap,rightGap) = count_end_gaps(a.align1[leftIx:rightIx+1])
		if (leftGap + rightGap == 0): break
		leftIx  += leftGap
		rightIx -= rightGap
		if (leftIx == rightIx): return None

	# now that we know the sub-alignment's bounds, we can convert back to
	# sequence positions

	b.start1 = column_to_pos(a.start1,a.end1,"+",a.align1,leftIx)
	b.end1   = column_to_pos(a.start1,a.end1,"+",a.align1,rightIx) + 1

	b.start2 = column_to_pos(a.start2,a.end2,a.strand2,a.align2,leftIx)
	b.end2   = column_to_pos(a.start2,a.end2,a.strand2,a.align2,rightIx) + 1
	if (b.strand2 == "-"):
		(b.start2,b.end2) = (b.end2-1,b.start2+1)

	b.align1 = a.align1[leftIx:rightIx+1]
	b.align2 = a.align2[leftIx:rightIx+1]

	return b


# restrict_species2_to--
#	Restrict an alignment to a subinterval in the second species

def restrict_species2_to(a,start2,end2):
	b = Alignment()
	b.lineNumber = a.lineNumber
	b.line       = a.line
	b.name1      = a.name1
	b.name2      = a.name2
	b.strand2    = a.strand2

	# if there is no actual restriction just return a copy of the alignment

	start2 = max(a.start2,start2)
	end2   = min(a.end2,  end2)

	if (start2 == a.start2) and (end2 == a.end2):
		b.start1 = a.start1
		b.end1   = a.end1
		b.align1 = a.align1
		b.start2 = a.start2
		b.end2   = a.end2
		b.align2 = a.align2
		return b

	# locate the alignment columns corresponding to this subinterval
	# nota bene: we intend to keep columns leftIx..rightIx, inclusive

	if (b.strand2 == "+"):
		leftIx  = pos_to_column(a.start2,a.end2,"+",a.align2,start2)
		rightIx = pos_to_column(a.start2,a.end2,"+",a.align2,end2-1)
	else: # if (b.strand2 == "-"):
		rightIx = pos_to_column(a.start2,a.end2,a.strand2,a.align2,start2)
		leftIx  = pos_to_column(a.start2,a.end2,a.strand2,a.align2,end2-1)

	if ("pos_to_column" in debug):
		print >>stderr, "%d %d to start %d" % (a.start2,a.end2,start2)
		print >>stderr, " %s" % a.align2
		print >>stderr, "%*s^ %d" % (1+rightIx," ",rightIx)
		print >>stderr, "%*s^ %d" % (1+leftIx, " ",leftIx)

	# if the sub-alignment ends in gaps in the other species, we need to
	# trim the sub-alignment to remove them;  this is an interative process,
	# trimming for one species than the other, until we no longer have any end
	# gaps or until the sub-alignment is empty

	while (True):
		(leftGap,rightGap) = count_end_gaps(a.align1[leftIx:rightIx+1])
		if (leftGap + rightGap == 0): break
		leftIx  += leftGap
		rightIx -= rightGap
		if (rightIx == leftIx): return None

		(leftGap,rightGap) = count_end_gaps(a.align2[leftIx:rightIx+1])
		if (leftGap + rightGap == 0): break
		leftIx  += leftGap
		rightIx -= rightGap
		if (rightIx == leftIx): return None

	# now that we know the sub-alignment's bounds, we can convert back to
	# sequence positions

	b.start1 = column_to_pos(a.start1,a.end1,"+",a.align1,leftIx)
	b.end1   = column_to_pos(a.start1,a.end1,"+",a.align1,rightIx) + 1

	if (b.strand2 == "+"):
		b.start2 = column_to_pos(a.start2,a.end2,a.strand2,a.align2,leftIx)
		b.end2   = column_to_pos(a.start2,a.end2,a.strand2,a.align2,rightIx) + 1
	else: # if (b.strand2 == "-"):
		b.start2 = column_to_pos(a.start2,a.end2,a.strand2,a.align2,rightIx)
		b.end2   = column_to_pos(a.start2,a.end2,a.strand2,a.align2,leftIx) + 1

	b.align1 = a.align1[leftIx:rightIx+1]
	b.align2 = a.align2[leftIx:rightIx+1]

	return b


# pos_to_column--
#	Return the index of the alignment column corresponding to a sequence
#	position

def pos_to_column(start,end,strand,alignmentText,pos):
	if (strand == "+"):
		ntPos = start
		for (colIx,nt) in enumerate(alignmentText):
			if (nt == "-"): continue
			if (ntPos == pos): return colIx
			ntPos += 1
		assert (False), "internal error"
	elif (strand == "-"):
		ntPos = end - 1
		for (colIx,nt) in enumerate(alignmentText):
			if (nt == "-"): continue
			if (ntPos == pos): return colIx
			ntPos -= 1
		assert (False), "internal error"
	else:
		assert (False), "internal error, bad strand: \"%s\"" % strand


# column_to_pos--
#	Return the sequence position corresponding to an index of an alignment
#	column

def column_to_pos(start,end,strand,alignmentText,colIx):
	assert (0 <= colIx < len(alignmentText)), "internal error"
	if (strand == "+"):
		ntPos = start
		for (ix,nt) in enumerate(alignmentText[:colIx]):
			if (nt != "-"): ntPos += 1
	elif (strand == "-"):
		ntPos = end - 1
		for (ix,nt) in enumerate(alignmentText[:colIx]):
			if (nt != "-"): ntPos -= 1
	else:
		assert (False), "internal error, bad strand: \"%s\"" % strand

	return ntPos


# count_end_gaps--
#	Count the number of gap characters at either end of some alignment text

def count_end_gaps(alignmentText):
	leftGap = None
	for (ix,nt) in enumerate(alignmentText):
		if (nt != "-"):
			leftGap = ix
			break
	if (leftGap == None):
		return (len(alignmentText),0)

	if (leftGap == 0):
		for (rightGap,nt) in enumerate(alignmentText[::-1]):
			if (nt != "-"): break
	else:
		for (rightGap,nt) in enumerate(alignmentText[:leftGap-1:-1]):
			if (nt != "-"): break

	return (leftGap,rightGap)


# match_mismatch_gap--
#	count matches, mismatches, and gaps in aligned strings

def match_mismatch_gap(nucs1,nucs2):
	assert (len(nucs1) == len(nucs2))

	m = mm = g1 = g2 = 0

	for (ch1,ch2) in zip(nucs1.upper(),nucs2.upper()):
		if   (ch1 == "-"): g1 += 1
		elif (ch2 == "-"): g2 += 1
		elif (ch1 == ch2): m  += 1
		else:              mm += 1

	return (m,mm,g1,g2)


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
