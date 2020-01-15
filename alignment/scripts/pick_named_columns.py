#!/usr/bin/env python
"""
Pick columns, by name, from a whitespace-delimited file
"""

from sys import argv,stdin,stderr,exit


def usage(s=None):
	message = """
usage: cat tabular file | pick_named_columns <column_names> [options]
  --except                 keep all columns *except* those listed
  --header                 produce header with columns labeled appropriately
                           (default is to output no header)
  --copyall                copy the entire input line following the picked
                           columns
  --separator=<separator>  specify separator for input and output fields
                           (default is tabs)
  --extra                  include a separator at the end of the line
                           (default is separator between fields only)
  --head=<number>          limit the number of input lines
  --progress=<number>      periodically report how many lines we've read"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

	# parse the command line

	columnsOfInterest  = []
	exceptThoseColumns = False
	outputHeader       = False
	copyAll            = False
	fieldSeparator     = "\t"
	extraSeparator     = False
	headLimit          = None
	headLimitQuiet     = False
	reportProgress     = None
	debug              = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg == "--except"):
			exceptThoseColumns = True
		elif (arg == "--header"):
			outputHeader = True
		elif (arg == "--copyall"):
			copyAll = True
		elif (arg.startswith("--separator=")) or (arg.startswith("--sep=")):
			fieldSeparator = argVal
			if   (fieldSeparator == "tab"):        fieldSeparator = "\t"
			elif (fieldSeparator == "space"):      fieldSeparator = " "
			elif (fieldSeparator == "none"):       fieldSeparator = ""
			elif (fieldSeparator == "whitespace"): fieldSeparator = None
		elif (arg == "--whitespace"):
			fieldSeparator = None
		elif (arg == "--tab") or (arg == "--tabs"):
			fieldSeparator = "\t"
		elif (arg == "--tab+") or (arg == "--tabs+"):
			fieldSeparator = "\t"
			extraSeparator = True
		elif (arg == "--space") or (arg == "--spaces"):
			fieldSeparator = " "
		elif (arg == "--space+") or (arg == "--spaces+"):
			fieldSeparator = " "
			extraSeparator = True
		elif (arg == "--extra"):
			extraSeparator = True
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
			headLimitQuiet = False
		elif (arg.startswith("--head:quiet=")):
			headLimit = int_with_unit(argVal)
			headLimitQuiet = True
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			columnsOfInterest += [arg]

	if (columnsOfInterest == []):
		usage("you must choose at least one column")

	if (fieldSeparator != None): fieldJoiner = fieldSeparator
	else:                        fieldJoiner = "\t"

	# process the lines

	nameToColumn = None
	numColumns   = None

	lineNumber = 0
	for line in stdin:
		lineNumber += 1
		if (headLimit != None) and (lineNumber > headLimit):
			if (not headLimitQuiet):
				print >>stderr, "limit of %s lines reached" % (commatize(headLimit))
			break
		if (reportProgress != None):
			if (lineNumber == 1) or (lineNumber % reportProgress == 0):
				print >>stderr, "reading line %s" % commatize(str(lineNumber))

		#line = line.strip()
		if (line == ""):
			continue
		if (line.startswith("#")):
			if (nameToColumn != None): continue
			fields = line[1:].strip().split(fieldSeparator)
			numColumns = len(fields)
			nameToColumn = {}
			for (colIx,name) in enumerate(fields):
				assert (name not in nameToColumn), \
				       "multiple columns are named \"%s\" (%d and %d)" \
				     % (name,nameToColumn[name]+1,colIx+1)
				nameToColumn[name] = colIx
				if ("header" in debug):
					print >>stderr, "column %d: \"%s\"" % (colIx,name)
			for name in columnsOfInterest:
				assert (name in nameToColumn), \
				       "input contains no column named \"%s\"" % name
			if (exceptThoseColumns):
				columnsOfDisinterest = columnsOfInterest
				columnsOfInterest = []
				for name in fields:
					if (name in columnsOfDisinterest): continue
					columnsOfInterest += [name]
				assert (columnsOfInterest != []), \
				       "for --except, after removing your columns there are none left"
			if (outputHeader):
				header = list(columnsOfInterest)
				if (copyAll): header += fields
				if (extraSeparator): header += [""]
				print "#%s" % fieldJoiner.join(header)
			continue

		assert (nameToColumn != None), \
		       "first line is not column headers (has to begin with #)"
		fields = line.split(fieldSeparator)
		assert (len(fields) == numColumns), \
		      "inconsistent number of fields at line %d (%d, expected %d)" \
		    % (lineNumber,len(fields),numColumns)

		outLine = []
		for name in columnsOfInterest:
			colIx = nameToColumn[name]
			outLine += [fields[colIx]]

		if (copyAll): outLine += fields
		if (extraSeparator): outLine += [""]
		print fieldJoiner.join(outLine)


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
