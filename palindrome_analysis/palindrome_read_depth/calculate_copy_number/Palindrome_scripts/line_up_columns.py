#!/usr/bin/env python

import sys


def parse_columns(s):
	if (".." in s):
		(colLo,colHi) = s.split("..",1)
		(colLo,colHi) = (int(colLo),int(colHi))
		return [col-1 for col in xrange(colLo,colHi+1)]
	else:
		cols = s.split(",")
		return [int(col)-1 for col in cols]


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


spacer       = " "
inSeparator  = None
outSeparator = None
passComments = []
bunchExtras  = None  # column number after which we no longer align
rightJustify = []    # column numbers to right-justify
columnWidth  = {}    # map from column number to width
commas       = []    # column numbers to commatize values (numeric string is assumed)
doubleSep    = []    # column numbers to give a double output separator before
barSep       = []    # column numbers to give a special output separator before

for arg in sys.argv[1:]:
	if ("=" in arg):
		argVal = arg.split("=",1)[1]

	if (arg.startswith("--spacer=")):
		spacer = argVal
	elif (arg.startswith("--separator=")):
		inSeparator = argVal
	elif (arg == "--tabs"):
		inSeparator = "\t"
	elif (arg == "--bars"):
		outSeparator = " |"
	elif (arg == "--passcomments"):
		passComments += ["#"]
	elif (arg.startswith("--passcomments=")):
		passComments += [argVal]
	elif (arg.startswith("--bunch=")):
		bunchExtras = int(argVal)
	elif (arg.startswith("--rightjustify=")) or (arg.startswith("--right=")):
		rightJustify += parse_columns(argVal)
	elif (arg.startswith("--width=")):
		for spec in argVal.split(","):
			assert (":" in spec), "can't understand: %s" % arg
			(col,w) = spec.split(":",1)
			(col,w) = (int(col)-1,int(w))
			assert (col not in columnWidth)
			columnWidth[col] = w
	elif (arg.startswith("--commas=")) or (arg.startswith("--commatize=")):
		commas += parse_columns(argVal)
	elif (arg.startswith("--double=")):
		doubleSep += parse_columns(argVal)
	elif (arg.startswith("--bar=")):
		barSep += parse_columns(argVal)
	elif (arg.startswith("--")):
		assert (False), "unknown argument: %s" % arg
	else:
		assert (False), "unknown argument: %s" % arg

if (passComments == []): passComments = None

colToWidth = {}
for col in columnWidth:  colToWidth[col] = columnWidth[col]

lines = []
for line in sys.stdin:
	line  =  line.strip()
	lines += [line]

	if (line == ""): continue

	if (passComments != None):
		isComment = False
		for prefix in passComments:
			if (line.startswith(prefix)):
				isComment = True
				break
		if (isComment):
			continue

	if (inSeparator == None): columns = line.split()
	else:                     columns = line.split(inSeparator)
	if (bunchExtras != None) and (len(columns) > bunchExtras):
		columns = columns[:bunchExtras] + [" ".join(columns[bunchExtras:])]

	for col in range(len(columns)):
		if (col in columnWidth): continue
		if (col in commas): w = len(commatize(columns[col]))
		else:               w = len(columns[col])
		if (col not in colToWidth): colToWidth[col] = w
		else:                       colToWidth[col] = max(w,colToWidth[col])

for line in lines:
	if (line == ""):
		print
		continue

	if (passComments != None):
		isComment = False
		for prefix in passComments:
			if (line.startswith(prefix)):
				isComment = True
				break
		if (isComment):
			print line
			continue

	if (inSeparator == None): columns = line.split()
	else:                     columns = line.split(inSeparator)
	if (bunchExtras != None) and (len(columns) > bunchExtras):
		columns = columns[:bunchExtras] + [" ".join(columns[bunchExtras:])]

	for col in range(len(columns)):
		if (col in rightJustify): fmt = "%*s"
		else:                     fmt = "%-*s"
		val = columns[col]
		if (col in commas): val = commatize(val)
		val = fmt % (colToWidth[col],val)
		if (outSeparator != None):
			val += outSeparator
			if (col+1 in barSep):    val += "|"
			if (col+1 in doubleSep): val += outSeparator.strip()
		else:
			if (col+1 in barSep): val += " |"
		columns[col] = val
	line = spacer.join(columns)
	print line.rstrip()

