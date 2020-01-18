#!/usr/bin/env python
"""
Read a tab-delimited file and compute expressions over the columns.

References:
  [1] Using eval() safely in python (lybniz2.sourceforge.net/safeeval.html)
"""

from __future__ import division  # make division non-truncating
from sys        import argv,stdin,stderr,exit
from math       import *
from random     import seed as random_seed,random as unit_random,randint,choice

md5_new = None
if (md5_new == None):
	try:                from hashlib import md5 as md5_new
	except ImportError: pass
if (md5_new == None):
	try:                from md5     import new as md5_new
	except ImportError: pass


def usage(s=None):
	message = """
usage: cat file | colrithmetic [options]
  <formula>               something to compute;  see description below
  --require:<criterion>   (cumulative) ignore any input lines that don't
                          satisfy the specified criterion;  this is a
                          statement that evaluates to true if the line should
                          be kept, or false if the line should be discarded
  --prohibit:<criterion>  (cumulative) ignore any input lines that satisfy
                          the specified criterion
  --header                the first row of the file is column headers, from 
                          which we can pull variable names
  --header=<char>         the header row begins with a particular charater,
                          which is not considered part of the first variable
                          name
  --remove=<c1,c2,...>    (cumulative) remove columns;  note that the columns
                          are listed like "c1", not "1".
  --names=<x:c1,y:c2,...> (cumulative) give names for incoming columns
  --divbyzero=<text>      text for divide by zero (default is "NA")
  --spaces                output separator is spaces
  --seed=<string>         set random seed; only useful if a formula uses the
                          random() function

Column names are used as variable names in expressions.  These can be given in
a header line.  By default, though, we assume no names are given, and we use
"c1", "c2", etc.

Each formula consists of a name, expression, destination column, and format.
These are all optional except for the expression.  So, at a minimum, the
formula could be like "c1+c3", in which case column 1 will be added to column
3, and the sum will be added as a new column at the end of the line, with
column header "c1+c3", and in whatever format is deemed approriate.

The full form of a formula is

	[<name>:][{<format>}][<dest>=]<expression>

<name> can be anything as long as it doesn't begin with "--".  Names are
intended for use as column headers.  By default, any input line that begins
with # is considered to be a list of column headers, and the <name> is written
as the new column header.  If <name> is absent, the <expression> is used as a
name.  Note that if your input has no column headers, you probably won't care
what name is assigned to your new column.

Formulas that want to refer to an incoming column with a header that doesn't
begin with a letter or underscore should add an underscore to the variable name
in the formula.  The program automatically introduces and underscore to the
name before assigning it a value.

<format> determines how values in the new column will be written to the output.
This is one of
	f           floating point, default precision
	f.<digits>  floating point, with specified precision
	%           floating point percentage, default precision
	%.<digits>  floating point percentage, with specified precision
	x           hexadecimal, with only the necessary digits
	x.<digits>  hexadecimal, with specified number of minimum digits

<dest> is a column name, e.g. "c5" (or if we have names from a header line,
then the name of the column).  The new value will be inserted in the output in
front of the fifth column.  This is always counted relative to the input
columns.  For example, even if the incoming c1 is to be deleted, a <dest> of c5
will always go in front of the value that was in the fifth column of input.  If
more than one <dest> is c5, they will all go in front of the fifth column (in
the order in which they appear on the command line).

<expression> is a restricted python expression.  The column names can be used
as input variables.  Output variables can be used in expressions, but only if
the variable is assigned in an earlier expression on the command line.  Also
note that if a column name is not a valid python variable name, using it in an
expession will not succeed."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

	# parse the command line

	formulas        = None
	requirements    = None
	prohibitions    = None
	useHeader       = False
	headerPrefix    = "#"
	removeColumns   = None
	nameColumns     = None
	divByZero       = "NA"
	separator       = "\t"
	debug           = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--require:")):
			if (requirements == None): requirements = []
			argVal = arg.split(":",1)[1]
			requirements += [argVal.strip()]
		elif (arg.startswith("--prohibit:")) or (arg.startswith("--reject:")):
			if (prohibitions == None): prohibitions = []
			argVal = arg.split(":",1)[1]
			prohibitions += [argVal.strip()]
		elif (arg.startswith("--remove=")):
			if (removeColumns == None): removeColumns = []
			removeColumns += argVal.split(",")
		elif (arg.startswith("--names=")) or (arg.startswith("--name=")):
			if (nameColumns == None): nameColumns = []
			nameColumns += argVal.split(",")
		elif (arg.startswith("~")):
			argVal = arg[1:]
			if (nameColumns == None): nameColumns = []
			nameColumns += argVal.split(",")
		elif (arg == "--header"):
			useHeader    = True
			headerPrefix = None
		elif (arg.startswith("--header=")):
			useHeader    = True
			headerPrefix = argVal
		elif (arg.startswith("--divbyzero=")):
			divByZero = argVal
		elif (arg == "--spaces"):
			separator = " "
		elif (arg.startswith("--seed=")):
			random_seed(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			if (formulas == None): formulas = []
			formulas += [arg]

	if (formulas == None) and (requirements == None) and (removeColumns == None):
		usage("you must give me some action to perform!")

	if (requirements == None): requirements = []
	if (prohibitions == None): prohibitions = []

	# parse the formulas

	formulaText = formulas if (formulas != None) else []
	formulas = []
	for arg in formulaText:
		if (useHeader): formula = parse_formula(arg,destIsColumn=False)
		else:           formula = parse_formula(arg,destIsColumn=True)
		if (formula == None):
			usage("unrecognized option: %s" % arg)
		(name,format,dest,expression) = formula
		if (name == None): name = expression
		formulas += [(name,format,dest,expression)]

	if (nameColumns != None):
		for spec in nameColumns:
			(name,colName) = spec.split(":")
			formulas = [(name,None,None,colName)] + formulas
			if (removeColumns == None): continue
			if (colName not in removeColumns): removeColumns += [colName]

	if ("formulas" in debug):
		print >>stderr, "=== formulas ==="
		for (name,format,dest,expression) in formulas:
			print >>stderr, "%s %s %s %s" % (name,format,dest,expression)

	if ("formulas" in debug) and (requirements != []):
		print >>stderr, "=== requirements ==="
		for criterion in requirements:
			print >>stderr, "%s" % (criterion)

	if ("formulas" in debug) and (prohibitions != []):
		print >>stderr, "=== prohibitions ==="
		for criterion in prohibitions:
			print >>stderr, "%s" % (criterion)

	# determine the minimum number of input columns we'll need

	maxCol = -1
	if (removeColumns == None): removeColumns = []

	if (useHeader):
		for (ix,name) in enumerate(removeColumns):
			assert (name not in removeColumns[:ix]), \
			       "%s appears more than once in --remove" % name

		neededColumns = None  # (we'll fill this in when we read the header)
	else:
		removalNames = removeColumns
		removeColumns = []
		for name in removalNames:
			col = parse_column(name)
			assert (col not in removeColumns), \
			       "%s appears more than once in --remove" % name
			removeColumns += [col]

		if (removeColumns != []):
			maxCol = max(removeColumns)

		neededColumns = []
		for (name,format,dest,expression) in formulas:
			if (dest != None) and (dest > maxCol): maxCol = dest
			for col in columns_needed(expression):
				if (col > maxCol): maxCol = col
				neededColumns += [col]

	columnsNeeded = maxCol+1

	# process the lines

	numColumns = None

	lineNum = 0
	for line in stdin:
		lineNum += 1
		line = line.strip()

		isHeader = False
		if (useHeader) and (lineNum == 1):
			isHeader = True
			if (headerPrefix != None):
				assert (line.startswith(headerPrefix)), \
				       "expected first line to begin with \"%s\"" % headerPrefix
				line = line[len(headerPrefix):].strip()
		elif (headerPrefix != None) and (line.startswith(headerPrefix)):
			isHeader = True
			line = line[len(headerPrefix):].strip()

		if (line == ""):
			if   (not isHeader):         print line
			elif (headerPrefix != None): print headerPrefix + line
			else:                        print line
			continue

		# split line into columns

		columns = line.split()

		if (useHeader) and (numColumns == None):
			# header required and this is the first line
			assert (len(columns) >= columnsNeeded), \
			       "not enough columns at line %d (%d, expected at least %d)" \
			     % (lineNum,len(columns),columnsNeeded)
			numColumns = len(columns)
			columnNames = columns
			# $$$ should each of columNames need to be a valid python variable?

			if ("header" in debug):
				print >>stderr, columnNames

			nameToColumn = {}
			for (col,name) in enumerate(columnNames):
				assert (name not in nameToColumn), \
					   "%s appears more than once in the header" % name
				nameToColumn[name] = col

			if (removeColumns != []):
				removalNames = removeColumns
				removeColumns = []
				for name in removalNames:
					assert (name in nameToColumn), \
						   "%s doesn't appear in the header" % name
					removeColumns += [nameToColumn[name]]

			for (ix,formula) in enumerate(formulas):
				(name,format,dest,expression) = formula
				if (dest == None): continue
				assert (dest in nameToColumn), \
					   "%s doesn't appear in the header" % dest
				dest = nameToColumn[dest]
				formulas[ix] = (name,format,dest,expression)

		elif (useHeader) and (isHeader):
			# header required and this is a later header line
			assert (len(columns) == numColumns), \
			       "inconsistent number of columns at line %d (%d, expected %d)" \
			     % (lineNum,len(columns),numColumns)
			assert (columns == columnNames), \
			       "inconsistent column headers at line %d" \
			     % lineNum
		elif (numColumns == None):
			# no header required and this is the first line
			assert (len(columns) >= columnsNeeded), \
			       "not enough columns at line %d (%d, expected at least %d)" \
			     % (lineNum,len(columns),columnsNeeded)
			numColumns = len(columns)
		else:
			# non-header
			assert (len(columns) == numColumns), \
			       "inconsistent number of columns at line %d (%d, expected %d)" \
			     % (lineNum,len(columns),numColumns)

		# create the expression dictionary

		if (not isHeader):
			if (useHeader):
				for (col,name) in enumerate(columnNames):
					if (name[0] != "_") and (not name[0].isalpha()):
						name = "_" + name
					safeDict[name] = to_value(columns[col])
				for (col,name) in enumerate(columnNames):
					colName = "c%d"%(1+col)
					if (colName not in columnNames):
						safeDict[colName] = to_value(columns[col])
			else:
				for col in neededColumns:
					colName = "c%d"%(1+col)
					safeDict[colName] = to_value(columns[col])

		# filter, by evaluating any requirements and prohibitions
		# $$$ we'd like to ignore failures here, and then have a second filtering
		# $$$ .. pass after we assign variables

		reject = False
		for criterion in requirements:
			if (isHeader): continue
			try:
				val = eval(criterion,{"__builtins__":None},safeDict)
			except NameError:
				print >>stderr, "failed to evaluate requirement \"%s\"" % criterion
				raise
			if (val == False):
				reject = True
				break
			if (val != True):
				assert (False), "requirement \"%s\" evaluates to \"%s\" on record at line %d\n%s" \
				              % (criterion,val,lineNumber,line)
		if (reject): continue

		for criterion in prohibitions:
			if (isHeader): continue
			try:
				val = eval(criterion,{"__builtins__":None},safeDict)
			except NameError:
				print >>stderr, "failed to evaluate prohibition \"%s\"" % criterion
				raise
			if (val == True):
				reject = True
				if ("evaluation" in debug): print >>stderr, "  (rejected)"
				break
			if (val != False):
				assert (False), "prohibition \"%s\" evaluates to \"%s\" on record at line %d\n%s" \
				              % (criterion,val,lineNumber,line)
		if (reject): continue

		# evaluate the formulas

		newText = []

		for (name,format,dest,expression) in formulas:
			if (isHeader):
				newText += [(dest,name)]
			else:
				try:
					val = eval(expression,{"__builtins__":None},safeDict)
					if (type(val) == list):
						valStr = [apply_format(format,v) for v in val]
						valStr = "[" + ",".join(valStr) + "]"
					elif (type(val) == tuple):
						valStr = (apply_format(format,v) for v in val)
						valStr = "(" + ",".join(valStr) + ")"
					else:
						valStr = apply_format(format,val)
				except ZeroDivisionError:
					valStr = val = divByZero
				newText += [(dest,valStr)]
				safeDict[name] = val

		# reconstruct the line

		numToText = {}

		for (col,text) in newText:
			if (col == None): col = numColumns
			if (col not in numToText): numToText[col] = []
			numToText[col] += [text]

		for (col,text) in enumerate(columns):
			if (col in removeColumns): continue
			if (col not in numToText): numToText[col] = []
			numToText[col] += [text]

		columns = []
		for col in xrange(numColumns+1):
			if (col not in numToText): continue
			columns += numToText[col]

		if (isHeader) and (headerPrefix != None):
			columns[0] = headerPrefix + columns[0]
		print separator.join(columns)


# parse_formula--
#	Parse a string, of the form [<name>:][{<format>}][<dest>=]<expression>
#
# returns (name,format,dest,expression), absent fields are None
# if the parse fails, we return None

def parse_formula(s,destIsColumn=True):
	s = s.strip()

	if (":" not in s):
		name = None
	else:
		(name,s) = s.split(":",1)
		name = name.strip()
		s    = s.strip()

	if ("{" not in s):
		if ("}" in s): return None
		format = None
	else:
		if (not s.startswith("{")) or ("}" not in s): return None
		(format,s) = s[1:].split("}",1)
		format = format.strip()
		s      = s.strip()

	dest = None
	if ("=" in s):
		assignIx      = s.find("=")
		equalIx       = s.find("==")
		unequalIx     = s.find("!=")
		lessOrEqualIx = s.find("<=")
		moreOrEqualIx = s.find(">=")
		if    (assignIx != equalIx) \
		  and (assignIx != unequalIx+1) \
		  and (assignIx != lessOrEqualIx+1) \
		  and (assignIx != moreOrEqualIx+1):
			(dest,s) = s.split("=",1)
			dest = dest.strip()
			s    = s.strip()

	expression = s

	if (destIsColumn) and (dest != None):
		try: dest = parse_column(dest)
		except ValueError: return None

	return (name,format,dest,expression)


# columns_needed--
#	Figure out what input columns are in an expression

def columns_needed(expression):
	columns = {}
	for varName in variable_names(expression):
		try:               columns[parse_column(varName)] = True
		except ValueError: pass
	return [col for col in columns]


def variable_names(expression):
	varFirst = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_"
	varAfter = varFirst + "0123456789"

	first = True
	for ch in expression:
		if (first):
			if (ch in varFirst): (first,varName) = (False,[ch])
		elif (ch in varAfter):
			varName += [ch]
		else:
			yield "".join(varName)
			first = True

	if (not first): yield "".join(varName)


# apply_format--

def apply_format(format,val):
	if (format == None):  return str(val)

	precision = None
	if ("." in format):
		(format,precision) = format.split(".",1)
		format    = format.strip()
		precision = int(precision)

	if (format == "f"):
		if (precision == None): fmt = "%f"
		else:                   fmt = "%%.%df" % precision
		return fmt % val

	if (format == "%"):
		if (precision == None): fmt = "%f%%"
		else:                   fmt = "%%.%df%%%%" % precision
		return fmt % (100*val)

	if (format == "x"):
		if (precision == None): fmt = "0x%X"
		else:                   fmt = "0x%%0%dX" % precision
		return fmt % val

	assert (False), "unsupported format: \"%s\"" % format


# to_value--
#	Parse a value string, allowing an integer, fraction, floating point, or
#	complex value.  Failing any of those, just return the string

def to_value(s):
	if ("/" in s):
		(numer,denom) = s.split("/",1)
		try: return float(numer)/float(denom)
		except ValueError: pass
		except ZeroDivisionError: pass
		try: return int_with_unit(numer)/float(int_with_unit(denom))
		except ValueError: pass
		except ZeroDivisionError: pass
		try: return complex(numer)/complex(denom)
		except ValueError: pass
		except ZeroDivisionError: pass
	else:
		try: return int(s)
		except ValueError: pass
		try: return float(s)
		except ValueError: pass
		try: return int_with_unit(s)
		except ValueError: pass
		if (s in ["J","j"]): return s
		try: return complex(s)
		except ValueError: pass
	return s


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


# parse_column--

def parse_column(s):
	if (not s.startswith("c")): raise ValueError
	col = int(s[1:]) - 1
	if (col < 0): raise ValueError
	return col


# (see reference [1], lybniz2.sourceforge.net/safeeval.html)

def round_int(x): return int(round(x))

def floor_int(x): return int(floor(x))

def ceil_int(x):  return int(ceil(x))

def log2(x):      return log(x) / log(2.0)

def my_random(u=None,v=None):
	# my_random()             --> real value in 0..1
	# my_random(int u)        --> integer value in 1..u
	# my_random(int u, int v) --> integer value in u..v
	# my_random(str u)        --> char value in u
	# my_random(str u, str v) --> str value u or v
	# my_random(list u)       --> choice from u
	# my_random(tuple u)      --> choice from u
	if (u == None):
		return unit_random()
	if (v == None):
		if   (type(u) == str):   return choice(u)
		elif (type(u) == list):  return choice(u)
		elif (type(u) == tuple): return choice(u)
		else:                    return randint(1,u)
	if (type(u) == str) and (type(v) == str): return choice([u,v])
	else:                                     return randint(u,v)

def ord_ascii(x):
	if (len(x) == 1): return ord(x)
	else:             return int(x,16)

def inverse_ord_ascii(x,offset=0):
	if (type(x) != int): return None
	else:                return chr(x+offset % 0x100)

def upper_ascii(x,segment=None):
	if (segment == None): return x.upper()
	else:                 return x[:segment].upper() + x[segment:]

def lower_ascii(x,segment=None):
	if (segment == None): return x.lower()
	else:                 return x[:segment].lower() + x[segment:]

def md5_hash(x,modulus=0x100000000,iterations=1):
	assert (md5_new != None), "wasn't able to import md5"
	for i in xrange(iterations):
		hashVal = md5_new()
		hashVal.update(str(x))
		hashVal = int(hashVal.hexdigest()[:25],16)
		if (modulus != None):
			hashVal = hashVal % modulus
		x = hashVal
	return hashVal


safeList = ["acos", "asin", "atan",  "atan2", "cos", "e",     "exp",
            "fabs", "fmod", "frexp", "ldexp", "log", "log10", "modf",
            "pi",   "pow",  "sin",   "sqrt",  "tan"]
safeDict = dict([(k,locals().get(k,None)) for k in safeList])
safeDict["abs"]     = abs
safeDict["int"]     = int
safeDict["float"]   = float
safeDict["complex"] = complex
safeDict["ceil"]    = ceil_int
safeDict["floor"]   = floor_int
safeDict["round"]   = round_int
safeDict["log2"]    = log2
safeDict["max"]     = max
safeDict["min"]     = min
safeDict["random"]  = my_random
safeDict["ord"]     = ord_ascii
safeDict["invord"]  = inverse_ord_ascii
safeDict["upper"]   = upper_ascii
safeDict["lower"]   = lower_ascii
safeDict["md5"]     = md5_hash


if __name__ == "__main__": main()
