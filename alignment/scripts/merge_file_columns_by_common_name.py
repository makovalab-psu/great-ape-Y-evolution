#!/usr/bin/env python
"""
Merge columns from multiple files
=================================

The key column can be set on a per-file basis, like this:

	merge_file_columns_by_common_name file1[7] file2[5]

Thus in file1 we will look for the keys in column 7, and in file2 we will look
for the keys in column 5.

"""

from sys import argv,stdin,stdout,stderr,exit
from re  import compile


filenameRe  = compile("^(?P<name>.+)\[(?P<key>[0-9]+)\]$")
filename2Re = compile("^(?P<name>.+)\[(?P<lowKey>[0-9]+):(?P<highKey>[0-9]+)\]$")


def usage(s=None):
	message = """
usage: merge_file_columns_by_common_name [options]
  <filename>[<column>]  (cumulative) filename to read some columns from; if
                        <column> is given it is used as the key column; if
                        <column> is not given --key is used
  --key=<column>        column to use as key for files for which key is not
                        specified
  --ignoreduplicates    ignore all but the first occurence of keys that appear
                        more than once in a file
                        (by default these are considered to be errors)
  --processduplicates   (not yet implemented)
  --keepcomments        treat comment lines the same as any other; comments
                        are lines beginning with #
  --removekey           remove all key columns and add the key as the first
                        output column
  --tabdelimited        output is tab-delimited
  --spacedelimited      output is space-delimited
  --delimiter=<string>  field delimiter for output
  --missing=<string>    replace missing columns with this string

Merge columns from multiple files.

The key column can be set on a per-file basis, like this:
	merge_file_columns_by_common_name file1[7] file2[5]

Thus in file1 we will look for the keys in column 7, and in file2 we will look
for the keys in column 5."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse command line

	masterKeyColumn = 0
	handleExtraKeys = "abort"
	ignoreComments  = True
	removeKey       = False
	delimiter       = "\t"
	inputFilenames  = []
	fillMissingWith = None

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--key=")) or (arg.startswith("--column=")):
			masterKeyColumn = argVal
			if (":" in masterKeyColumn):
				masterKeyColumn = masterKeyColumn.split(":",1)
				masterKeyColumn = (int(masterKeyColumn[0])-1,int(masterKeyColumn[1])-1)
				assert (0 <= masterKeyColumn[0] < masterKeyColumn[1])
			else:
				masterKeyColumn = int(masterKeyColumn) - 1
				assert (masterKeyColumn >= 0)
		elif (arg == "--ignoreduplicates"):
			handleExtraKeys = "discard"
		elif (arg == "--processduplicates"):
			handleExtraKeys = "process"
		elif (arg == "--keepcomments"):
			ignoreComments = False
		elif (arg == "--removekey"):
			removeKey = True
		elif (arg == "--tabdelimited"):
			delimiter = "\t"
		elif (arg == "--spacedelimited"):
			delimiter = " "
		elif (arg.startswith("--delimiter=")):
			delimiter = argVal
		elif (arg.startswith("--missing=")):
			fillMissingWith = argVal
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			inputFilenames += [arg]

	if (inputFilenames == []):
		inputFilenames += [""]	# ("" means read from stdin)

	if (fillMissingWith != None):
		fillMissingWith = delimiter.join(fillMissingWith.split())

	# process the files/lines

	nameOrder   = []
	namesSeen   = {}
	nameToText  = None
	totalFields = 0

	for (fIx,filename) in enumerate(inputFilenames):
		if (type(masterKeyColumn) == tuple):
			(keyColumnLo,keyColumnHi) = masterKeyColumn
		else:
			keyColumnLo = keyColumnHi = masterKeyColumn

		if (filename == ""):
			(f,filename) = (stdin,"stdin")
			closeFile = False
		else:
			m = filenameRe.match(filename)
			if (m != None):
				filename  = m.group("name")
				keyColumnLo = keyColumnHi = int(m.group("key")) - 1
				assert (keyColumnLo >= 0)
			else:
				m = filename2Re.match(filename)
				if (m != None):
					filename  = m.group("name")
					keyColumnLo = int(m.group("lowKey" )) - 1
					keyColumnHi = int(m.group("highKey")) - 1
					assert (0 <= keyColumnLo < keyColumnHi)

			f = file(filename,"rt")
			closeFile = True

		nameToNewText = {}
		numFields = None
		lineNumber = 0
		for line in f:
			lineNumber += 1
			line = line.strip()
			if (ignoreComments) and (line.startswith("#")): continue
	
			fields = line.split()
			numFields = len(fields)
			if (numFields <= keyColumnLo): continue
			if (keyColumnLo == keyColumnHi):
				name = fields[keyColumnLo]
			else:
				name = delimiter.join(fields[keyColumnLo:keyColumnHi+1])

			if (handleExtraKeys == "discard"):
				if (name in nameToNewText): continue
			elif (handleExtraKeys == "abort"):
				assert (name not in nameToNewText), \
				       "%s appears more than once in %s" % (name,filename)
			else:
				assert (False), \
				       "sorry, handleExtraKeys == \"%s\" not yet supported" % handleExtraKeys

			if (removeKey):
				fields = fields[:keyColumnLo] + fields[keyColumnHi+1:]
				numFields = len(fields)
				line = delimiter.join(fields)

			if (name not in namesSeen):
				nameOrder += [name]
				namesSeen[name] = True

			nameToNewText[name] = line
			#print "nameToNewText[%s] = %s" % (name,line)

		if (closeFile): f.close()

		if (nameToText == None):
			nameToText = {}
			for name in nameOrder:
				nameToText[name] = [nameToNewText[name]]
			totalFields += numFields
			continue

		for name in nameOrder:
			if (name in nameToNewText):
				line = nameToNewText[name]
			else:
				if (fillMissingWith == None):
					assert (name in nameToNewText), \
					       "%s does not appear in %s" % (name,filename)
				line = delimiter.join([fillMissingWith] * numFields)
			if (name in nameToText): nameToText[name] += [line]

		for name in nameToNewText:
			if (name in nameToText): continue
			if (fillMissingWith == None):
				assert (False), \
					   "%s does not appear in %s" % (name,inputFilenames[0])
			nameToText[name] = ([fillMissingWith] * totalFields) + [nameToNewText[name]]

		totalFields += numFields

	for name in nameOrder:
		if (removeKey): text = [name] + nameToText[name]
		else:           text =          nameToText[name]
		print delimiter.join(text)


if __name__ == "__main__": main()
