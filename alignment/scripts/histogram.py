#!/usr/bin/env python
"""
Make a histogram of the values on stdin.
"""

from sys  import argv,stdin,stderr,exit
from math import floor,ceil,log10


def usage(s=None):
	message = """
usage: cat values | histogram [options] > histogram_file
  --bucket=<size>         set size of histogram buckets
  --center=<value>        set the center of some bucket
  --log                   bucketization is logarithmic
  --accuracy=<digits>     set the number of digits for bucket ranges
  --percentage[=<digits>] show bucket counts as percentage of total
  --label=<value>         label for bucket counts
  --input:counts          input is <count>/<value> (instead of just <value>)
  --noempties             don't output empty bins

  The input file (stdin) is usually just one number per line.  But it the
  --input:counts option is used, the input is comprised of a list of item
  counts, like  this (count and number are separated by a tab character)

    (count) (value)  <-- this line not included
 	127     19
 	416     20
 	336     21"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug,showEmptyBins

	# parse args

	bucketSize    = 1.0
	bucketCenter  = None
	logScale      = False
	accuracy      = 3
	pctgDigits    = None
	label         = None
	inputAsCounts = False
	showEmptyBins = True
	debug         = []

	# pick off options

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--bucket=")):
			bucketSize = float_with_unit(argVal)
		elif (arg.startswith("--center=")):
			bucketCenter = float_with_unit(argVal)
		elif (arg == "--log"):
			logScale = True
		elif (arg.startswith("--accuracy=")):
			accuracy = int(argVal)
		elif (arg == "--percentage"):
			pctgDigits = 1
		elif (arg.startswith("--percentage=")):
			pctgDigits = int(argVal)
		elif (arg.startswith("--label=")):
			label = argVal
		elif (arg == "--input:counts"):
			inputAsCounts = True
		elif (arg == "--noempties"):
			showEmptyBins = False
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (bucketCenter == None): bucketCenter = (bucketSize/2.0)

	bucketLeft = bucketCenter - (bucketSize/2.0)

	if (pctgDigits != None):
		pctgFmt = "%%%d.%df%%%%" % (pctgDigits+2,pctgDigits)

	# collect values into buckets

	histogram = {}
	total     = 0

	lineNum = 0
	for line in stdin:
		lineNum += 1
		line = line.strip()

		if (inputAsCounts):
			fields = line.split(None,1)
			count = int(fields[0])
			val   = float(fields[1])
			if (logScale): val = log10(val)
			bucket = floor((val-bucketLeft) / bucketSize)
			if (bucket not in histogram): histogram[bucket] =  count
			else:                         histogram[bucket] += count
			total += count
		else:
			for v in line.split():
				try:
					val = float(v)
				except ValueError:
					assert (False), "invalid number: %s (line %d)" % (v,lineNum)
				if (logScale): val = log10(val)

				bucket = int((val-bucketLeft) / bucketSize)
				if (bucket not in histogram): histogram[bucket] =  1
				else:                         histogram[bucket] += 1
				total += 1

	# print the buckets

	if (label != None):
		print "\t%s" % label

	for (b,count) in all_buckets(histogram):
		if (pctgDigits != None):
			count = pctgFmt % ((100.0*count)/total)
		bucketVal = bucketCenter + b*bucketSize
		if (logScale): bucketVal = 10**bucketVal
		print "%.*f\t%s" % (accuracy,bucketVal,count)


def all_buckets(histogram):
	buckets = [b for b in histogram]
	buckets.sort()
	prevB = None
	for b in buckets:
		if (prevB != None) and (showEmptyBins):
			for bb in xrange(int(prevB)+1,int(b)):
				yield (bb,0)
		yield (b,histogram[b])
		prevB = b


# float_with_unit--
#	Parse a string as an integer, allowing unit suffixes

def float_with_unit(s):
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

	return float(s) * multiplier


if __name__ == "__main__": main()
