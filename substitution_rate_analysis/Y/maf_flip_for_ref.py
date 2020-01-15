#!/usr/bin/env python
"""
Application to convert MAF so that the reference species is always + stranded
-----------------------------------------------------------------------------

:Author: Bob Harris (rsharris@bx.psu.edu)
:Version: $Revision: $

The application reads a MAF file from standard input and writes an equivalent
file to standard out.  Any blocks for which the reference species is on the
negative strand are flipped.  Any other blocks, including those not containing
the reference species, are left unchanged.
"""

import sys
import bx.align.maf

def usage(s=None):
	message = """
maf_flip_for_ref ref_species < maf_file > maf_file

If the input maf has sequence names of the form <species>.<chromosome>,
ref_species should be <species> *and* the dot that follows it.

If the input maf has sequence names of the form <species> and nothing else,
then ref_species should be <species> (without a dot).

If the input maf has sequence names that are just chromosomes, scaffolds, or
contigs, ref_species should be the chromosome, scaffold, or contig.
"""
	if (s == None): sys.exit (message)
	else:           sys.exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	if (len(sys.argv) != 2):
		usage("wrong number of command line arguments")

	refSpecies = sys.argv[1]

	# read the alignments and other info
	#
	# nota bene: get_component_by_src_start(refSpecies) returns the first
	# component in the block that has a prefix matching refSpecies

	out = bx.align.maf.Writer(sys.stdout)

	for mafBlock in bx.align.maf.Reader(sys.stdin):

		c = mafBlock.get_component_by_src_start(refSpecies)
		if (c != None):
			if (c.strand == "-"): flip_block (mafBlock)

		out.write (mafBlock)


def flip_block(mafBlock):
	for ix in range(0,len(mafBlock.components)):
		c = mafBlock.components[ix]
		mafBlock.components[ix] = c.reverse_complement()


if __name__ == "__main__": main()

