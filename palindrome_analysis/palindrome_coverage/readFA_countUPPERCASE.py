#!/usr/bin/env python

"""\
Usage: python readFA_countUC.py <FASTA>


"""\

from Bio import SeqIO
import sys

if len(sys.argv[1:]) != 1:
  sys.exit(__doc__)

#FASTA 
fasta=sys.argv[1]
refFH=open(fasta,"rU")

for record in SeqIO.parse(refFH, "fasta"):
	print sum(1 for nt in str(record.seq) if nt.isupper()) 
