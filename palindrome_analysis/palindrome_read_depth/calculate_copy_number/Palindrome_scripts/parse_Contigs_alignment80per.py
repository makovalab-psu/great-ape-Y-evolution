#!/usr/bin/env python
"""\
Identifies windows mapping to palindromes with >80% alignment from windows_for_identity.dat file.
usage: parse_Contigs_alignment80per.py <ContigNamesPalindrome> <windows_for_identity> <output>

Example input files:
<ContigNamesPalindrome>
Contig1000_0_1365
Contig1023_29000_30000
Contig1038_0_1447
Contig1047_0_1000
Contig106_0_1000
Contig106_1000_2359

<windows_for_identity>
#contig start   end     Ns      alignments      m       mm      gapGor  gapHum
Contig4681      1       1123            0       0       0       0       0
Contig5332      1       1277            0       0       0       0       0
Contig2574      1       1000    0       0       0       0       0       0


For each contig the code obtains the lenght of the contig which is end - start (+1) and divides the number of m by lenght of contig.
Identify contigs with >80% alignment.
Check if number of N's is less than 200 assuming the size of the window is always 1000.
"""
import sys

if len(sys.argv[1:]) != 3:
  sys.exit(__doc__)

Contigs=sys.argv[1]
#creating a dictionary of contig names
contigs_dic=dict()
with open(Contigs, "rU") as r:
	for line in r:
		n=line.rstrip('\n')
		if n in contigs_dic:
			contigs_dic[n]+= 1
		else:
			contigs_dic[n]=0

r.close()

alignment=sys.argv[2]
filtered=sys.argv[3]

with open(alignment, "rU") as s,open(filtered, "w") as w:
	s.next()
	for line in s:
		col=line.split("\t")
		id=str(col[0])+"_"+str(int(col[1])-1)+"_"+str(col[2])
		if id in contigs_dic:
			M=float(col[5])
			size=int(col[2])-int(col[1])+1
			Per=M/size*100
			if Per >= 80:
				if col[3] == '':
					continue
				elif int(col[3]) <= 200:
					w.write(line)

w.close()
s.close()