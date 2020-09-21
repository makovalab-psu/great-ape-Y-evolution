#!/usr/bin/python

"""\
Read the FASTA file, For each gene print FASTA file ready for GENECONV

Usage: parse_FASTA_byGENE_individual_blocks.py 
#########################################
Example input:
# >hg.chrX|NLGN4Y;chrY:14796175-14796182;5934358|7|1|156040895;
# AACCCCT
# >hg_Y.chrY|NLGN4Y;chrY:14796175-14796182;42431233|7|-1|57227415;
# AACCCCT
# >panTro_Y.chrY|NLGN4Y;chrY:14796175-14796182;22994382|7|1|26350515;
# AACCCCT
# >panTro.chrX|NLGN4Y;chrY:14796175-14796182;5744150|7|1|151576176;
# AACCCCT
# >panPan_Y.Contig695|NLGN4Y;chrY:14796175-14796182;115478|7|-1|203731;
# AACCCCT
# >panPan.chrX|NLGN4Y;chrY:14796175-14796182;5797621|7|1|155392059;
# AACCCCT
# >gorGor_Y.QB.A_1|NLGN4Y;chrY:14796175-14796182;165826|7|-1|362376;
# AACCCCT
# >gorGor.chrX|NLGN4Y;chrY:14796175-14796182;5585293|7|1|156331669;
# AACCCCT
# >ponAbe_Y.QB.A_125|NLGN4Y;chrY:14796175-14796182;437|7|-1|11920;
# AACTCTT
# >ponAbe.chrX|NLGN4Y;chrY:14796175-14796182;2298411|7|1|151242693;
# AACCCCT


Header FASTA
>hg_Y.chrY|NLGN4Y;chrY:14793201-14793542;42433873|341|-1|57227415;
>chromosome|gene;id                     ;start|size|strand|chrSize

For each gene filter blocks greater than 50bp.
For each block check if there is one copy of X and Y.
Print individual blocks into 

"""\

import sys
from Bio import SeqIO

fasta='XDG_gene_alignments.fa'
refFH=open(fasta,"rU")


gene_object={}
##Gene->ID->XY pair present than assign FASTA
for seq in SeqIO.parse(refFH, "fasta"):
	temp=seq.id.split(";")
	chr=temp[0].split("|")[0]
	if 'hg_Y' in chr:
		chr_type=("Human","Y")
	elif 'hg.chrX' in chr:
		chr_type=("Human","X")
	elif 'panTro_Y' in chr:
		chr_type=("Chimp","Y")
	elif 'panTro.chrX' in chr:
		chr_type=("Chimp","X")
	elif 'gorGor_Y' in chr:
		chr_type=("Gorilla","Y")
	elif 'gorGor.chrX' in chr:
		chr_type=("Gorilla","X")
	elif 'panPan_Y' in chr:
		chr_type=("Bonobo","Y")
	elif 'panPan.chrX' in chr:
		chr_type=("Bonobo","X")
	elif 'ponAbe_Y' in chr:
		chr_type=("Orang","Y")
	elif 'ponAbe.chrX' in chr:
		chr_type=("Orang","X")
	else:
		print "ERRROR: Species and chromosome name not defined."
		break
	gene=temp[0].split("|")[1]
	id=temp[1]
	metadata=temp[2].split("|") #start|size|strand|chrSize
	if gene in gene_object:
		if id in gene_object[gene]:
			if chr_type[0] in gene_object[gene][id]:
				if chr_type[1] in gene_object[gene][id][chr_type[0]]:
					gene_object[gene][id][chr_type[0]][chr_type[1]].append(seq)
				else:
					gene_object[gene][id][chr_type[0]][chr_type[1]]=[seq]
			else:
				gene_object[gene][id][chr_type[0]]={chr_type[1]:[seq]}
		else:
			gene_object[gene][id]={chr_type[0]:{chr_type[1]:[seq]}}
	else:
		gene_object[gene]={id:{chr_type[0]:{chr_type[1]:[seq]}}}


#gene_object['PRKY']['chrY:7276340-7276346']['Human']



#create a directory
#mkdir gene_msablock_fasta

#break fasta by block
for gene in gene_object:
	sum=0
	for loc in gene_object[gene].keys():
		diff=(int(loc.split(':')[1].split('-')[1])-(int(loc.split(':')[1].split('-')[0])))
		if diff>50:
			sum+=diff
			out=open("gene_msablock_fasta/"+gene+"_"+loc+"_XY.fasta","w")
			for spe in gene_object[gene][loc].keys():
				if len(gene_object[gene][loc][spe]) == 2 : 	
					for chr in gene_object[gene][loc][spe]:
						for seq in gene_object[gene][loc][spe][chr]:
							#seq.id=str(spe)+"_"+str(chr)
							#out.write(">"+seq.id+"\n")
							#out.write(seq.seq+"\n")
                                                        SeqIO.write(seq, out, "fasta")
							#print seq.id
	#print gene, sum	




