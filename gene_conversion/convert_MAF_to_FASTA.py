#!/usr/bin/python

"""\
Read the MAF file, if line is human specific. Obtain the blocks which fall within the gene as FASTA.

Usage: convert_MAF_to_FASTA.py <SRY>
#########################################
Example input:
# a
# s hg.chrX             6505847   14    + 156040895 ccctctctggaaca
# s hg_Y.chrY           15921882  14    + 57227415  ccctctctggaaca
# s panTro_Y.chrY       23563878  14    + 26350515  ccctctctggaaca
# s panTro.chrX         6307695   14    + 151576176 tcctctctggaaca
# s panPan_Y.Contig3719 106804    14    - 113384    ccctctctggaaca
# s panPan.chrX         6360357   14    + 155392059 tcctctctggaaca
# s gorGor_Y.QB.A_150   30506     14    + 33903     ccctctctggaaca
# s gorGor_Y.QB.A_52    83736     14    - 89995     ccctctctggaaca
# s gorGor.chrX         6167006   14    + 156331669 ccctctctggaaca
# s ponAbe_Y.QB.A_154   7367      14    + 8956      ccctctctggaaca
# s ponAbe.chrX         2871750   14    + 151242693 ccctctctggaaca


Header FASTA
>hg_Y.chrY|NLGN4Y;chrY:14793201-14793542;42433873|341|-1|57227415;
>chromosome|gene;id                     ;start|size|strand|chrSize

"""\

import sys
from Bio import AlignIO


#XDG coordinates on hg38
# SRY	2786855	2787699
# AMELY	6865918	6874027
# DBY	12904831	12920478
# EIF1AY	20575725	20593154
# NLGN4Y	14523746	14843726
# PRKY	7273972	7381547
# KDM5D	19705417	19744939
# TBL1Y	6910686	7091683
# TMSB4Y	13703567	13706024
# USP9Y	12701231	12860839
# UTY	13248379	13480673
# ZFY	2935477	2982506
Gene_coordinates={
"SRY":[('2786855', '2787699')],
"AMELY":[('6865918', '6874027')],
"DBY":[('12904831', '12920478')],
"EIF1AY":[('20575725', '20593154')],
"NLGN4Y":[('14523746', '14843726')],
"PRKY":[('7273972', '7381547')],
"KDM5D":[('19705417', '19744939')],
"TBL1Y":[('6910686', '7091683')],
"TMSB4Y":[('13703567', '13706024')],
"USP9Y":[('12701231', '12860839')],
"UTY":[('13248379', '13480673')],
"ZFY":[('2935477', '2982506')]
}

#List of species names in the MAF file
species_list=["panPan_Y", "panTro_Y", "gorGor_Y", "ponAbe_Y", "hg_Y"]

#all X and Y chromosome alignment file.
#/nfs/brubeck.bx.psu.edu/scratch6/rahul/GeneConversion/analysis/cactus_align/
input_handle = open("XY_human/alignment.hg_centric.chrXY.guided.20191201.maf", "rU") 
output_handle = open("XDG_gene_alignments.fa", "w")


#File handle of MAF
alignments = AlignIO.parse(input_handle, "maf")


#Output dictionary to which the block sizes are added each time
Gene={}
overall_size={"panPan_Y":0, "panTro_Y":0, "gorGor_Y":0, "ponAbe_Y":0, "hg_Y":0}	
overallnonRepeat_size={"panPan_Y":0, "panTro_Y":0, "gorGor_Y":0, "ponAbe_Y":0, "hg_Y":0}	
#Reading each block in the MAF file
for msa in alignments :
	for gene in Gene_coordinates.keys(): #Obtain the location of gene
		region_start=Gene_coordinates[gene][0][0]
		region_end=Gene_coordinates[gene][0][1]
		#temp dict to store the block size of each species within the block
		temp_MAF={}
		species_blockSize={}
		species_notRepeatSize={}
		for row in msa:
			if row.id== 'hg_Y.chrY':
				guide_row=row 
				#Convert - strand coordinates to positive strand coordinates 
				if int(guide_row.annotations['strand'])==1:
					msa_start=int(guide_row.annotations['start'])
					msa_end=int(guide_row.annotations['start'])+int(guide_row.annotations['size'])
				else:
					msa_start=int(guide_row.annotations['srcSize'])-int(guide_row.annotations['start'])-int(guide_row.annotations['size'])
					msa_end=int(guide_row.annotations['srcSize'])-int(guide_row.annotations['start'])
				#If the size of human alignment is above 0
				if int(guide_row.annotations['size'])>0:
				#Parse the gene specific blocks
					if ((int(msa_start)) >= int(region_start)) and ((int(msa_end)) <= int(region_end)): 
						for record in msa:
							#for each row in the block save the maximum block size value.
							output_handle.write(">"+str(record.id)+"|"+str(gene)+";"+str("chrY")+":"+str(msa_start)+"-"+str(msa_end)+";"+str(record.annotations['start'])+"|"+str(record.annotations['size'])+"|"+str(record.annotations['strand'])+"|"+str(record.annotations['srcSize'])+";"+"\n")
							output_handle.write(str(record.seq)+"\n")


output_handle.close()
input_handle.close()