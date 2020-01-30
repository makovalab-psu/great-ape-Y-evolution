#!/usr/bin/python

"""\
Read the MAF file, first line is human specific. Obtain the blocks which fall within the palindrome.
For each block parse the size of the block and sum it for each species to get overall alignment count.
If there are multiple sequences from same species within block, longest alignment is considered.

Usage: palindrome_coverage_human.py <P1-8>
#########################################
Example input:
#a score=0
#s hg_Y.chrY         15879049 152 + 57227415 atctgtagtcccagctactcgagaggctgaggcaggagaatctcttgaacctgggag
#gctgaggttgcagtgagctgagatcctgccactgcactcctgcctgggtgacagagcgagactccgtctcaaatatatatatatatatattatat
#s hg_Y.chrY         41326666 152 - 57227415 atctgtagtcccagctactcgagaggctgaggcaggagaatctcttgaacctgggag
#gctgaggttgcagtgagctgagatcctgccactgcactcctgcctgggtgacagagcgagactccgtctcaaatatatatatatatatattatat
#s ponAbe_Y.QB.A_3     888820 149 -  1144157 acctgtagtcccagctactcgagaggctgaggcaggagaatctcttgaacctggaag
#gctgaggttgcagtgagctgagatcatgccactgcactcctgcctgggtgacagagtgagactccatctcaaaaatatatatatatata---tat
#s gorGor_Y.QB.A_150    24897 135 -    33903 acctgtagtcccagctactcgaggggctgaggcaggagaatctcttgaacctgggag
#gttgaggttgcaatgagccgagatcctgccactgcactcctgcctgggtgacagagcgagactccatctcaaaTATAT-----------------
#s panTro_Y.chrY     23522068 135 + 26350515 acctgtagtcccagctactcgagaggctgaggcaggagaatctcttgaacctgggag
#gttgaggttgcagtgagccgagatcctgccactgcactcctgcctgggtgacagagcgagactccatctcaaatatat-----------------
#s panTro_Y.chrY      2807824 135 - 26350515 acctgtagtcccagctactcgagaggctgaggcaggagaatctcttgaacctgggag
#gttgaggttgcagtgagccgagatcctgccactgcactcctgcctgggtgacagagcgagactccatctcaaatatat-----------------

##Steps: For each palindrome if the first row overlaps with palindrom
##			Then obtain the size of the block for each species (example 152 is for human)
##			Obtain the over all sum of block sizes per each species	
"""\

import sys
from Bio import AlignIO

if len(sys.argv[1:]) != 1:
  sys.exit(__doc__)

#Set the palindrome here and run each time
inputP=sys.argv[1] #"P1"

#Palindrome arm coordinates on hg38
Palindrome_coordinates={"P1":[('23359067', '24822577')],
"P2":[('23061889', '23208197')],
"P3":[('21924954', '22208730')],
"P4":[('18450291', '18640356')],
"P5":[('17455877', '17951255')],
"P6":[('16159590', '16269541')],
"P7":[('15874906', '15883575')],
"P8":[('13984498', '14019652')]
}

#List of species names in the MAF file
species_list=["panPan_Y", "panTro_Y", "gorGor_Y", "ponAbe_Y", "hg_Y"]

#all Y chromosome alignment file.
#/nfs/brubeck.bx.psu.edu/scratch6/rahul/Bonobo_Y/analysis/palindrome/palindrome_coverage/multi_alignment_based
input_handle = open("msa/alignment.hg_Y_centric.20191126.maf", "rU") 
output_handle = open(inputP+"_sizeNonRepeatBlock.tab", "w")

#Obtaining the palindrome of interest
Palindrome=Palindrome_coordinates[inputP]

#File handle of MAF
alignments = AlignIO.parse(input_handle, "maf")

#Block overlap cutoff

percent_cutoff=0.95
pal_size=(int(Palindrome[0][1])-int(Palindrome[0][0]))+1
sequence_Bon=[0]*pal_size
sequence_Gor=[0]*pal_size
sequence_Orang=[0]*pal_size
sequence_Chimp=[0]*pal_size
#sequence_Human=[0]*pal_size
#Reading each block in the MAF file
for msa in alignments :
    for region in Palindrome: #Obtain the location of palindromes
        region_start=region[0]
        region_end=region[1]
        #temp dict to store the block size of each species within the block
        for row in msa:
            if row.id== 'hg_Y.chrY':
                guide_row=row 
                #Convert - strand coordinates to positive strand coordinates 
                if int(guide_row.annotations['strand'])==1:
                    msa_start=int(guide_row.annotations['start'])
                    msa_end=int(guide_row.annotations['start'])+(int(guide_row.annotations['size'])-1) # 1 based
                    human_seq=str(guide_row.seq)
                    #If the size of human alignment is above 0
                    if int(guide_row.annotations['size'])>0:
                        #Parse the palindrome specific blocks
                        if ((int(msa_start)) >= int(region_start)) and ((int(msa_end)) <= int(region_end)): 
                            for record in msa:
                                #for each row in the block save the maximum block size value.
                                #print (float(record.annotations['size'])/guide_row.annotations['size'])
                                if(record.id).split(".",1)[0]== "panPan_Y":
                                    if (float(record.annotations['size'])/guide_row.annotations['size']) > percent_cutoff :
                                        start=msa_start-int(Palindrome[0][0])
                                        seq_t=str(record.seq)
                                        for i in range(len(seq_t)):
                                            if human_seq[i].isupper():
                                                if seq_t[i] != 'N' and seq_t[i] != '-' :
                                                    sequence_Bon[i+start]=1
                                elif(record.id).split(".",1)[0]== "gorGor_Y":
                                    if (float(record.annotations['size'])/guide_row.annotations['size']) > percent_cutoff :
                                        start=msa_start-int(Palindrome[0][0])
                                        seq_t=str(record.seq)
                                        for i in range(len(seq_t)):
                                            if human_seq[i].isupper():
                                                if seq_t[i] != 'N' and seq_t[i] != '-' :
                                                    sequence_Gor[i+start]=1
                                elif(record.id).split(".",1)[0]== "ponAbe_Y":
                                    if (float(record.annotations['size'])/guide_row.annotations['size']) > percent_cutoff :
                                        start=msa_start-int(Palindrome[0][0])
                                        seq_t=str(record.seq)
                                        for i in range(len(seq_t)):
                                            if human_seq[i].isupper():
                                                if seq_t[i] != 'N' and seq_t[i] != '-' :
                                                    sequence_Orang[i+start]=1	
                                elif(record.id).split(".",1)[0]== "panTro_Y":
                                    if (float(record.annotations['size'])/guide_row.annotations['size']) > percent_cutoff :
                                        start=msa_start-int(Palindrome[0][0])
                                        seq_t=str(record.seq)
                                        for i in range(len(seq_t)):
                                            if human_seq[i].isupper():
                                                if seq_t[i] != 'N' and seq_t[i] != '-' :
                                                    sequence_Chimp[i+start]=1														
                else:
                    msa_start=int(guide_row.annotations['srcSize'])-((int(guide_row.annotations['start'])-1)+(int(guide_row.annotations['size'])-1))
                    msa_end=int(guide_row.annotations['srcSize'])-(int(guide_row.annotations['start'])-1)
                    human_seq=str(guide_row.seq)[::-1]
                    #If the size of human alignment is above 0
                    if int(guide_row.annotations['size'])>0:
                        #Parse the palindrome specific blocks
                        if ((int(msa_start)) >= int(region_start)) and ((int(msa_end)) <= int(region_end)): 
                            for record in msa:
                                #for each row in the block save the maximum block size value.
                                #print (float(record.annotations['size'])/guide_row.annotations['size'])
                                if(record.id).split(".",1)[0]== "panPan_Y":
                                    if (float(record.annotations['size'])/guide_row.annotations['size']) > percent_cutoff :
                                        start=msa_start-int(Palindrome[0][0])
                                        seq_t=str(record.seq)[::-1]
                                        for i in range(len(seq_t)):
                                            if human_seq[i].isupper():
                                                if seq_t[i] != 'N' and seq_t[i] != '-' :
                                                    sequence_Bon[i+start]=1
                                elif(record.id).split(".",1)[0]== "gorGor_Y":
                                    if (float(record.annotations['size'])/guide_row.annotations['size']) > percent_cutoff :
                                        start=msa_start-int(Palindrome[0][0])
                                        seq_t=str(record.seq)[::-1]
                                        for i in range(len(seq_t)):
                                            if human_seq[i].isupper():
                                                if seq_t[i] != 'N' and seq_t[i] != '-' :
                                                    sequence_Gor[i+start]=1
                                elif(record.id).split(".",1)[0]== "ponAbe_Y":
                                    if (float(record.annotations['size'])/guide_row.annotations['size']) > percent_cutoff :
                                        start=msa_start-int(Palindrome[0][0])
                                        seq_t=str(record.seq)[::-1]
                                        for i in range(len(seq_t)):
                                            if human_seq[i].isupper():
                                                if seq_t[i] != 'N' and seq_t[i] != '-' :
                                                    sequence_Orang[i+start]=1							
                                elif(record.id).split(".",1)[0]== "panTro_Y":
                                    if (float(record.annotations['size'])/guide_row.annotations['size']) > percent_cutoff :
                                        start=msa_start-int(Palindrome[0][0])
                                        seq_t=str(record.seq)[::-1]
                                        for i in range(len(seq_t)):
                                            if human_seq[i].isupper():
                                                if seq_t[i] != 'N' and seq_t[i] != '-' :
                                                    sequence_Chimp[i+start]=1	


output=inputP+"_Orangutan_alignment_vector.txt"
with open(output, "w") as file:  
	file.write("Orangutan\n")
	for index in range(len(sequence_Orang)):
		file.write(str(sequence_Orang[index])+ "\n")

output=inputP+"_Bonobo_alignment_vector.txt"
with open(output, "w") as file:  
	file.write("Bonobo\n")
	for index in range(len(sequence_Bon)):
		file.write(str(sequence_Bon[index])+ "\n")

output=inputP+"_Gorilla_alignment_vector.txt"
with open(output, "w") as file:  
	file.write("Gorilla\n")
	for index in range(len(sequence_Gor)):
		file.write(str(sequence_Gor[index])+ "\n")

output=inputP+"_Chimpanzee_alignment_vector.txt"
with open(output, "w") as file:  
	file.write("Chimpanzee\n")
	for index in range(len(sequence_Chimp)):
		file.write(str(sequence_Chimp[index])+ "\n")