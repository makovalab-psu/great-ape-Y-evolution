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

#Palindrome arm coordinates on chimpanzee
# Palindrome	Start	End	HumanHomologs	Amplicon	Approx. arm length	Left_arm_end (Approx)			
# C1	1759451	2,053,069		Pink	146809	1906260			
# C2	2298081	2,984,818		Blue	343368.5	2641450			
# C3	3587737	3,925,944	P1 and P2	Red	169103.5	3756841			
# C4	4669973	5,444,881	P5	Gold	387454	5057427	#First Gold is from 699000 on dotplot	699687	1082500
# C5	8651546	9,099,775		Turquoise	224114.5	8875661			
# C6	9099776	9,383,863		Pink	142043.5	9241820			
# C7	9383864	9,832,093		Turquoise	224114.5	9607979			
# C8	9832094	10,116,181		Pink	142043.5	9974138			
# C9	10116182	10,564,411		Turquoise	224114.5	10340297			
# C10	10564412	10,851,248		Pink	143418	10707830			
# C11	11060051	11,674,261		Blue	307105	11367156			
# C12	12651797	12,963,129	P1 and P2	Red	155666	12807463			
# C13	13340000	14,084,660	P5	Gold	372330	13712330			
# C14	14798116	15,024,316		Pink	113100	14911216			
# C15	15493746	16,056,815		Blue	281534.5	15775281			
# C16	16477310	16,791,733		Pink	157211.5	16634522			
# C17	21591500	21,671,300	P8		39900	21631400			
# C18	23517939	23,546,819	P7		14440	23532379			
# C19	23807052	24,577,396	P6		385172	24192224			

Palindrome_coordinates={"C1":[('1759451', '1906260')],
"C2":[('2298081', '2641450')],
"C3":[('3587737', '3756841')],
"C4":[('4669973', '5057427')],
"C5":[('8651546', '8875661')],
"C6":[('9099776', '9241820')],
"C7":[('9383864', '9607979')],
"C8":[('9832094', '9974138')],
"C9":[('10116182', '10340297')],
"C10":[('10564412', '10707830')],
"C11":[('11060051', '11367156')],
"C12":[('12651797', '12807463')],
"C13":[('13340000', '13712330')],
"C14":[('14798116', '14911216')],
"C15":[('15493746', '15775281')],
"C16":[('16477310', '16634522')],
"C17":[('21591500', '21631400')],
"C18":[('23517939', '23532379')],
"C19":[('23807052', '24192224')]
}

#List of species names in the MAF file
species_list=["panPan_Y", "panTro_Y", "gorGor_Y","ponAbe_Y", "hg_Y"]

#all Y chromosome alignment file.
#/nfs/brubeck.bx.psu.edu/scratch6/rahul/Bonobo_Y/analysis/palindrome/palindrome_coverage/multi_alignment_based
input_handle = open("msa/alignment.panTro_Y_centric.20191126.maf", "rU") 

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
#sequence_Chimp=[0]*pal_size
sequence_Human=[0]*pal_size
#Reading each block in the MAF file
for msa in alignments :
    for region in Palindrome: #Obtain the location of palindromes
        region_start=region[0]
        region_end=region[1]
        #temp dict to store the block size of each species within the block
        for row in msa:
            if row.id== 'panTro_Y.chrY':
                guide_row=row 
                #Convert - strand coordinates to positive strand coordinates 
                if int(guide_row.annotations['strand'])==1:
                    msa_start=int(guide_row.annotations['start'])
                    msa_end=int(guide_row.annotations['start'])+(int(guide_row.annotations['size'])-1) # 1 based
                    chimp_seq=str(guide_row.seq)
                    #If the size of human alignment is above 0
                    if int(guide_row.annotations['size'])>0:
                        #Parse the palindrome specific blocks
                        if ((int(msa_start)) >= int(region_start)) and ((int(msa_end)) <= int(region_end)): 
                            for record in msa:
                                #for each row in the block save the maximum block size value.
                                if(record.id).split(".",1)[0]== "panPan_Y":
                                    if (float(record.annotations['size'])/guide_row.annotations['size']) > percent_cutoff :
                                        start=msa_start-int(Palindrome[0][0])
                                        seq_t=str(record.seq)
                                        for i in range(len(seq_t)):
                                            if chimp_seq[i].isupper():  ###Checking if human sequence has a repeat and not the species sequences
                                                if seq_t[i] != 'N' and seq_t[i] != '-' :
                                                    sequence_Bon[i+start]=1
                                elif(record.id).split(".",1)[0]== "gorGor_Y":
                                    if (float(record.annotations['size'])/guide_row.annotations['size']) > percent_cutoff :
                                        start=msa_start-int(Palindrome[0][0])
                                        seq_t=str(record.seq)
                                        for i in range(len(seq_t)):
                                            if chimp_seq[i].isupper():
                                                if seq_t[i] != 'N' and seq_t[i] != '-' :
                                                    sequence_Gor[i+start]=1
                                elif(record.id).split(".",1)[0]== "ponAbe_Y":
                                    if (float(record.annotations['size'])/guide_row.annotations['size']) > percent_cutoff :
                                        start=msa_start-int(Palindrome[0][0])
                                        seq_t=str(record.seq)
                                        for i in range(len(seq_t)):
                                            if chimp_seq[i].isupper():
                                                if seq_t[i] != 'N' and seq_t[i] != '-' :
                                                    sequence_Orang[i+start]=1
                                elif(record.id).split(".",1)[0]== "hg_Y":
                                    if (float(record.annotations['size'])/guide_row.annotations['size']) > percent_cutoff :
                                        start=msa_start-int(Palindrome[0][0])
                                        seq_t=str(record.seq)
                                        for i in range(len(seq_t)):
                                            if chimp_seq[i].isupper():
                                                if seq_t[i] != 'N' and seq_t[i] != '-' :
                                                    sequence_Human[i+start]=1
                else:
                    msa_start=int(guide_row.annotations['srcSize'])-((int(guide_row.annotations['start'])-1)+(int(guide_row.annotations['size'])-1))
                    msa_end=int(guide_row.annotations['srcSize'])-(int(guide_row.annotations['start'])-1)
                    chimp_seq=str(guide_row.seq)[::-1]
                    #If the size of human alignment is above 0
                    if int(guide_row.annotations['size'])>0:
                        #Parse the palindrome specific blocks
                        if ((int(msa_start)) >= int(region_start)) and ((int(msa_end)) <= int(region_end)): 
                            for record in msa:
                                #for each row in the block save the maximum block size value.
                                if(record.id).split(".",1)[0]== "panPan_Y":
                                    if (float(record.annotations['size'])/guide_row.annotations['size']) > percent_cutoff :
                                        start=msa_start-int(Palindrome[0][0])
                                        seq_t=str(record.seq)[::-1]
                                        for i in range(len(seq_t)):
                                            if chimp_seq[i].isupper():  ###Checking if human sequence has a repeat and not the species sequences
                                                if seq_t[i] != 'N' and seq_t[i] != '-' :
                                                    sequence_Bon[i+start]=1
                                elif(record.id).split(".",1)[0]== "gorGor_Y":
                                    if (float(record.annotations['size'])/guide_row.annotations['size']) > percent_cutoff :
                                        start=msa_start-int(Palindrome[0][0])
                                        seq_t=str(record.seq)[::-1]
                                        for i in range(len(seq_t)):
                                            if chimp_seq[i].isupper():
                                                if seq_t[i] != 'N' and seq_t[i] != '-' :
                                                    sequence_Gor[i+start]=1
                                elif(record.id).split(".",1)[0]== "ponAbe_Y":
                                    if (float(record.annotations['size'])/guide_row.annotations['size']) > percent_cutoff :
                                        start=msa_start-int(Palindrome[0][0])
                                        seq_t=str(record.seq)[::-1]
                                        for i in range(len(seq_t)):
                                            if chimp_seq[i].isupper():
                                                if seq_t[i] != 'N' and seq_t[i] != '-' :
                                                    sequence_Orang[i+start]=1
                                elif(record.id).split(".",1)[0]== "hg_Y":
                                    if (float(record.annotations['size'])/guide_row.annotations['size']) > percent_cutoff :
                                        start=msa_start-int(Palindrome[0][0])
                                        seq_t=str(record.seq)[::-1]
                                        for i in range(len(seq_t)):
                                            if chimp_seq[i].isupper():
                                                if seq_t[i] != 'N' and seq_t[i] != '-' :
                                                    sequence_Human[i+start]=1
                        


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
		
output=inputP+"_Human_alignment_vector.txt"
with open(output, "w") as file:  
	file.write("Human\n")
	for index in range(len(sequence_Human)):
		file.write(str(sequence_Human[index])+ "\n")