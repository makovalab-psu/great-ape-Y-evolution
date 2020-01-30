
#source activate /galaxy/home/biomonika/conda/GApy2.7
#python

#Input file provided by Monika.
#QB.A_154/0_8927 hg38.chrY 0 1 -314 98.4615 40401582 40401647 57227415 6509 6574 8927 1334
#QB.A_494/0_2540 hg38.chrY 0 1 -280 100 49364982 49365038 57227415 1897 1953 2540 1142
#QB.A_653/0_1970 hg38.chrY 0 0 -9041 96.1013 23111531 23113503 57227415 0 1965 1970 41393
#QB.A_784/0_1681 hg38.chrY 0 0 -220 100 1456944 1456988 57227415 1634 1678 1681 890
from Bio import SeqIO

contigs="list_YmappedContigs/Yhits.bonobo.allreads.jelly.out.fasta.reference.mapping.txt"
#contigs="list_YmappedContigs/Yhits.sorang_metassembler.fasta.reference.mapping.txt"

blastp="blastp_completeGene_above80identifty.outfmt6"

blast_hit_gene={}
with open(blastp,"r") as r:
	for line in r :
		col=line.rstrip('\n').split("\t")
		if float(col[2])>80:
			if col[0].split(".")[0] in blast_hit_gene:
				blast_hit_gene[col[0].split(".")[0]].append(col[1])
			else:
				blast_hit_gene[col[0].split(".")[0]]=[col[1]]


human_chimp_hit_contigs={}
with open(contigs,"r") as r:
	for line in r :
		col=line.rstrip('\n').split(" ")
		if col[0].split("/")[0] in human_chimp_hit_contigs:
			human_chimp_hit_contigs[col[0].split("/")[0]]+=1
		else:
			human_chimp_hit_contigs[col[0].split("/")[0]]=0
		

#RUN THIS PART FOR BONOBO
#Bonobo		
fasta_nt='panPan.msY.makovalab.ver1.smsk.codingseq_complete'

refFH_nt=open(fasta_nt,"rU")
nseq_record={}
for chr_record in SeqIO.parse(refFH_nt, "fasta"):
	temp=chr_record.id.split(".")
	#id=temp[0]+"."+temp[1] #for orangutan
	id=temp[0] #for bonobo
	if id in human_chimp_hit_contigs:
		#print chr_record.id.split(".")[1] #[2] for orangutan
		if chr_record.id.split(".")[1] in blast_hit_gene:
			print chr_record.id.split(".")[1], blast_hit_gene[chr_record.id.split(".")[1]] 
		


print "Looking at denovo genes in list"
denovo='possible_denovo.fa'
denFH_nt=open(denovo,"rU")
for chr_record in SeqIO.parse(denFH_nt, "fasta"):
	temp=chr_record.id.split(".")
	#id=temp[0]+"."+temp[1] #for orangutan
	id=temp[0] #for bonobo
	if id in human_chimp_hit_contigs:
		#print chr_record.id.split(".")[1] #[2] for orangutan
		if chr_record.id.split(".")[1] in blast_hit_gene:
			print chr_record.id.split(".")[1], blast_hit_gene[chr_record.id.split(".")[1]] 
		

		
		
		
#### RUN THIS PARTT FFOR ORANGUTAN	
#Orangutan
fasta_nt='ponAbe.msY.makovalab.ver3.smsk.codingseq_complete'
refFH_nt=open(fasta_nt,"rU")
nseq_record={}
for chr_record in SeqIO.parse(refFH_nt, "fasta"):
	temp=chr_record.id.split(".")
	id=temp[0]+"."+temp[1] #for orangutan
	#id=temp[0] #for bonobo
	if id in human_chimp_hit_contigs:
		print chr_record.id.split(".")[2] #[2] for orangutan
		if chr_record.id.split(".")[2] in blast_hit_gene:
			print chr_record.id.split(".")[2], blast_hit_gene[chr_record.id.split(".")[2]] 
			
print "Looking at denovo genes in list"
denovo='possible_denovo.fa'
denFH_nt=open(denovo,"rU")
for chr_record in SeqIO.parse(denFH_nt, "fasta"):
	temp=chr_record.id.split(".")
	#id=temp[0]+"."+temp[1] #for orangutan
	id=temp[0] #for bonobo
	if id in human_chimp_hit_contigs:
		#print chr_record.id.split(".")[1] #[2] for orangutan
		if chr_record.id.split(".")[1] in blast_hit_gene:
			print chr_record.id.split(".")[1], blast_hit_gene[chr_record.id.split(".")[1]] 
		




