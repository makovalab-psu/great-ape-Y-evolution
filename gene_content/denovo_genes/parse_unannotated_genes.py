
#source activate /galaxy/home/biomonika/conda/GApy2.7
#python

#read balstp out put and parse sequences which are not annotated in the results
#Sample input
#g84.t1  K1920_HUMAN     95.498  658     422     0.0     N/A     N/A     N/A     N/A     K1920_HUMAN Putative uncharacterized protein KIAA1920
#g85.t1  DYN1_BOVIN      75.309  174     81      6.81e-35        N/A     N/A     N/A     N/A     DYN1_BOVIN Dynamin-1
#g87.t1  RBY1C_HUMAN     95.833  154     72      1.85e-37        N/A     N/A     N/A     N/A     RBY1C_HUMAN RNA-binding motif protein, Y chromosome, family 1 member C

from Bio import SeqIO

blastp="blastp_completeGene.outfmt6"
parsed_seq="possible_denovo.fa"

blast_hit_gene={}
with open(blastp,"r") as r:
	for line in r :
		col=line.rstrip('\n').split("\t")
		if col[0].split(".")[0] in blast_hit_gene:
			blast_hit_gene[col[0].split(".")[0]]+=1
		else:
			blast_hit_gene[col[0].split(".")[0]]=0
		





#fasta_nt='panPan.msY.makovalab.ver1.smsk.codingseq_complete'
fasta_nt='ponAbe.msY.makovalab.ver3.smsk.codingseq_complete'
refFH_nt=open(fasta_nt,"rU")
nseq_record={}
for chr_record in SeqIO.parse(refFH_nt, "fasta"):
	print chr_record.id.split(".")[2] #[1] for bonobo; [2] for orangutan
	nseq_record[str(chr_record.id.split(".")[2])]=chr_record #[1] for bonobo; [2] for orangutan

nucl_seq = []
for key in nseq_record.keys():
	if key in blast_hit_gene:
		continue
	else:
		nucl_seq.append(nseq_record[key])

print("Found %i complete genes" % len(nucl_seq))

SeqIO.write(nucl_seq, parsed_seq, "fasta")
