
#source activate /galaxy/home/biomonika/conda/GApy2.7
#python


# Example input file format.
# Contig170       AUGUSTUS        gene    5789    6403    1       -       .       g2
# Contig170       AUGUSTUS        transcript      5789    6403    .       -       .       g2.t1
# Contig170       AUGUSTUS        stop_codon      5789    5791    .       -       0       transcript_id "g2.t1"; gene_id "g2";
# Contig170       AUGUSTUS        CDS     5789    6403    .       -       0       transcript_id "g2.t1"; gene_id "g2";
# Contig170       AUGUSTUS        start_codon     6401    6403    .       -       0       transcript_id "g2.t1"; gene_id "g2";
from Bio import SeqIO

gene_complete={}
#with open("panPan.msY.makovalab.ver1.smsk.gff","r") as r:
with open("ponAbe.msY.makovalab.ver3.smsk.gff","r") as r:
	for line in r :
		col=line.rstrip('\n')
		if col[0] == "#":
			continue
		else:
			temp=col.split("\t")
			if temp[2]=="gene":
				#print temp
				if temp[8] in gene_complete:
					print str(temp[8])+"  Found twice in file. Code cannot hande"
					break
				else:
					gene_complete[temp[8]]={}
					gene_complete[temp[8]]['chr']=temp[0]
					gene_complete[temp[8]]['start']=temp[3]
					gene_complete[temp[8]]['end']=temp[4]
					#gene_complete[temp[2]]['start_codon']=""
					#gene_complete[temp[2]]['stop_codon']=""
			elif temp[2]=="stop_codon":
					gene_complete[temp[8].split(";")[1].split(" ")[2].replace("\"","")]['stop_codon']=temp[3]
			elif temp[2]=="start_codon":
					gene_complete[temp[8].split(";")[1].split(" ")[2].replace("\"","")]['start_codon']=temp[3]
			#elif temp[2]=="CDS":
					#gene_complete[temp[8].split(";")[1].split(" ")[2].replace("\"","")]['CDS'].append(str(temp[3]+"-"+temp[4]))








# fasta_aa='panPan.msY.makovalab.ver1.smsk.aa'
# output_aa="panPan.msY.makovalab.ver1.smsk.aa_complete"

fasta_aa='ponAbe.msY.makovalab.ver3.smsk.aa'
output_aa="ponAbe.msY.makovalab.ver3.smsk.aa_complete"


refFH_aa=open(fasta_aa,"rU")
aseq_record={}
for chr_record in SeqIO.parse(refFH_aa, "fasta"):
	#print chr_record.id.split(".")[0]
	aseq_record[str(chr_record.id.split(".")[0])]=chr_record

coding_seq = []
for key in gene_complete.keys():
	if "start_codon" in gene_complete[key] and "stop_codon" in gene_complete[key]:
		coding_seq.append(aseq_record[key])

print("Found %i complete genes" % len(coding_seq))

SeqIO.write(coding_seq, output_aa, "fasta")



# fasta_nt='panPan.msY.makovalab.ver1.smsk.codingseq'
# output_nt="panPan.msY.makovalab.ver1.smsk.codingseq_complete"

fasta_nt='ponAbe.msY.makovalab.ver3.smsk.codingseq'
output_nt="ponAbe.msY.makovalab.ver3.smsk.codingseq_complete"

refFH_nt=open(fasta_nt,"rU")
nseq_record={}
for chr_record in SeqIO.parse(refFH_nt, "fasta"):
	print chr_record.id.split(".")[2]   #[1] for bonobo; [2] for orangutan
	nseq_record[str(chr_record.id.split(".")[2])]=chr_record  #[1] for bonobo; [2] for orangutan

nucl_seq = []
for key in gene_complete.keys():
	if "start_codon" in gene_complete[key] and "stop_codon" in gene_complete[key]:
		nucl_seq.append(nseq_record[key])

print("Found %i complete genes" % len(nucl_seq))

SeqIO.write(nucl_seq, output_nt, "fasta")