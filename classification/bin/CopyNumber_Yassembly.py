#!/usr/bin/python
import sys
import re
import pandas
import numpy
import subprocess
import copy
import argparse
import os.path
import pysam

import warnings
warnings.filterwarnings("error")

parser=argparse.ArgumentParser(
    description='''Ampliconic Gene Copy Number Estimator (AmpliCoNE) : Estimates the copy number of the 9 ampliconic gene families on Human Y chromosome. ''',
    epilog="""Email: v.rahul.simham@gmail.com for errors and bugs""")
parser.add_argument('--BAM', '-b', required=True, help=' Indexed BAM file. (Aligned using BWA MEM) ', metavar='<BAM>')
parser.add_argument('--CHR', '-c', required=True, help=' Chromosome annotation as found in the BAM header file. ', choices=['Y','chrY'])
parser.add_argument('--LENGTH','-l', nargs='?', type=int, default=57227415, help='Length of the Y chromosome in the reference (hg38).(default: %(default)s) ')
parser.add_argument('--READ','-r', nargs='?', default="PAIRED", help='The reads are paired end or single end, if paired we filter for proper read pairs. (default: %(default)s)', choices=['PAIRED','SINGLE'])

#PATH TO TOOLS and REFERENCES
#GFile has gene family definition
##GENE	START	END	TYPE	COPIES_inREF_99perIDENTITY	MAPPABILITY
##GENE - GENE SYMBOL
##START,END - GENE location based on NCBI RefSeq
##TYPE - FAMILY NAME or "Control" All the genes within a family should have same name and Control genes should always be named Control(CASE SENSITIVE)
##COPIES_inREF_99perIDENTITY - Number of regions on reference with great than 99% identity of the whole gene.
##MAPPABILITY - Mappability score pre set per family
#GFile="/nfs/brubeck.bx.psu.edu/scratch5/rahul/GTEx/TESTIS/REF/MAP_CUSTOM/Gene_Definition.tab" 

parser.add_argument('--GFile', '-g', required=True, help=' Gene family definition ', metavar='<TAB>')

#SFile has Y chromosome specific definition HG38 version. Number of rows is equal to length of the chromosome. Position,chrY_M,chrY_RM,chrY_TRF,chrY_GC,chrY_Map
##Position - Each nucleotide position on the chromosome of interest(Y)
##Mappability - Mappability score of that position based on 100 base pairs downstream
##RepeatMasker - Presence or absence of repeats based on repeatmasker output from UCSC reference database values are 0 and 1.
##TandemRepeatFinder - Presence or absence of tandemrepeats based on TandemRepeatFinder output from UCSC reference database values are 0 and 1.
##GC content - 250bp window with this position as center the percentage of G and C in it.
##MAP_Count - Information of sites that unique to gene family. 
#SFile="/nfs/brubeck.bx.psu.edu/scratch5/rahul/GTEx/TESTIS/REF/MAP_CUSTOM/Summary_hg38_chrY_Map_RM_TRF_GC_MapC.tab"	
parser.add_argument('--SFile', '-s', required=True, help=' Y chromosome specific definition HG38 version ', metavar='<TAB>') #Y chromosome specific definition HG38 version. Number of rows is equal to length of the chromosome. Position,chrY_M,chrY_RM,chrY_TRF,chrY_GC,chrY_Map


args=parser.parse_args()

BAM=args.BAM #input MUST be a BWA aligned sorted indexed BAM file
CHR=args.CHR #Y chromosome as defined in the BAM header Y or chrY
len_chr=args.LENGTH #Length of Y chromosome in reference
read_type=args.READ # if single the does not check for proper read pairs use SINGLE
GFile=args.GFile
SFile=args.SFile

#CHECK IF FILES EXIST
if os.path.exists(BAM):
	print "Using "+str(BAM)+" as input BAM file\n"
else:
	print "ERROR: Cannot find input BAM file."
	sys.exit(0)

if os.path.exists(GFile):
	print "Specified path to Gene definition file : "+str(GFile)+"\n"
else:
	print"ERROR: Cannot find Gene definition file "+str(GFile)+" . Please check the code to update the correct path/file name."
	sys.exit(0)

if os.path.exists(SFile):
	print "Specified path to Chromosome Summary file : "+str(SFile)+"\n"
else:
	print"ERROR: Cannot find the file with Chromosome summary file: "+str(SFile)+" with mappability, GCcontent and repeats information. Please check the code to update the correct path/file name."
	sys.exit(0)



def Get_Read_Length(bam_file):
	"""Reads first 1000 reads and obtain the common read length of the sample """
	readlength=[]
	with pysam.AlignmentFile(bam_file, "rb") as bamfile:
		size=1000
		count=0
		for read in bamfile.fetch():
			readlength+=[read.query_length]
			count+=1
			if count>=100:
				break
	return max(set(readlength), key=readlength.count)
			

if Get_Read_Length(BAM) < 100:
	print 'ERROR: Read length of the Input BAM is less than required length of 100 bases. (Tested first 1000 reads)'
	sys.exit(0)

	
#Open the sam file with proper paired reads and filter the reads by alignment
def Count_Matches_CIGAR(cigar_char,cigar_val):
	"""Function to read the parsed CIGAR characters to calculate the number of matches in the first 90 positions """
	i=0 #iterator for while loop
	position=0 #positions in the alignment upto which the cigar string was read
	M=0 #number of matches
	last_val=[] #list to store the position and number of matches before updating in current iteration
	while position < 90 and i <len(cigar_char):
		last_val=[position,M]
		position+=int(cigar_val[i])
		if cigar_char[i]=="M":
			M+=int(cigar_val[i])
		if position>=90: #when the position in alignment crosses the 90 point we want to trim to remove post90 alignment
			extra=position-90 #number of extra positions after 90 
			diff=int(cigar_val[i])-extra #subtract the observed cigar with the extra positions 
			#position=last_val[0]+diff #update the position so we have alignment for first 90 positions
			if cigar_char[i]=="M":
				M=last_val[1]+diff
		i+=1
	return [M,position]

def Count_MisMatch_Deletions_MDZtag(mismatch_char,mismatch_val):
	"""Function to read the parsed MDZtag characters to calculate number of perfect matches in the first 90 positions """
	i=0 #iterator for while loop
	position=0 #positions in the alignment upto which the MDZ tag was read
	MM=0 #number of mismatches and deletions
	while position < 90 and i <len(mismatch_char):
		last_val=[MM,position] #the number of mismatch&deletion before updating in this iteration
		if len(mismatch_char[i])>0:
			if "^" in mismatch_char[i]:
				len_nonM=int(len(mismatch_char[i])-1) #deletion leads with ^, we subtract 1 to ignore the ^
			else:
				len_nonM=int(len(mismatch_char[i]))
		else:
			len_nonM=0
		position+=int(mismatch_val[i])+len_nonM 
		MM+=len_nonM
		if position>=90:
			position=last_val[1]+int(mismatch_val[i])
			if position < 90:
				position=position+len_nonM
				extra=position-90
				diff=len_nonM-extra
				MM=last_val[0]+diff
			else:
				MM=last_val[0]
		i+=1
	return [MM,position]

print "\rFiltering reads for perfect matches"	
#Read the input bam file and parse for proper read pairs and then look for reads with atleast 88 perfect matches in the first 90 base pairs of the read
filtered_read_start_position=BAM+"_"+CHR+"_alignmentSTARTPosition.tab"
with pysam.AlignmentFile(BAM, "rb") as bamfile, open(filtered_read_start_position, "w") as w:
	j,i=0,0
	for read in bamfile.fetch(CHR):
		if read_type=="PAIRED":
			if read.flag == 99 or read.flag == 163 or read.flag == 83 or read.flag == 147:
				cigar_char=re.split('\d+',read.cigarstring)[1:] #parse the character
				cigar_val=re.split('\D+',read.cigarstring)[:-1] #parse the number
				mismatch=read.get_tag('MD')
				mismatch_char=re.split('\d+',mismatch)[1:] #parse the character
				mismatch_val=re.split('\D+',mismatch)	   #parse the number
				if int(cigar_val[0])>=90 and int(mismatch_val[0])>=90:
					w.write(str(read.query_name)+"\t"+str(read.reference_start+1)+"\n")
				elif int(cigar_val[0])>=90 and int(mismatch_val[0])< 90:
					MM=Count_MisMatch_Deletions_MDZtag(mismatch_char,mismatch_val)
					if MM[0] <=2 and MM[1]>=90:
						w.write(str(read.query_name)+"\t"+str(read.reference_start+1)+"\n")				
				elif int(cigar_val[0])<90:
					M=Count_Matches_CIGAR(cigar_char,cigar_val)
					if M[0]>=88 and M[1]>=90:
						if int(mismatch_val[0])>=88:
							w.write(str(read.query_name)+"\t"+str(read.reference_start+1)+"\n")
						else:
							MM=Count_MisMatch_Deletions_MDZtag(mismatch_char,mismatch_val)
							if MM[0] <=2 and MM[1]>=90:
								w.write(str(read.query_name)+"\t"+str(read.reference_start+1)+"\n")		
				i+=1
				if i==1000000:
					i=0
					j+=1
					print "\rProcessed "+str(j)+"000000 lines"
		elif read_type=="SINGLE":
			cigar_char=re.split('\d+',read.cigarstring)[1:] #parse the character
			cigar_val=re.split('\D+',read.cigarstring)[:-1] #parse the number
			mismatch=read.get_tag('MD')
			mismatch_char=re.split('\d+',mismatch)[1:] #parse the character
			mismatch_val=re.split('\D+',mismatch)	   #parse the number
			if int(cigar_val[0])>=90 and int(mismatch_val[0])>=90:
				w.write(str(read.query_name)+"\t"+str(read.reference_start+1)+"\n")
			elif int(cigar_val[0])>=90 and int(mismatch_val[0])< 90:
				MM=Count_MisMatch_Deletions_MDZtag(mismatch_char,mismatch_val)
				if MM[0] <=2 and MM[1]>=90:
					w.write(str(read.query_name)+"\t"+str(read.reference_start+1)+"\n")				
			elif int(cigar_val[0])<90:
				M=Count_Matches_CIGAR(cigar_char,cigar_val)
				if M[0]>=88 and M[1]>=90:
					if int(mismatch_val[0])>=88:
						w.write(str(read.query_name)+"\t"+str(read.reference_start+1)+"\n")
					else:
						MM=Count_MisMatch_Deletions_MDZtag(mismatch_char,mismatch_val)
						if MM[0] <=2 and MM[1]>=90:
							w.write(str(read.query_name)+"\t"+str(read.reference_start+1)+"\n")		
			i+=1
			if i==1000000:
				i=0
				j+=1
				print "\rProcessed "+str(j)+"000000 lines"

print "\rFinished filtering reads"
print "\rObtaining chromosome wide read start counts."

##Get the read start counts from alignment start positions
counts={}
with open(filtered_read_start_position,"r") as r:
	for line in r :
		col=line.rstrip('\n').split("\t")
		if col[1] in counts:
			counts[col[1]]+=1
		else:
			counts[col[1]]=1
	


##For every position on Y chromosome, if no reads starting then fill with zeros
##Output file is a one-based position specific counts

temp=0
val=[]
for pos in sorted(counts.keys(),key=int):
	diff=int(pos)-temp
	if diff==1:
		val+=[counts[pos]]
	else:
		val+=[0]*(diff-1)
		val+=[counts[pos]]
	temp=int(pos)

#####len_chr=57227415	#Y=57227415
tail_cov=len_chr-temp
val+=[0]*(tail_cov)


output=BAM+"_"+CHR+"_ReadStartCount.txt"
with open(output, "w") as file:  
	file.write("StartCount\n")
	for index in range(len(val)):
		file.write(str(val[index])+ "\n")

# val_db=list() #temporary reading instead of creating the file for the debugging
# with open("sim_to_38.bam_Y_ReadStartCount.txt","r") as rsc:
# 	rsc.next()
# 	for line in rsc:
# 		col=line.rstrip('\n').split("\t") 
# 		val_db.append(col[0])

	
##Get the list of genes whose read depth is needed to be calculated		
print "\rObtaining the gene list"

Gene_list={}
Family_list={}
with open(GFile, "rU") as g: 
  #Expand the cigar string 
  header=g.next()	#GENE	START	END	TYPE	COPIES_inREF_99perIDENTITY	MAPPABILITY
  for line in g:
	col=line.rstrip('\n').split("\t") 
	Gene_list[col[0]]=col
	if col[3] in Family_list:
		Family_list[col[3]]+=1
	else:
		Family_list[col[3]]=0


##Parse repeat region and add mappability values to RSCcounts
print "\rLoading the Read start counts (RSC) and Mappability values"
Data={} # Dictionary to temporarily store values and convert to data-frame later . Saves runtime
with open(SFile, "rU") as s: 
	header=s.next() #'Position\tMappability\tRepeatMasked\tTRF\tGC\tMapCount\n'
	for row in s:
		col=row.rstrip('\n').split("\t") 
		if int(col[2])==1 and int(col[3])==1 :
			if float(col[4])>0 and float(col[1])>0:
				#Data[int(col[0])]={'Position': int(col[0]), 'Mappability':float(col[1]), 'GC':float(col[4]), 'RSC':int(val[int(col[0])-1]), 'MapCount':int(col[5]) } #Converting 1 based to 0 based by subtracting 1
				Data[int(col[0])]=[int(col[0]), float(col[1]), float(col[4]), int(val[int(col[0])-1]),int(col[5])]


Summary_data= pandas.DataFrame.from_dict(Data, orient="index")				
Summary_data.columns=['Position', 'Mappability', 'GC', 'RSC','MapCount']
Control_data=Summary_data.copy()
Control_data=Control_data.loc[Control_data['Mappability']==1]


Summary_data=Summary_data.values
Control_data=Control_data.values
whole_M=numpy.mean(Control_data[:,3])

print "\rPerforming the GC correction"
###CONTROL & GC correction
GCmean=numpy.empty((0, 1))
for i in range(1,100):
	counts_temp=Control_data[((Control_data[:,2]>=i)==(Control_data[:,2]<(i+1))).nonzero()][:,3]
	if len(counts_temp) > 0 :
		gcm=numpy.mean(counts_temp)
	else:
		gcm=0
	#print i, gcm
	GCmean = numpy.append(GCmean,gcm)

#This step below is to make sure there are no 0 to divide with in the next step.
GCmean[GCmean==0]=whole_M
Correction=(whole_M/GCmean)
for i in range(len(Correction)):
	if numpy.isnan(Correction[i]):
		Correction[i]=1
	if numpy.isinf(Correction[i]):
		Correction[i]=1

GCcor_Summary_data=Summary_data.copy()
for i in range(1,100):
	id=((GCcor_Summary_data[:,2]>=(i))==(GCcor_Summary_data[:,2]<(i+1))).nonzero()
	GCcor_Summary_data[id,3]=GCcor_Summary_data[id][:,3]*Correction[i-1]
	
####IMPORTANT ndarray[r,c]=assign the values, ndarray[r][:,c] does not work

##################################################PSEUDOCODE
# #Location of gene on Y chromosome
# Gene_location=(start,end)
# #Mappability value based on number of gene copies on the reference
# Expected_mappability= MAP 
# #Matrix with mappability and coverage values for each position on Y chr 
# Ychr_summary=[Position,Mappability,Coverage]
#
# Control region coverage(Ychr_summary) 
# {
	# unique_sites= empty list
	# #all the positions that are unique to Y chromosome
	# For each position in the Ychr
# {
	# If ( Ychr_summary[position,2] == 1 )  
# {
		# unique_sites.append(Ychr_summary[position,])
		# }
	# }
	# Return average(unique_sites)
# }
##############################################
def Control_region_coverage(Y_Summary_data) : # Ychr_summary=['Position', 'Mappability', 'GC', 'RSC','MapCount']
		Control_region=Y_Summary_data[(Y_Summary_data[:,1]==1).nonzero()]
		try:
			mean_of_control_region=numpy.mean(Control_region[:,3])
		except RuntimeWarning:
			print("Oops! The mean of control region could not be calculated.")
		return mean_of_control_region
	

##################################################PSEUDOCODE
# Gene Coverage(Gene_location,Expected_mappability,Ychr_summary) 
# {
# Family_specific_sites= empty list
# For each position in the Ychr
# {	
# #if position overlaps with gene location
# If (position>= Gene_location[1] and position < Gene_location[2]) 
# {
	# #mappabillity of the position is same as expected mappability
	# If  (Ychr_summary[position,2] == Expected_mappability) 
# {
		# #obtain the positions that are unique to the gene family
 # Family_specific_sites.append(Ychr_summary[position,])
	# }
# }
# Return average(Family_specific_sites[,3])
# }
#################################################
def Get_window_coverage(Gene_info,Y_Summary_data) : # Y_Summary_data=['Position', 'Mappability', 'GC', 'RSC','MapCount']
	Gene_Summary_data=Y_Summary_data[((Y_Summary_data[:,0]>int(Gene_info[1]))==(Y_Summary_data[:,0]<int(Gene_info[2]))).nonzero()] #parse region of the gene (Gene_info[1] is start,Gene_info[2] is end)
	return [Gene_info[0],Gene_info[3],numpy.mean(Gene_Summary_data[:,3])]  #return(gene name,family_type(type),mean coverage)


print "\rObtaining the gene level RSC"
Control_coverage=Control_region_coverage(GCcor_Summary_data)
Temp_Coverage={}
for gene in Gene_list:
	#print gene
	tmp_window_coverage=Get_window_coverage(Gene_list[gene],GCcor_Summary_data)
	if (numpy.isnan(tmp_window_coverage[2])):
		print(gene + " is being excluded due to the unassigned value for the GC content.")
	else:
		Temp_Coverage[gene]=tmp_window_coverage #only assign if the gene/window has the value, nan is excluded

Gene_coverage= pandas.DataFrame.from_dict(Temp_Coverage, orient="index")				
Gene_coverage.columns=['Gene', 'Family', 'RSCDepth']

XDG_Genes=Gene_coverage.values[(Gene_coverage.values[:,1]=="Control").nonzero()]

try:
	XDG_Control_coverage=numpy.mean(XDG_Genes[:,2])
except ZeroDivisionError:
	print("Oops! XDG_Control_coverage could not be calculated.")
	print(XDG_Genes[:,2])


XDG_Control_coverage=numpy.mean(XDG_Genes[:,2])

XDG_CopyNumber=XDG_Genes.copy()
XDG_CopyNumber[:,2]=XDG_CopyNumber[:,2]/Control_coverage

Allgene=Gene_coverage.values
Allgene_CopyNumber=Allgene.copy()
Allgene_CopyNumber[:,2]=Allgene_CopyNumber[:,2]/Control_coverage
Allgene_CopyNumber[Allgene_CopyNumber[:,1].argsort()]

print "\rCaliculating CN values and printing"
Ampliconic_Genes=Gene_coverage.values[(Gene_coverage.values[:,1]!="Control").nonzero()]
AG_CopyNumber=Ampliconic_Genes.copy()
AG_CopyNumber[:,2]=AG_CopyNumber[:,2]/Control_coverage
AG_CopyNumber_XDGbased=Ampliconic_Genes.copy()
AG_CopyNumber_XDGbased[:,2]=AG_CopyNumber_XDGbased[:,2]/XDG_Control_coverage

AG_out=BAM+"window_CN_summary.txt"
with open(AG_out, "w") as AGfile:  
	AGfile.write("GeneFamily\tCopyNumber(MAP=1)\tCopyNumber(XDG)\n")
	for family in sorted(Family_list.iterkeys()):
		if family == "Control": 
			continue
		Family_set=AG_CopyNumber[(AG_CopyNumber[:,1]==family).nonzero()]
		Family_setXDG=AG_CopyNumber_XDGbased[(AG_CopyNumber_XDGbased[:,1]==family).nonzero()]
		#print family+"\t"+str(numpy.sum(Family_set[:,2]))+"\t"+str(numpy.sum(Family_setXDG[:,2]))
		AGfile.write(family+"\t"+str(numpy.sum(Family_set[:,2]))+"\t"+str(numpy.sum(Family_setXDG[:,2]))+ "\n")


XDG_out=BAM+"XDG_CopyNumber.txt"
with open(XDG_out, "w") as XDGfile:  
	XDGfile.write("Gene\tCopyNumber\n")		
	for i in XDG_CopyNumber:
		#print i[0]+"\t"+str(i[2])
		XDGfile.write(i[0]+"\t"+str(i[2])+"\n")

ALL_out=BAM+"Window_CopyNumber.txt"
with open(ALL_out, "w") as ALLfile:  
	ALLfile.write("Gene\tCopyNumber\n")		
	for i in Allgene_CopyNumber[Allgene_CopyNumber[:,1].argsort()]:
		#print i[0]+"\t"+str(i[2])
		ALLfile.write(i[0]+"\t"+str(i[2])+"\n")
