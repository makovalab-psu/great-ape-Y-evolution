#!/bin/bash
#SBATCH --job-name=evaluate
#SBATCH --output=evaluate-%j.out
#SBATCH --error=evaluate-%j.err
#SBATCH --mem=150G
#SBATCH --ntasks=32


#load conda environment
source activate /galaxy/home/biomonika/conda/3Dgenome > /dev/null

#set -e
set -x

while getopts :i:a:r:t:s:x:y: option
do
case "${option}"
in
i) index=${OPTARG};; #index="pacbio_assembly_HQ_Y_scaffolds.fa"
a) assembly_name=${OPTARG};; #assembly_name="pacbio_HQ_Y_GIAB" #assembly_name="hg38"
r) reads=${OPTARG};; #reads="HG003_1"
t) lowerCNthreshold=${OPTARG};; #copy number for X-degenerate regions (<)
s) higherCNthreshold=${OPTARG};; #copy number for ampliconic regions (>)
x) ampliconic_fasta=${OPTARG};; #fasta file with species-specific ampliconic sequences
y) XDEG_fasta=${OPTARG};; #fasta file with species-specific XDEG sequences
esac
done

bam_file=${reads}_on_${assembly_name}

#RUN THE CLASSIFICATION INTO AMPLICONIC VS OTHER REGIONS
if [ ! -f ${assembly_name}/"other.txt" ]; then
directory=$(pwd) #current working directory
Rscript --vanilla bin/classify.R ${assembly_name}/${index}".MAPQ.txt" ${bam_file}".bamwindow_CN_summary.txt" ${directory}"/"${assembly_name} ${lowerCNthreshold} ${higherCNthreshold}
mv Rplots.pdf ${assembly_name}.ColoredThresholds.pdf
echo "Classification into ampliconic and other regions finished."

#sort the files
sort ${assembly_name}/other.txt -o ${assembly_name}/other.txt
sort ${assembly_name}/ampliconic.txt -o ${assembly_name}/ampliconic.txt

seqtk subseq ${index} ${assembly_name}/ampliconic.txt >${assembly_name}/ampliconic.fasta
seqtk subseq ${index} ${assembly_name}/other.txt >${assembly_name}/other.fasta

RepeatMasker -species primates -pa 16 ${assembly_name}/ampliconic.fasta &
RepeatMasker -species primates -pa 16 ${assembly_name}/other.fasta &
echo "Repeatmasking of classified scaffolds finished."
fi

#EVALUATE IF THE CLASSIFICATION WORKED USING WHOLE SCAFFOLDS
if [ ! -f ${assembly_name}/Xdeg_in_species.psl ]; then
	blat ${index} ${ampliconic_fasta} ${assembly_name}/amplicon_in_species.psl #first 5 lines are a header and need to be skipped
	blat ${index} ${XDEG_fasta} ${assembly_name}/Xdeg_in_species.psl #first 5 lines are a header and need to be skipped

	#require at least half of the query to be present and the identity over 95%
	#psl explanation at https://useast.ensembl.org/info/website/upload/psl.html
	for a in ${assembly_name}/*.psl; do cat $a | tail -n +6 | sort -rgk1 | awk '{if ($1>=30) {print $0}}' | awk '{if (($1/$11)>=0.2) {print $0}}' | awk '{if (($1/($1+$2))>=0.99) print;}' | cut -f14 | sort | uniq >${a}.txt; done;

	echo "FALSE;   Ampliconic classified as other:"
	comm -12 ${assembly_name}/amplicon_in_species.psl.txt ${assembly_name}/other.txt | wc -l

	echo "CORRECT; Ampliconic classified as ampliconic:"
	comm -12 ${assembly_name}/amplicon_in_species.psl.txt ${assembly_name}/ampliconic.txt | wc -l

	echo "CORRECT; Xdeg classified as other:"
	comm -12 ${assembly_name}/Xdeg_in_species.psl.txt ${assembly_name}/other.txt | wc -l

	echo "FALSE;   Xdeg classified as ampliconic:"
	comm -12 ${assembly_name}/Xdeg_in_species.psl.txt ${assembly_name}/ampliconic.txt | wc -l

	echo "Scaffold evaluation finished."

fi

wait
EVALUATE IF THE CLASSIFICATION WORKED USING WINDOWS
if [ ! -f ${assembly_name}/Xdeg_in_species.psl.gff ]; then
	echo "Evaluation using windows:"
	for a in ${assembly_name}/*.psl; do cat $a | sort -rgk1 | awk '{if ($1>=30) {print $0}}' | awk '{if (($1/$11)>=0.2) {print $0}}' | awk '{if (($1/($1+$2))>=0.99) print;}' | awk -v OFS='\t' '{print $14,$10,"gene_hit",$16<$17?$16:$17,$17<$16?$16:$17,($1/$11),".",".",$14}' | bedtools sort >${a}.gff; done #create gff file for genes
	for a in ${assembly_name}/*.psl.gff; do echo $a; bedtools intersect -a ${assembly_name}/classification.gff -b $a | bedtools groupby -i stdin -g 9 -c 6 -o mean; done #intersect classification results for windows with gene coordinates
fi
echo "Evaluation finished. Done."
echo "#######################"