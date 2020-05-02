#!/bin/bash
#SBATCH --job-name=run_classification
#SBATCH --output=run_classification-%j.out
#SBATCH --error=run_classification-%j.err
#SBATCH --mem=500G
#SBATCH --ntasks=64


#load conda environment
source activate /galaxy/home/biomonika/conda/3Dgenome > /dev/null

set -e
#set -x

while getopts :i:a:r:f:g:m:l: option
do
case "${option}"
in
i) index=${OPTARG};; #index="pacbio_assembly_HQ_Y_scaffolds.fa"
a) assembly_name=${OPTARG};; #assembly_name="pacbio_HQ_Y_GIAB" #assembly_name="hg38"
r) reads=${OPTARG};; #reads="HG003_1"
f) female=${OPTARG};; #provide female genome WITHOUT the Y chromosome
g) gene_file=${OPTARG};; #Trusted_genes
m) GEM_FOLDER=${OPTARG};; #"/galaxy/home/rxv923/src/GEM-binaries-Linux-x86_64-core_2-20130406-045632/bin"
l) fmappability=${OPTARG};; #mappability mode; either "female" or "Y"
esac
done

echo "index: ${index} assembly_name: ${assembly_name} reads: ${reads} female: ${female} GEM_FOLDER: ${GEM_FOLDER} "
GEM_FOLDER="/galaxy/home/rxv923/src/GEM-binaries-Linux-x86_64-core_2-20130406-045632/bin"
REPEATMASKER_LIB_DIR=/nfs/brubeck.bx.psu.edu/scratch6/monika/project_Y_chromosome/repeatmasker_library

#check if the software exists
hash bwa 2>/dev/null || { echo "I require bwa in path but it's not installed.  Aborting."; exit 1; }
hash blat 2>/dev/null || { echo "I require blat in path but it's not installed.  Aborting."; exit 1; }
hash seqtk 2>/dev/null || { echo "I require seqtk in path but it's not installed.  Aborting."; exit 1; }
hash samtools 2>/dev/null || { echo "I require samtools in path but it's not installed.  Aborting."; exit 1; }
hash RepeatMasker 2>/dev/null || { echo "I require RepeatMasker in path but it's not installed.  Aborting."; exit 1; }
hash wig2bed 2>/dev/null || { echo "I require wig2bed in path but it's not installed.  Aborting."; exit 1; }
hash bedtools 2>/dev/null || { echo "I require bedtools in path but it's not installed.  Aborting."; exit 1; }

#check if the python scripts are available
if [ ! -f "bin/analyzeRepeats.py" ]; then
    echo "bin/analyzeRepeats.py not found!"
    exit 1;
fi
if [ ! -f "bin/generate_positional_information.py" ]; then
    echo "bin/generate_positional_information.py not found!"
    exit 1;
fi
if [ ! -f "bin/analyzeGCcontent.py" ]; then
    echo "bin/analyzeGCcontent.py not found!"
    exit 1;
fi
if [ ! -f "bin/createSummary.R" ]; then
    echo "bin/createSummary.R not found!"
    exit 1;
fi
if [ ! -f ${female} ]; then
    echo "Female genome not found!"
    exit 1;
fi

################################# PREPARE THE GENE FILES #################################

#create a folder for the file creation
mkdir -p ${assembly_name}

if [ ! -f ${assembly_name}/${index}.fai ]; then
    samtools faidx ${index}
    mv ${index}.fai ${assembly_name}/${index}.fai
    awk -v OFS='\t' {'print $1,$2'} ${assembly_name}/${index}.fai > ${assembly_name}/${index}.txt
    echo -e "NA\t0\t0" >${assembly_name}/${index}.offset.txt
    awk '{sum+=$2; printf $0"\t""%.0f\n", sum}' ${assembly_name}/${index}.txt >>${assembly_name}/${index}.offset.txt
fi

if [ ! -f ${assembly_name}/${assembly_name}.windows.bed ]; then
    #create windows from the reference
    bedtools makewindows -g ${assembly_name}/${index}.txt -w 5000 -s 2000 >${assembly_name}/${assembly_name}.windows.bed 
    #to compare all 5-kb sequence segments, in 2-kb steps, to the entire remainder of the MSY sequence
fi


if [ ! -f ${assembly_name}/${assembly_name}.nsq ]; then
    echo "Blast database has not been created yet."
    makeblastdb -in ${index} -parse_seqids -dbtype nucl -out ${assembly_name}
    mv ${assembly_name}.nhr ${assembly_name}.nin ${assembly_name}.nog ${assembly_name}.nsd ${assembly_name}.nsi ${assembly_name}.nsq ${assembly_name}
fi

if [ ! -f ${assembly_name}/${gene_file}.${assembly_name}.bed ]; then
    echo "Trusted genes have not yet been blasted against the Y reference."
    #the length of at least 100bp, identity at least 95% and query coverage at least 85%
    blastn -query ${gene_file}.fasta -db ${assembly_name}/${assembly_name} -outfmt 6 -perc_identity 95 -qcov_hsp_perc 85 -evalue 0.01 | awk '{if ($4>99) print}' >${assembly_name}/${gene_file}.${assembly_name}.outfmt 
    cat ${assembly_name}/${gene_file}.${assembly_name}.outfmt | awk '{print "chrY" "\t" $9 "\t" $10}' | awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t' | sort -k1,1n -k2,2n -k3,3n >${assembly_name}/${gene_file}.${assembly_name}.bed
fi

if [ ! -f ${assembly_name}/window_definition.${assembly_name}.tab ]; then
    echo "The windows file with controls has not yet been created."
    awk '{print "chrY\t"$2"\t"$3"\t"$1"\t"1"\t1"}' ${assembly_name}/${assembly_name}.windows.bed >${assembly_name}/wfile.bed
    cat ${assembly_name}/wfile.bed | awk '{print $1"\t"$2"\t"$3}' >${assembly_name}/wfile_small.bed 
    #NOW INTERSECT GENES WITH WINDOWS TO DETERMINE CONTROLS
    bedtools intersect -u -a ${assembly_name}/wfile_small.bed -b ${assembly_name}/${gene_file}.${assembly_name}.bed >${assembly_name}/intersection.bed
    python bin/insert_intersection_as_control.py ${assembly_name}/intersection.bed ${assembly_name}/wfile.bed ${assembly_name}
fi

#now switch to the continous definitions
if [ ! -f ${assembly_name}/"window_definition."${assembly_name}".continous.tab" ]; then
    python bin/switch_to_continous.py ${assembly_name}/"window_definition."${assembly_name}".tab" ${assembly_name}/${index}.offset.txt >${assembly_name}/"window_definition."${assembly_name}".continous.tab" 
    #FINAL OUTPUT
fi

wait
#rm ${assembly_name}/wfile_small.bed ${assembly_name}/wfile.bed ${assembly_name}/intersection.bed ${assembly_name}/window_definition.${assembly_name}.tab

echo "Preparation done."



################################# RUN THE CALCULATION #################################

bam_file=${reads}_on_${assembly_name}
window_file="window_definition."${assembly_name}".continous.tab" #"100.wfile.controls.formatted.txt" #"wfile.controls.formatted.txt" #"#window_file="/nfs/brubeck.bx.psu.edu/scratch5/rahul/GTEx/TESTIS/REF/MAP_CUSTOM/Gene_Definition.tab"
reference_length=`awk '{s+=$2} END {printf "%.0f", s}' ${assembly_name}/${index}.txt`



if [ ! -f ${assembly_name}/${assembly_name}.continous.fa ]; then
	#create continous referemce
	echo ">"${assembly_name} >${assembly_name}/${assembly_name}.continous.fa
	cat ${index} | bioawk -c fastx '{ print $seq }' | tr '\n' ' ' | tr -d ' ' >>${assembly_name}/${assembly_name}.continous.fa
	echo "" >>${assembly_name}/${assembly_name}.continous.fa
fi

chromosome_name=`bioawk -c fastx '{ print $name }' ${assembly_name}/${assembly_name}.continous.fa`

if [ ! -f ${bam_file}.bam ]; then
    #map reads

    #we need to concatenate female genome with Y chromosome 
    cat ${female} ${assembly_name}/${assembly_name}.continous.fa >${assembly_name}/female_genome_and_${assembly_name}.continous.fa

    if [ ! -f ${assembly_name}/female_genome_and_${assembly_name}.continous.fa.bwt ]; then
        bwa index ${assembly_name}/female_genome_and_${assembly_name}.continous.fa #only generate bwt index if it doesn't exist yet
    fi

    #map reads to the concatenation of female and Y
    #bwa mem -t 63 ${assembly_name}/female_genome_and_${assembly_name}.continous.fa ${reads}1.fastq ${reads}2.fastq >${bam_file}.sam
    bwa mem -t 63 ${assembly_name}/female_genome_and_${assembly_name}.continous.fa ${reads}1.fastq >${bam_file}.sam


    samtools view -F 4 -bhS ${bam_file}.sam >${bam_file}.bam #exclude unmapped reads and convert to bam
    samtools sort --threads 64 ${bam_file}.bam >${bam_file}.sorted.bam
    mv ${bam_file}.sorted.bam ${bam_file}.bam

    #subset to the Y chromosome/assembly
    samtools index ${bam_file}.bam
    samtools view -hb ${bam_file}.bam ${assembly_name} >${bam_file}.Y.bam
    #awk -v chr="$assembly_name" '{if ($3==chr) print;}' ${bam_file}.sam >${bam_file}.sam.Y
    samtools sort --threads 64 ${bam_file}.Y.bam >${bam_file}.Y.sorted.bam
    mv ${bam_file}.Y.sorted.bam ${bam_file}.bam
    samtools index ${bam_file}.bam
    #rm ${bam_file}.sam ${bam_file}.Y.bam
fi

if [ ! -f ${assembly_name}/female_genome_and_${assembly_name}.continous.fa_MAPPABILITY.bed ]; then
    echo "Mappability track (for female+Y) not found!"

    #source gem
    PATH=${GEM_FOLDER}:$PATH

    #should the mappability be calculated using the female+Y or only Y chromosome
    if [ "${fmappability}" = "female" ]; then
        cat ${female} ${assembly_name}/${assembly_name}.continous.fa >${assembly_name}/female_genome_and_${assembly_name}.continous.fa #we need to concatenate female genome with Y chromosome for the mappability calculation
    else
        cat ${assembly_name}/${assembly_name}.continous.fa >${assembly_name}/female_genome_and_${assembly_name}.continous.fa #only Y assembly is used for the mappability calculation
    fi

    gem-indexer -i ${assembly_name}/female_genome_and_${assembly_name}.continous.fa -o ${assembly_name}/female_genome_and_${assembly_name}.continous.fa --complement emulate --verbose
    gem-mappability -I ${assembly_name}/female_genome_and_${assembly_name}.continous.fa.gem -l 101 -o ${assembly_name}/female_genome_and_${assembly_name}.continous.fa_MAPPABILITY -m 2 -e 2
    gem-2-wig -I ${assembly_name}/female_genome_and_${assembly_name}.continous.fa.gem -i ${assembly_name}/female_genome_and_${assembly_name}.continous.fa_MAPPABILITY.mappability -o ${assembly_name}/female_genome_and_${assembly_name}.continous.fa_MAPPABILITY
    wig2bed < ${assembly_name}/female_genome_and_${assembly_name}.continous.fa_MAPPABILITY.wig > ${assembly_name}/female_genome_and_${assembly_name}.continous.fa_MAPPABILITY.bed
fi


if [ ! -f ${assembly_name}/${assembly_name}.continous.fa.out ]; then
    echo "Repeatmasker has not been called before.!"
    RepeatMasker -species primates -pa 63 ${assembly_name}/${assembly_name}.continous.fa
fi

if [ ! -f ${assembly_name}/"chrY_"${assembly_name}"_REPEAT_MASKED.txt" ]; then
    echo "Repeatmasked file not found!"
    #skip header and reformat to .bed
    tail -n +4 ${assembly_name}/${assembly_name}.continous.fa.out | sed 's/^/ /' | tr -s " " | awk '{ print $5"\t"$6"\t"$7"\t"$10"\t"$1}' > ${assembly_name}/RM_${assembly_name}_chrY.bed
    python bin/analyzeRepeats.py ${assembly_name}/RM_${assembly_name}_chrY.bed ${assembly_name} ${reference_length} #outputs "chrY_"assembly_name"_REPEAT_MASKED.txt" 
fi

if [ ! -f ${assembly_name}/"chrY_MAPPABILITY_"${assembly_name}"_bypos.txt" ]; then
    echo "chrY_MAPPABILITY_"${assembly_name}"_bypos.txt not found!"
    #this should subset from female+Y to the Y chromosome only
    grep "^${assembly_name}" ${assembly_name}/female_genome_and_${assembly_name}.continous.fa_MAPPABILITY.bed >${assembly_name}/chrY_${assembly_name}_MAPPABILITY.bed
    python bin/generate_positional_information.py ${assembly_name}/chrY_${assembly_name}_MAPPABILITY.bed ${assembly_name} ${reference_length} & #generates "chrY_MAPPABILITY_"+args.assembly_name+"_bypos.txt"
fi

if [ ! -f ${assembly_name}/"chrY_GC_Content"${assembly_name}".txt" ]; then
    echo "chrY_GC_Content"${assembly_name}".txt not found!"
    python bin/analyzeGCcontent.py ${assembly_name}/${assembly_name}.continous.fa ${assembly_name} & #generates "chrY_GC_Content" + args.assembly_name + ".txt"
fi

#if [ ! -f trfMask_chrY.bed ]; then
#    echo "chrY_MAPPABILITY.bed not found!"
#    trf ${index} 2 7 7 80 10 50 2000 -l 6 -f -d -h -ngs >trfMask_chrY.bed &
#fi

wait
#AT THIS POINT, ALL NEEDED FILES SHOULD BE ALREADY CREATED
#WE CAN GENERATE THE SUMMARY FILE


if [ ! -f ${assembly_name}/"Summary_"${assembly_name}"_chrY_Map_RM_TRF_GC_MapC.tab" ]; then
    echo "Summary_"${assembly_name}"_chrY_Map_RM_TRF_GC_MapC.tab not found!"
    Rscript --vanilla bin/createSummary.R ${assembly_name} ${assembly_name}/"chrY_MAPPABILITY_"${assembly_name}"_bypos.txt" ${assembly_name}/"chrY_"${assembly_name}"_REPEAT_MASKED.txt" ${assembly_name}/"chrY_GC_Content"${assembly_name}".txt" ${reference_length} #generates "Summary_"${assembly_name}"_chrY_Map_RM_TRF_GC_MapC.tab"
    echo "Summary file generated."
fi


if [ ! -f ${bam_file}.bam"window_CN_summary.txt" ]; then
    #START COPY NUMBER ANALYSIS
    python bin/set_mappability_to_0_for_boundaries.py ${assembly_name}/Summary_${assembly_name}_chrY_Map_RM_TRF_GC_MapC.tab ${assembly_name}/${index}.offset.txt #set mappability to 0 for the scaffold boundaries
    #python CopyNumber_Yassembly.py --BAM ${bam_file}.bam --CHR ${chromosome_name} --GFile ${assembly_name}/${window_file} --SFile ${assembly_name}/Summary_${assembly_name}_chrY_Map_RM_TRF_GC_MapC.tab --READ PAIRED
    python CopyNumber_Yassembly.py --BAM ${bam_file}.bam --CHR ${chromosome_name} --GFile ${assembly_name}/${window_file} --SFile ${assembly_name}/Summary_${assembly_name}_chrY_Map_RM_TRF_GC_MapC.tab --READ SINGLE
    echo "Copy number calculated."
fi

if [ ! -f ${assembly_name}/${index}".MAPQ.txt" ]; then
    #GENERATE THE MAPPING FILE TO OBTAIN MAPQ VALUES FOR EACH WINDOW
    bedfile=${assembly_name}/${assembly_name}.windows.bed

    cp ${index} ${assembly_name}/${index}
    fastaFromBed -fi ${assembly_name}/${index} -bed ${bedfile} -fo ${bedfile}.fasta #create sequences from the intervals
    bwa index ${assembly_name}/${index}
    bwa mem ${assembly_name}/${index} ${bedfile}.fasta >${bedfile}.sam
    samtools view -bhS ${bedfile}.sam >${bedfile}.bam
    samtools sort --threads 63 ${bedfile}.bam >${bedfile}.sorted.bam
    mv ${bedfile}.sorted.bam ${bedfile}.bam
    samtools index ${bedfile}.bam

    echo -e "window\tflag\tcontig\tposition\tMAPQ\tCIGAR" >${assembly_name}/${index}.MAPQ.txt #add header
    samtools view ${bedfile}.bam | cut -f1-6 | sed 's/:/_/g' | sed 's/-/_/g' >>${assembly_name}/${index}.MAPQ.txt
    rm ${bedfile}.sam
    echo "Assembly windows mapped and MAPQ.txt generated."
fi

Rscript --vanilla bin/pickThreshold.R ${assembly_name}/${index}".MAPQ.txt" ${bam_file}".bamwindow_CN_summary.txt" ${directory}"/"${assembly_name}
mv Rplots.pdf ${assembly_name}.pickThreshold.pdf

echo "Pipeline finished. Done."
echo "#######################"
