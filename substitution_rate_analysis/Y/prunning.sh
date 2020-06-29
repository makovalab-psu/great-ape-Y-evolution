#!/bin/bash
#SBATCH --job-name=prunning
#SBATCH --output=prunning-%j.out
#SBATCH --error=prunning-%j.err
#SBATCH --mem=50G
#SBATCH --ntasks=1

source /opt/anaconda/etc/profile.d/conda.sh
conda activate /galaxy/home/biomonika/conda/GApy2.7 > /dev/null 2>&1

set -e
set -x

maf_file=$(basename -- "$1")
maf_file="${maf_file%.*}"

#GET THE DATA READY FOR PRUNING
echo "##maf version=1 scoring=none" >${maf_file}.header
cat ${maf_file}.header ${maf_file}.maf >>${maf_file}.header.maf

/galaxy/home/biomonika/MAF/mafTools/bin/mafDuplicateFilter --maf ${maf_file}.header.maf >${maf_file}.prunedNOTflipped.maf
python maf_flip_for_ref.py Anc0. < ${maf_file}.prunedNOTflipped.maf >${maf_file}.pruned.maf  #human-centric

rm ${maf_file}.header
rm ${maf_file}.header.maf
rm ${maf_file}.prunedNOTflipped.maf
