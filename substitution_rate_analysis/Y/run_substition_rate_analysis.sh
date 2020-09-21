#!/bin/bash
#!/bin/bash
#SBATCH --job-name=run_substition_rate_analysis
#SBATCH --output=run_substition_rate_analysis-%j.out
#SBATCH --error=run_substition_rate_analysis-%j.err
#SBATCH --mem=150G
#SBATCH --ntasks=1

set -x

source /opt/anaconda/etc/profile.d/conda.sh > /dev/null 2>&1

srun --nodes=1 --ntasks=1 prunning.sh chrY.ancestor.maf #remove the duplicate regions from the MAF file
wait

conda activate /galaxy/home/biomonika/conda/GApy2.7 > /dev/null 2>&1
srun --nodes=1 --ntasks=1 python parse_cactus_Y.py #only keep parts of the alignment where all 5 species are present
wait
conda activate /galaxy/home/biomonika/conda/mutations > /dev/null 2>&1
srun run_phylofit_5species_gamma.sh #run phylofit to estimate the substitution rates
wait
echo "Done."