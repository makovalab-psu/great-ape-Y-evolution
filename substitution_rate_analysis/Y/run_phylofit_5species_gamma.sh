#!/bin/bash
#!/bin/bash
#SBATCH --job-name=phyloFit_nrates
#SBATCH --output=phyloFit_nrates-%j.out
#SBATCH --error=phyloFit_nrates-%j.err
#SBATCH --mem=50G
#SBATCH --ntasks=1

set -x

conda activate /galaxy/home/biomonika/conda/mutations > /dev/null 2>&1

#REV
#Y CHROMOSOME PRUNED ANCESTOR
#/galaxy/home/biomonika/phast/bin/phyloFit -E -Z -D 42 --subst-mod REV --nrates 4 --tree "(((hg_Y,panTro_Y),gorGor_Y),ponAbe_Y)" --out-root chrY.EZD42.subst-modREV.nrates4.speciesHCGS.ancestor.pruned.5species.maf --msa-format MAF chrY.ancestor.pruned.5species.maf

#AUTOSOMES PRUNED ANCESTOR
#/galaxy/home/biomonika/phast/bin/phyloFit -E -Z -D 42 --subst-mod REV --nrates 4 --tree "(((hg,panTro),gorGor),ponAbe)" --out-root autosomes.EZD42.subst-modREV.nrates4.speciesHCGS.ancestor.pruned.5species.maf --msa-format MAF autosomes.ancestor.pruned.5species.maf


#REV
#Y CHROMOSOME PRUNED ANCESTOR
/galaxy/home/biomonika/phast/bin/phyloFit -E -Z --subst-mod REV --nrates 4 --tree "(((hg_Y,(panTro_Y,panPan_Y)),gorGor_Y),ponAbe_Y)" --out-root chrY.EZD42.subst-modREV.nrates4.speciesHCBGS.ancestor.pruned.5species.maf --msa-format MAF chrY.ancestor.pruned.5species.maf

#AUTOSOMES PRUNED ANCESTOR
#/galaxy/home/biomonika/phast/bin/phyloFit -E -Z -D 42 --subst-mod REV --nrates 4 --tree "(((hg,(panTro,panPan)),gorGor),ponAbe)" --out-root autosomes.EZD42.subst-modREV.nrates4.speciesHCBGS.ancestor.pruned.5species.maf --msa-format MAF autosomes.ancestor.pruned.5species.maf


echo ".mod calculated. Done."