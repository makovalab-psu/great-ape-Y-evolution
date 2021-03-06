# multiple alignment

All multiple alignment relating to this section of the repository was done
using progressiveCactus.

## contents

fasta_file.py&mdash;
Support for reading FASTA files.

five${C}.unguided.seqFile (for C in 1,2,...,22,X,Y)&mdash;
Control file for progressiveCactus for the 5-species alignment of chromosome C.
No phylogenetic tree is specified, so progressiveCactus will assume a star
tree (a single root with all leaves connected to it).

fiveXY.guided.seqFile&mdash;
Control file for progressiveCactus for the 5-species alignment of chromosomes X
and Y. A phylogenetic tree is specified, a root with five branches, each
splitting into branches for X and Y.

maf_blocks_to_subset_base_counts.py&mdash;
Read alignments in maf format and, conceptually, partition blocks by the set of
species present. Collect and report species-specific stats within each subset.

maf_filter_by_species_set.py&mdash;
Read alignments in maf format and output those blocks that have a specified
set of species, and no other species.

maf_reader.py&mdash;
Support for reading MAF files.

maf_to_pairwise_identity.py&mdash;
Read alignments in maf format and output pairwise identity stats, and other
alignment stats.

multi_fasta_to_pairwise_identity.py&mdash;
Read an alignment in multi-fasta format and output pairwise identity stats.

subset_base_counts_plot.r&mdash;
R script to plot counts of aligned bases, by species subsets. It assumes an
input table as created by maf_blocks_to_subset_base_counts.py.

## mini-pipelines

Run progressiveCactus to create the 5-species alignment of chromosome Y&mdash;

```bash  
progressiveCactus \
  --maxThreads=24 \
  fiveY.unguided.seqFile \
  work_dir \
  fiveY.hal
```

Run progressiveCactus to create the 5-species alignment of chromosomes X and
Y&mdash;

```bash  
progressiveCactus \
  --maxThreads=24 \
  fiveXY.guided.seqFile \
  work_dir \
  fiveXY.hal
```

Run progressiveCactus to create the 5-species alignment of chromosome X&mdash;

```bash  
progressiveCactus \
  --maxThreads=24 \
  fiveX.unguided.seqFile \
  work_dir \
  fiveX.hal
```

Run progressiveCactus to create the 5-species alignment of chromosome C&mdash;

```bash  
progressiveCactus \
  --maxThreads=24 \
  five${C}.unguided.seqFile \
  work_dir \
  five${C}.hal
```

Convert the 5-species alignment to MAF format, with the putative ancestral
sequence Anc0 as the reference&mdash;

```bash  
hal2maf fiveY.hal \
      --refGenome Anc0 \
      --maxRefGap 100 \
      --maxBlockLen 10000 \
      /dev/stdout \
  | gzip \
  > fiveY.Anc0_centric.maf.gz
```

Convert the 5-species alignment to MAF format, with species S1 as the
reference&mdash;

```bash  
hal2maf fiveY.hal \
      --noAncestors --refGenome ${S1} \
      --refGenome Anc0 \
      --maxRefGap 100 \
      --maxBlockLen 10000 \
      /dev/stdout \
  | gzip \
  > fiveY.${S1}_centric.maf.gz
```

Collect by-subset stats from the 5-species alignment&mdash;

Note that the ancestor-centric MAF can lack species-specific portions of each
species, so the counts for singleton sets in the following mini-pipeline are
replaced by another mini-pipeline (described elsewhere in this README).

```bash  
gzip -dc fiveY.Anc0_centric.maf.gz \
  | maf_blocks_to_subset_base_counts \
      --species=hg_Y,panTro_Y,panPan_Y,gorGor_Y,ponAbe_Y \
  > fiveY.Anc0_centric.subset_base_counts
```

Identify and collect stats for species-specific segments for species S1, from
the 5-species alignment&mdash;

Note that the ancestor-centric MAF can lack species-specific portions of each
species, so the results of the following mini-pipeline replace the counts for
singleton sets in fiveY.Anc0_centric.subset_base_counts.

Species-specific segments are created as a byproduct of this pipeline, in BED
format.

```bash  
gzip -dc fiveY.${S1}_Y_centric.maf.gz \
  | maf_filter_by_species_set --species=${S1}_Y \
  | grep ${S1}_Y \
  | awk '{
         if ($5 == "+") print $2,$3,$3+$4;
         else           print $2,$6-($3+$4),$6-$3;
         }' \
  | sed "s/^${S1}_Y\.//" \
  | genodsp --novalue --uncovered:hide --nooutputvalue \
      --chromosomes=${S1}_Y.lengths \
      = mask ${S1}_Y.N_intervals \
  | tee fiveY.${S1}_Y_centric.${S1}_Y_specific.bed \
  | awk '{ t+=$3-$2; }
     END { printf("%s %d\n",species,t); }' species=${S1} \
  > fiveY.${S1}_Y_centric.${S1}_Y_specific

gzip -dc fiveY.${S1}_Y_centric.maf.gz \
  | maf_filter_by_species_set --species=${S1}_Y \
      --discard:duplicates \
  | grep ${S1}_Y \
  | awk '{
         if ($5 == "+") print $2,$3,$3+$4;
         else           print $2,$6-($3+$4),$6-$3;
         }' \
  | sed "s/^${S1}_Y\.//" \
  | genodsp --novalue --uncovered:hide --nooutputvalue \
      --chromosomes=${S1}_Y.lengths \
      = mask ${S1}_Y.N_intervals \
  | awk '{ t+=$3-$2; }
     END { printf("%s %d\n",species,t); }' species=${S1} \
  > fiveY.${S1}_Y_centric.${S1}_Y_unique

cat fiveY.*_Y_specific > temp.specific
cat fiveY.*_Y_unique   > temp.unique
merge_file_columns_by_common_name --removekey \
    temp.specific \
    temp.unique \
  | awk 'BEGIN { print "#species","specific","unique" } { print $0 }' \
  | colrithmetic --header "duplicate:specific-unique"
```

Extract alignment blocks that contain all five species&mdash;

```bash  
gzip -dc fiveY.Anc0_centric.maf.gz \
  | maf_filter_by_species_set \
      --species=hg_Y,panTro_Y,panPan_Y,gorGor_Y,ponAbe_Y \
  | gzip \
  > fiveY.Anc0_centric.all_five .maf.gz
```

Compute average identity between each pair of species, in alignment blocks that
contain all five species&mdash;

chrY.ancestor.pruned.5species.fas contains a combined alignment of the
five-species alignment blocks, in multi-fasta format. Its creation is described
elsewhere in this repository.

```bash  
cat chrY.ancestor.pruned.5species.fas \
  | multi_fasta_to_pairwise_identity
```

Compute the amount of each species aligned to each other species&mdash;

The output of maf_to_pairwise_identity is a table with one row for each
species pair. The "aligned" column is the amount of one species aligned to
the other.

```bash  
gzip -dc fiveY.Anc0_centric.maf.gz \
  | maf_to_pairwise_identity \
      hg_Y     --fasta:hg_Y=hg38.msY.smsk.fa \
      panTro_Y --fasta:panTro_Y=panTro6.msY.smsk.fa \
      panPan_Y --fasta:panPan_Y=panPan.msY.makovalab.ver1.smsk.fa \
      gorGor_Y --fasta:gorGor_Y=gorGor.msY.makovalab.ver3.smsk.fa \
      ponAbe_Y --fasta:ponAbe_Y=ponAbe.msY.makovalab.ver3.smsk.fa \
      --discard:weeds
```
