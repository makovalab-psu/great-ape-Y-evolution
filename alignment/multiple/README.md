# multiple alignment
Support for multiple alignment tasks.

All multiple alignment relating to this section of the repository was done
using progressiveCactus.

## contents

maf_blocks_to_subset_base_counts.py&mdash;
Read alignments in maf format and, conceptually, partition blocks by the set of
species present. Collect and report species-specific stats within each subset.

maf_to_pairwise_identity.py&mdash;
Read alignments in maf format and output pairwise identity stats, and other
alignment stats.

multi_fasta_to_pairwise_identity.py&mdash;
Read an alignment in multi-fasta format and output pairwise identity stats.

fiveY.seqFile&mdash;
Control file for progressiveCactus for the 5-species Y alignment. This
describes the location of the input sequence files. No phylogenetic tree is
specified, so progressiveCactus will assume a star-tree (a single root with all
leaves connected to it).

## mini-pipelines

Run progressiveCactus to create the 5-species Y alignment&mdash;

```bash  
progressiveCactus \
  --maxThreads=24 \
  fiveY.seqFile \
  work_dir \
  fiveY.hal
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

```bash  
gzip -dc fiveY.Anc0_centric.maf.gz \
  | maf_blocks_to_subset_base_counts \
      --species=hg_Y,panTro_Y,panPan_Y,gorGor_Y,ponAbe_Y \
  > fiveY.Anc0_centric.subset_base_counts
```

Compute average identity between each pair of species, in alignment blocks that
contain all five species&mdash;

chrY.ancestor.pruned.5species.fas contains a combined alignment of the
five-species alignment blocks, in multi-fasta format. Its creation is described
elsewhere in this repository.

```bash  
cat chrY.ancestor.pruned.5species.fas \
  | multi_fasta_to_pairwise_identity \
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
