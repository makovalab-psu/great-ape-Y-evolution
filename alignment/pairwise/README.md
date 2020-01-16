# pairwise alignment

All pairwise alignment relating to this section of the repository was done
using either lastz, or by extracting pairwise alignments from progressiveCactus
multiple alignments.

## contents

alignments_to_best.py&mdash;
From an alignment, collect intervals that are "best" as measured by some
column.

human_primate.q&mdash;
Scoring and parameters file used for lastz primate pairwise alignments.

## mini-pipelines

Record the length of scaffolds/contigs in species S1's assembly&mdash;

```bash  
cat ${S1}.smsk.fa | fasta_lengths > ${S1}.lengths
```

Record the runs of Ns in species S1's assembly&mdash;

```bash  
cat ${S1}.smsk.fa | fasta_N_intervals > ${S1}.N_intervals
```

Record the nucleotide content (counts) of scaffolds/contigs in species S1's
assembly&mdash;

```bash  
cat ${S1}.smsk.fa | fasta_content > ${S1}.content
```

Record the number of non-N bases in species S1's assembly&mdash;

```bash  
cat ${S1}.content \
  | pick_named_columns --sep=whitespace nonN \
  | awk '{ nonNTotal += $0 }
     END { printf("%d\n",nonNTotal) }' \
  > ${S1}.nonN
```

Convert species S1's assembly to 2bit format&mdash;

```bash  
faToTwoBit ${S1}.smsk.fa ${S1}.smsk.2bit
```

Align species S1 and S2&mdash;

```bash  
# perform alignment

jobId=${S1}_${S2}
lastz \
      ${S1}.smsk.2bit[multiple,unmask] \
      ${S2}.smsk.2bit[unmask] \
      --scores=human_primate.q \
      --seed=match12 \
      --allocate:traceback=200M \
      --format=general:name1,zstart1,end1,name2,strand2,zstart2+,end2+,score,cov%,con%,id%,text1,text2 \
      --markend \
  | gzip \
  > ${jobId}.dat.gz

# sort alignments by decreasing identity

gzip -dc ${jobId}.dat.gz \
  | tr -d \"%\" \
  | awk '/^#/  { if (NR==1) print $0 }
         !/^#/ { print $0 | "sort -k 1,1d -k 11,11nr" }' \
  | gzip \
  > ${jobId}.sorted.gz
```

Compute the amount of each species S1 aligned to S2&mdash;

```bash  
jobId=${S1}_${S2}
gzip -dc ${jobId}.sorted.gz \
  | awk '{ print $1,$2,$3 }' \
  | genodsp --novalue --uncovered:show \
      --chromosomes=${S1}_Y.lengths \
  | gzip \
  > ${jobId}.cov.dat.gz

nonN=`cat ${S1}_Y.nonN`
gzip -dc ${jobId}.cov.dat.gz \
  | awk '{ if ($4 != 0) bp += $3-$2 }
     END { printf("%s\t%d\t%d\t%0.2f%%\n",jobId,bp,nonN,100.0*bp/nonN) }' \
     jobId=${jobId} nonN=${nonN} \
  > ${jobId}.cov.bp.dat
```

Find the highest identity alignment between species S1 and S2 at each position
in S1&mdash;

At the same time, discard short alignments (fewer than 30 alignment columns).

```bash  
minAlignment=30    # discard alignments with fewer columns than this
text1Col=12

jobId=${S1}_${S2}
gzip -dc ${jobId}.sorted.gz \
  | awk '/^#/  { print $0 }
         !/^#/ { if (length($text1Col)>=M) print $0 }' \
               text1Col=${text1Col} \
               M=${minAlignment} \
  | alignments_to_best --format=automatic \
      --interval=sequence1 --key=id \
  | gzip \
  > ${jobId}.id.best.gz
```

Make a histogram of the distribution of identity between species S1 and
S2&mdash;

```bash  
# find intervals that have no alignment, excluding N intervals

jobId=${S1}_${S2}
gzip -dc ${jobId}.id.best.gz \
  | awk '{ print $1,$2,$3 }' \
  | genodsp --novalue --uncovered:hide --nooutputvalue \
      --chromosomes=${S1}_Y.lengths \
      = clip --max=1 \
      = invert 0.5 \
      = mask ${S1}_Y.N_intervals \
  | gzip \
  > temp.${jobId}.unaligned.gz

# compute histogram

gzip -dc ${jobId}.id.best.gz \
         temp.${jobId}.unaligned.gz \
  | grep -v "^#" \
  | awk '{ if (NF==3) print $1,$2,$3,"0.0";
                 else print $1,$2,$3,$idCol; }' idCol=${idCol} \
  | awk '{ bp[$4] += $3-$2 }
     END { for (id in bp) printf("%s %d\n",id,bp[id]) }' \
  | sort -nr \
  | awk 'BEGIN { print "id","bases" } { print $1,$2 }' \
  > ${jobId}.id.hist_bins

# cleanup

rm temp.${jobId}.unaligned.gz
```

Compute average identity between species S1 and S2&mdash;

```bash  
jobId=${S1}_${S2}
cat ${jobId}.id.hist_bins \
  | grep -v "^id" \
  | awk '{
         id = $1;  bases = $2;
         if (id != 0)
           { numer += id*bases;  denom += bases; }
         }
     END {
         printf("%s %s %.2f%% %d\n",
                S1,S2,numer/denom,denom);
         }' S1=${S1} S2=${S2} \
  | awk 'BEGIN { print "#S1 S2 id.avg bp.aligned" }
               { print $0; }' \
```

