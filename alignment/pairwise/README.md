# pairwise alignment
Support for pairwise alignment tasks.

All pairwise alignment relating to this section of the repository was done
using either lastz, or by extracting pairwise alignments from progressiveCactus
multiple alignments.

# Contents

human_primate.q&mdash;
scoring and parameters file used for lastz primate pairwise alignments.

# Mini-pipelines


TBD add the alignment steps
TBD add the selection of best alignments per base

Make a histogram of the distribution of identity between species S1 and
S2&mdash;

```bash  
jobId=${S1}_${S2}

# find intervals that have no alignment, excluding N intervals

gzip -dc alignments/${jobId}.id.best.gz \
  | awk '{ print $1,$2,$3 }' \
  | genodsp --novalue --uncovered:hide --nooutputvalue \
      --chromosomes=${S1}_Y.lengths \
      = clip --max=1 \
      = invert 0.5 \
      = mask ${S1}_Y.N_intervals \
  | gzip \
  > temp.${jobId}.unaligned.gz

# compute histogram

gzip -dc alignments/${jobId}.id.best.gz \
         temp.${jobId}.unaligned.gz \
  | grep -v "^#" \
  | awk '{ if (NF==3) print $1,$2,$3,"0.0";
                 else print $1,$2,$3,$idCol; }' idCol=${idCol} \
  | awk '{ bp[$4] += $3-$2 }
     END { for (id in bp) printf("%s %d\n",id,bp[id]) }' \
  | sort -nr \
  | awk 'BEGIN { print "id","bases" } { print $1,$2 }' \
  > alignments/${jobId}.id.hist_bins
```

Compute average identity between species S1 and S2&mdash;

```bash  
jobId=${S1}_${S2}
cat alignments/${jobId}.id.hist_bins \
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

