#!/usr/bin/bash

for file in gene_msablock_fasta/*.fasta; do 
    if [ -f "$file" ]; then 
		echo "$file" 
		geneconv $file  /w9  /lp -nolog
    fi 
done
