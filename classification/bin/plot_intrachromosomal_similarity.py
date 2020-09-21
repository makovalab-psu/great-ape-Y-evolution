#!/usr/bin/env python3
import os, errno
import argparse
import re
parser = argparse.ArgumentParser()                                               

parser.add_argument("--file", "-f", type=str, required=True)
args = parser.parse_args()

filepath = args.file
filename = os.path.basename(filepath)
print(filename)

output = open('Y_map_similarity' + filename + '.txt','w')

window_percentage = {} 

with open(filepath) as fp:  
   line = fp.readline()
   while line:
       array=line.split("\t")
       if (len(array)==12):
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = array

        if qseqid in window_percentage:
            max_percentage=window_percentage[qseqid]
            window_percentage[str(qseqid)] = max(max_percentage,pident)
        else:
            header=re.split('[- :]',qseqid) #check if this is self-alignment
            header_start=int(header[1])
            header_end=int(header[2])
            if (header_start!=int(sstart)):
                if (header_end!=int(send)):
                    window_percentage[str(qseqid)]=pident #UNIQUE MATCH
        line = fp.readline()
print("Finished writing the identities for the plotting.")
for key in sorted(window_percentage.iterkeys()):
    header=re.split('[- :]',key)
    output.write(header[0] + '\t' + header[1] + '\t' + header[2] + '\t' + window_percentage[key] + '\n')
output.close()
