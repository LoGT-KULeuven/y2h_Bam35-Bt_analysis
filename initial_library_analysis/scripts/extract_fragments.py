import os
import sys

import pandas as pd

from Bio import SeqIO

intersect_infile = sys.argv[1]
genome_file = sys.argv[2]
fragments_outfile = sys.argv[3]
fragments_rejected = sys.argv[4]
fragments_atg = sys.argv[5]
fragments_non_atg = sys.argv[6]

intersect = pd.read_csv(intersect_infile, sep="\t", header=0)
genome = list(SeqIO.parse(open(genome_file),'fasta'))

accession = { "CP050183" : 0,
              "CP050184" : 1,
              "CP050185" : 2,
              "CP050186" : 3 }

start = 0
stop = 0
sequences = []
rejected = []
for i in range(0, len(intersect), 2):  #len(intersect)):
    if(intersect.iloc[i][5] == "+"):
        start = intersect.iloc[i][1]
        stop = intersect.iloc[i+1][2]
        sequence = genome[accession[intersect.iloc[i][0]]][start:stop]
        sequence.description = "pos"
    else:
        start = intersect.iloc[i+1][1]
        stop = intersect.iloc[i][2]
        sequence = genome[accession[intersect.iloc[i][0]]][start:stop]
        sequence = sequence.reverse_complement()
        sequence.description = "neg"
    sequence.id = intersect.iloc[i][3]

    if(len(sequence.seq) > 1000):
        rejected.append(sequence)
    else:
        sequences.append(sequence)

SeqIO.write(sequences, fragments_outfile, "fasta")
SeqIO.write(rejected, fragments_rejected, "fasta")

# adhoc trick to filter out fragment that don't start with ATG (due to mapping)
FRAGIN = fragments_outfile
FRAGOUTATG = fragments_atg
FRAGOUTNONATG = fragments_non_atg

FRAG_records = SeqIO.parse(open(FRAGIN, "rt"), "fasta")

FRAG_ATG = open(FRAGOUTATG, 'wt')
FRAG_nonATG = open(FRAGOUTNONATG, 'wt')


for FRAG in FRAG_records:
    if not FRAG.seq.startswith("ATG"):
        FRAG_nonATG.write(FRAG.format("fasta"))
    else:
        FRAG_ATG.write(FRAG.format("fasta"))
        
FRAG_ATG.close()
FRAG_nonATG.close()
