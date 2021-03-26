import os
import sys

from Bio import SeqIO, bgzf

FRAGIN = sys.argv[1] 
FRAGOUTATG = sys.argv[2]
FRAGOUTNONATG = sys.argv[3]

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
