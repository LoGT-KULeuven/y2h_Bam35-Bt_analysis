import os
import sys

from Bio.Blast.Applications import NcbiblastpCommandline

from Bio import SeqIO

NUM_THREADS = 8
EVALUE = 1e-10

ffn_infile = sys.argv[1]
ffas_out = sys.argv[2]
blastp_results = sys.argv[3]
db = sys.argv[4]

ffns = list(SeqIO.parse(open(ffn_infile),'fasta'))
ffas = []

for ffn in ffns:
    ffa = ffn.translate()
    ffa.id = ffn.id
    ffa.description = ffn.description
    ffas.append(ffa)

#print(ffa)
SeqIO.write(ffas, ffas_out, "fasta")

blastp_cline = NcbiblastpCommandline(query=ffas_out,
                                     db=db,
                                     outfmt='"6 qacc sacc stitle qlen slen qstart qend length nident mismatch positive evalue"',
                                     max_target_seqs=1,
                                     num_threads=NUM_THREADS,
                                     evalue=EVALUE,
                                     out=blastp_results)
blastp_cline()

with open(blastp_results, 'r') as original: data = original.read()
headers = "qacc\tsacc\tstitle\tqlen\tslen\tqstart\tqend\tlength\tnident\tmismatch\tpositive\tevalue\n"
with open(blastp_results, 'w') as modified: modified.write(headers + data)
