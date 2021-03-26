import gzip
import os
import sys

from Bio import SeqIO, bgzf

R1IN = sys.argv[1] 
R2IN = sys.argv[2]
R1OUTCLEAN = sys.argv[3]
R2OUTCLEAN = sys.argv[4]
R1OUTREMATG = sys.argv[5]
R2OUTREMATG = sys.argv[6]
R1OUTREMPROT = sys.argv[7]
R2OUTREMPROT = sys.argv[8]
PROTEOME_FAA = sys.argv[9]

R1_records = SeqIO.parse(gzip.open(R1IN, "rt"), "fastq")
R2_records = SeqIO.parse(gzip.open(R2IN, "rt"), "fastq")

R1_pass = gzip.open(R1OUTCLEAN, 'wt')
R2_pass = gzip.open(R2OUTCLEAN, 'wt')
R1_removed_nonATG = gzip.open(R1OUTREMATG, 'wt')
R2_removed_nonATG = gzip.open(R2OUTREMATG, 'wt')
R1_removed_nonProt = gzip.open(R1OUTREMPROT, 'wt')
R2_removed_nonProt = gzip.open(R2OUTREMPROT, 'wt')

faa_file = open(PROTEOME_FAA)
proteome = faa_file.read()

for (R1_reads, R2_reads) in zip(R1_records, R2_records):

    if not R1_reads.seq.startswith("ATG"):
        R1_removed_nonATG.write(R1_reads.format("fastq"))
        R2_removed_nonATG.write(R2_reads.format("fastq"))
        continue


    # other filter initially considered, but not pursued
    R1_translated = R1_reads.translate()

    if str(R1_translated.seq)[0:15] in proteome:
        R1_pass.write(R1_reads.format("fastq"))
        R2_pass.write(R2_reads.format("fastq"))
    else:
        R1_removed_nonProt.write(R1_reads.format("fastq"))
        R2_removed_nonProt.write(R2_reads.format("fastq"))
        
R1_pass.close()
R2_pass.close()
R1_removed_nonATG.close()
R2_removed_nonATG.close()
R1_removed_nonProt.close()
R2_removed_nonProt.close()
