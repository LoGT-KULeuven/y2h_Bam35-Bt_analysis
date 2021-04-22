import os
import sys
import re

import pandas as pd

from Bio import SeqIO

fragments_nr_infile = sys.argv[1]
faas = sys.argv[2]
blastp_file = sys.argv[3]
cluster_file = sys.argv[4]
sample = sys.argv[5]

ffns = list(SeqIO.parse(open(fragments_nr_infile, "rt"), 'fasta'))
proteome = list(SeqIO.parse(open(faas, "rt"), 'fasta'))
blastp = pd.read_csv(blastp_file, sep="\t", header=0)

pattern = r">(.*)\.\.\." # capture read name in clstr file
cluster = []
results = {}
representative = ""
with open(cluster_file, 'r') as clusters:
    for line in clusters:
        if(line.startswith(">")):
            if(len(cluster) != 0):
                results[representative] = cluster[:]
            cluster = []
            representative = ""
        else:
            m = re.search(pattern, line)
            cluster.append(m.group(1))
            if(line.endswith("*\n")):
                representative = m.group(1)
    results[representative] = cluster[:]

def lcs(S,T):
    m = len(S)
    n = len(T)
    counter = [[0]*(n+1) for x in range(m+1)]
    longest = 0
    lcs_set = set()
    for i in range(m):
        for j in range(n):
            if S[i] == T[j]:
                c = counter[i][j] + 1
                counter[i+1][j+1] = c
                if c > longest:
                    lcs_set = set()
                    longest = c
                    lcs_set.add(S[i-c+1:i+1])
                elif c == longest:
                    lcs_set.add(S[i-c+1:i+1])
    return lcs_set

count = 0

for fn in ffns:
    # print("Fragment %i: %s" % (count, fn.description))
    faa = fn.translate()
    faa_seq = str(faa.seq)
    # print(faa_seq)

    # find blastp result for the fragment
    i = blastp.loc[blastp['qacc']==fn.id].index
    #print(i)
    hit_blastp = blastp.loc[i, 'sacc'].values[0]
    #print("%i\t%shit_blastp)

    # find that hit in the list of protein (sorted)
    index = 0
    for idx, val in enumerate(proteome):
        if str(val.id) == str(hit_blastp):
            index = idx
            break

    # actually should be careful if index < 3, but that doesn't
    # happen in our dataset
    
    for fa in proteome[index-3:index+4]:
        fa_seq = str(fa.seq)
        lcss = lcs(faa_seq, fa_seq)
        for val in lcss:
            if(len(val)> 10):
                print("%i\t%s\t%s\t%s\t%s\t%i\t%i\t%s" % (count, sample, fa.id, fn.description, hit_blastp, len(results[fn.id]), len(fn.seq), faa_seq))
    count += 1
    #if(count > 10): break
