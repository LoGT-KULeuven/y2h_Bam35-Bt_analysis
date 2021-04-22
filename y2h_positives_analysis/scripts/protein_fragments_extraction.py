import re
import sys

import pandas as pd

from Bio import SeqIO

NUM_THREADS = 8

blastp_infile = sys.argv[1]
clusters_infile = sys.argv[2]
fragments_infile = sys.argv[3]
fragments_outfile = sys.argv[4]

# organizing the cluster file as dictionary with
# key= representative sequence of cluster (read name)
# value = list of read names belonging to the cluster
pattern = r">(.*)\.\.\." # capture read name in clstr file²
cluster = []
results = {}
representative = ""
with open(clusters_infile, 'r') as clusters:
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

# filtering out fragments that did not blastp to proteome
results_filtered = {}
blastp = pd.read_csv(blastp_infile, sep="\t", header=0)
for i in range(0, len(blastp)):
    qacc = blastp.iloc[i]["qacc"]
    results_filtered[qacc] = results[qacc]

# concatenate all lists
fragments_out = []
for key in results_filtered:
    fragments_out += results_filtered[key]

print(len(fragments_out))
set_fragments = set(fragments_out)

# filtering fragments
fragments_in = list(SeqIO.parse(open(fragments_infile),'fasta'))
fragments_validated_out = []
for fragment in fragments_in:
    if(fragment.id in set_fragments):
        fragments_validated_out.append(fragment)
    
print(len(fragments_validated_out))
SeqIO.write(fragments_validated_out, fragments_outfile, "fasta")
#clusters_infile.close()
