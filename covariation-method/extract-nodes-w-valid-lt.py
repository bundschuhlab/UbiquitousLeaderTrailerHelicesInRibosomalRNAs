import os
import csv

all_seqids = []
rRNA = "16S"
target_dirs = [rRNA + "_RNAClust", "Ar_" + rRNA + "_RNAClust"]
out_dir = "check_seqclusters/" + rRNA + "-seq-clusters.txt"

for target_dir in target_dirs:
    for dir in os.listdir(target_dir):
        seq_clusters = target_dir + "/" + dir + "/mmseqs_cluster.tsv"

        with open(seq_clusters, "r") as f:
            for row in f:
                cluster, seqid = row.replace("\n","").split("\t")
                all_seqids.append(seqid)
        f.close()

with open(out_dir, "w") as f:
    #out = csv.writer(f)
    for seqid in all_seqids:
        f.write(seqid + "\n")
f.close()
