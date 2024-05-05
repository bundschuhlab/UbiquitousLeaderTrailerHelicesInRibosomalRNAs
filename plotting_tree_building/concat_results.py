import csv

"""
Bacteria:
For rRNA: 16S
# of Seqs with EITHER Predict: 17625
# of Seqs with BOTH Predict: 14666
# of Seqs with RNAClust NOT Single: 98
# of Seqs with Single NOT RNAClust: 2861

For rRNA: 23S
# of Seqs with EITHER Predict: 17488
# of Seqs with BOTH Predict: 14485
# of Seqs with RNAClust NOT Single: 136
# of Seqs with Single NOT RNAClust: 2867

Archaea:
For rRNA: 16S
# of Seqs with EITHER Predict: 415
# of Seqs with BOTH Predict: 338
# of Seqs with RNAClust NOT Single: 0
# of Seqs with Single NOT RNAClust: 77

For rRNA: 23S
# of Seqs with EITHER Predict: 420
# of Seqs with BOTH Predict: 341
# of Seqs with RNAClust NOT Single: 2
# of Seqs with Single NOT RNAClust: 77
"""


rRNA_type = "16S"

RNAClust_results = "../LT_new_approach/results/" + rRNA_type + "-summary.csv"
vsScramb_results = "../leader_trailer_v2/results_" + rRNA_type + "-p__proteobacteria/" + \
                    rRNA_type + "-p__proteobacteria-bio_vs_scramb_structs-10_8nt.csv"
#vsScramb_results = "../leader_trailer_v2/results_" + rRNA_type + "-d__archaea/" + \
#                    rRNA_type + "-d__archaea-bio_vs_scramb_structs-10_8nt.csv"

out_loc = "new_results/" + rRNA_type + "-concat.csv"

import os
os.makedirs("new_results/", exist_ok=True)

# Load RNAClust results for every seq (go cluster -> members)
adj_species_genome_dict = {}
acc_species_dict = {}
assembly_data_summary = "../LT_new_approach/resources/taxonomy-genome-summary.txt"
with open(assembly_data_summary, "r") as f:
    for row in f:
        row = row.replace("\n","").split("\t")
        species, accs = row

        accs = accs.split(";")
        acc = accs[0]
        adj_species = species.replace(" ","_").lower()

        adj_species_genome_dict[adj_species] = acc
        acc_species_dict[acc] = species
f.close()

rnaclust_acc_dict = {}
with open(RNAClust_results, "r") as f:
    reader = csv.reader(f)
    next(reader, None)

    for row in reader:
        clade, cluster, num_seqs, fract_base_pairs, \
                                fract_red, fract_lt_red, seq_len, \
                                number_leader_trailer_bps, \
                                num_sws_meeting_threshold_10_8, \
                                num_sws_meeting_threshold_9_7, \
                                num_sws_meeting_threshold_8_6 = row
        
        rnaclust_fldr = "../LT_new_approach/" + rRNA_type + "_RNAClust/"
        #rnaclust_fldr = "../LT_new_approach/Ar_" + rRNA_type + "_RNAClust/"
        root_fldr = rnaclust_fldr + clade + "/"
        cluster_dir = root_fldr + cluster + "/"

        with open(cluster_dir + cluster + "__" + "cluster-seqs.out", "r") as ff:
            for row in ff:
                seqid = row.replace("\n","")
                prior = seqid.split("_")[0]
                corrected_seqid = seqid.split(prior + "_")[1]
                
                counter = seqid.split("-")[-1]
                taxid = seqid.split(prior + "_")[1][:-(1+len(counter))]

                if taxid.lower() not in adj_species_genome_dict:
                    print("unable to ID", seqid, taxid)
                    continue
                acc = adj_species_genome_dict[taxid.lower()]

                rnaclust_acc_dict[acc + "-" + counter] = [clade, cluster, \
                                                fract_lt_red, num_sws_meeting_threshold_10_8]
        ff.close()
f.close()

# Load single seq vs scrambled results
vsScramb_acc_dict = {}
with open(vsScramb_results, "r") as f:
    reader = csv.reader(f)
    next(reader, None)
    for row in reader:
        acc, counter, z_score, bio_signal, mean, stdev, leader_len, trailer_len, joiner_len, gc_fract = row

        vsScramb_acc_dict[acc + "-" + counter] = [z_score, bio_signal, mean, stdev, leader_len, trailer_len]
f.close()

# Write out combined results
rnaclust_keys = set(rnaclust_acc_dict.keys())
vsScramb_keys = set(vsScramb_acc_dict.keys())
combined_keys = rnaclust_keys | vsScramb_keys

print("For rRNA:", rRNA_type)
print("# of Seqs with EITHER Predict:", len(combined_keys))
print("# of Seqs with BOTH Predict:", len(rnaclust_keys.intersection(vsScramb_keys)))
print("# of Seqs with RNAClust NOT Single:", len(rnaclust_keys - vsScramb_keys))
print("# of Seqs with Single NOT RNAClust:", len(vsScramb_keys - rnaclust_keys))


with open(out_loc, "w") as f:
    out = csv.writer(f)
    out.writerow(["Acc", "Species", "Counter", \
                    "RNACLUST: Clade", \
                    "RNACLUST: Cluster ID", \
                    "RNACLUST: Fract LT Red", \
                    "RNACLUST: Num SWs 8_6", \
                    "SINGLE: Z Score", \
                    "SINGLE: Bio Signal", \
                    "SINGLE: Mean Scrambled", \
                    "SINGLE: StDev Scrambled", \
                    "SINGLE: Leader Length", \
                    "SINGLE: Trailer Length"])

    for guid in combined_keys:
        acc, counter = guid.split("-")
        species = acc_species_dict[acc]
        
        if guid in rnaclust_acc_dict:
            clade, cluster, \
                fract_lt_red, num_sws_meeting_threshold_8_6 = rnaclust_acc_dict[guid]
        else:
            clade, cluster, \
                fract_lt_red, num_sws_meeting_threshold_8_6 = "", "", "", ""
        
        if guid in vsScramb_acc_dict:
            z_score, bio_signal, mean, stdev, leader_len, trailer_len = vsScramb_acc_dict[guid]
        else:
            z_score, bio_signal, mean, stdev, leader_len, trailer_len = "", "", "", "", "", ""
        
        out.writerow([acc, species, counter, \
                    clade, cluster, \
                    fract_lt_red, num_sws_meeting_threshold_8_6, \
                    z_score, bio_signal, mean, stdev, \
                    leader_len, trailer_len])
f.close()
