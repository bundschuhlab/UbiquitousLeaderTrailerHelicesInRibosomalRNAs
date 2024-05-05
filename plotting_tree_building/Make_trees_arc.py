import csv
import os
from ete3 import Tree


rRNA_type = "16S"

tree_fldr = "new-results-arc/trees_" + rRNA_type + "/"
os.makedirs(tree_fldr, exist_ok=True)
singleseq_fldr = tree_fldr + "singleseq/"
os.makedirs(singleseq_fldr, exist_ok=True)
rnaclust_fldr = tree_fldr + "rnaclust/"
os.makedirs(rnaclust_fldr, exist_ok=True)

clustered_results_loc = "new-results/" + rRNA_type + "-data.csv"

acc_fullnode_dict = {}
sp_clusters = "../leader_trailer/resources/temp/sp_clusters.tsv"
with open(sp_clusters, "r") as f:
    next(f)
    for row in f:
        data = row.split("\t")[0]
        
        if "GB_" in data:
            acc = data.split("GB_")[1]
        elif "RS_" in data:
            acc = data.split("RS_")[1]
        
        acc_fullnode_dict[acc] = data
f.close()

gtdb_full_tree_loc = "../codon_freqs/resources/ar53_r207.tree"
t = Tree(gtdb_full_tree_loc, format=1, quoted_node_names = True)

all_accs = set()
singleseq_acc_signal_dict = {}
rnaclust_acc_signal_dict = {}
singleseq_acc_taxonomy_dict = {}
rnaclust_acc_taxonomy_dict = {}
with open(clustered_results_loc, "r") as f:
    reader = csv.reader(f)
    next(reader, None)

    for row in reader:
        domain, species, acc, counter, \
            z_score, bio_signal, mean, stdev, \
            cluster, covariation, num_sws = row
        
        all_accs.add(acc)

        # Assign single-seq signals
        if z_score == "":
            singleseq_color = "#D3D3D3" # light grey
        else:
            z_score = float(z_score)
            if z_score <= 0.0:
                singleseq_color = "#00008B" # dark blue
            elif z_score <= 1.0:
                singleseq_color = "#ADD8E6" # light blue
            elif z_score <= 2.0:
                singleseq_color = "#FF9A98" # light red
            else:
                singleseq_color = "#FF0500" # dark red

        # Assign rnaclust signals
        if num_sws == "":
            rnaclust_color = "#D3D3D3" # light grey
        else:
            sw_num = float(num_sws)
            if covariation == "":
                rnaclust_color = "#D3D3D3" # light grey
                continue
            fract = 1.0 - float(covariation)

            if sw_num >= 15 and fract <= 0.75:
                rnaclust_color = "#FF0500" # dark red # STEM, Covariation
            elif sw_num >= 15 and fract > 0.75:
                rnaclust_color = "#FF9A98" # light red # STEM, No Covariation
            elif sw_num < 15 and fract <= 0.75:
                rnaclust_color = "#00008B" # dark blue # No Stem, Covariation
            elif sw_num < 15 and fract > 0.75:
                rnaclust_color = "#ADD8E6" # light blue # No Stem, No Covariation
            else:
                raise Exception("unable to assign rnaclust:", sw_num, fract, acc, species, counter)
        
        # Write to signal dict
        singleseq_acc_signal_dict[acc + "-" + counter] = singleseq_color
        rnaclust_acc_signal_dict[acc + "-" + counter] = rnaclust_color

        # capture taxonomy/labeling info
        if z_score == "":
            singleseq_sig_out = " [" + "NA" + "]"
        else:
            singleseq_sig_out = " [" + str(round(z_score, 2)) + "]"

        if covariation == "":
            rnaclust_sig_out = " [" + "NA" + "]"
        else:
            rnaclust_sig_out = " [stem-" + str(round(float(num_sws), 2)) + \
                " ,covariation-" + str(round(float(covariation), 2)) + "]"

        singleseq_acc_taxonomy_dict[acc + "-" + counter] = species + "-" + counter + \
                singleseq_sig_out
        rnaclust_acc_taxonomy_dict[acc + "-" + counter] = species + "-" + counter + \
                rnaclust_sig_out
f.close()

# search for each available phyla
clade_tree_search = t.search_nodes(name = "d__Archaea")[0] 
bacteria_node = clade_tree_search.detach()

# now prune the tree to only those nodes of interest
desired_nodes = set()
acc_node_map = {}
for acc in all_accs:
    node_name = acc_fullnode_dict[acc]

    acc_node_map[acc] = node_name
    desired_nodes.add(node_name)

all_phyla = {}
for node in bacteria_node.search_nodes():
    node_name = node.name

    if "p__" in node_name:
        all_phyla[node_name] = []

phyla_name_node_dict = {}
encountered_phyla = []
# track all other phyla, lump then together
phyla_name_node_dict["other"] = []
for phyla in all_phyla:
    if phyla not in encountered_phyla:
        phyla_name_node_dict["other"].append(phyla)

def make_trees(acc_signal_dict, phyla, phyla_tree_loc, acc_taxonomy_dict, add_rrna):
    """
    """
    if add_rrna == True:
        # add lower nodes for multiple rRNA copies per organism
        for guid in acc_signal_dict:
            acc, counter = guid.split("-")
            node = acc_node_map[acc]

            tree_search = downselected_node.search_nodes(name = node)
            if len(tree_search) == 0:
                continue
            
            tree_search = tree_search[0]
            tree_search.add_child(name = guid)

    downselected_node.write(format = 1, outfile = phyla_tree_loc + phyla + "-" "itol.tree")

    with open(phyla_tree_loc + "tree-labels.txt", "w") as f:
        f.write("TREE_COLORS" + "\n")
        f.write("SEPARATOR TAB" + "\n")
        f.write("DATA" + "\n")

        f.write("100.0_p__Halobacteriota" + "\t" + "range" + "\t" + "#FFFFFF" + "\t" + "Halobacteriota" + "\n")
        f.write("78.0_p__Thermoproteota" + "\t" + "range" + "\t" + "#FFFFFF" + "\t" + "Thermoproteota" + "\n")
        f.write("100.0_p__Methanobacteriota_B_ c__Thermococci" + "\t" + "range" + "\t" + "#FFFFFF" + "\t" + "Methanobacteriota_B" + "\n")
        f.write("100.0_p__Methanobacteriota_ c__Methanobacteria_ o__Methanobacteriales" + "\t" + "range" + "\t" + "#FFFFFF" + "\t" + "Methanobacteriota" + "\n")
        f.write("26.0_p__Methanobacteriota_A" + "\t" + "range" + "\t" + "#FFFFFF" + "\t" + "Methanobacteriota_A" + "\n")

        for acc in acc_signal_dict:
            node_name = acc
            color = acc_signal_dict[acc]

            f.write(node_name + "\tlabel\t" + str(color) + "\n")
            f.write(node_name + "\tbranch\t" + str(color) + "\tnormal\t1\n")
    f.close()

    with open(phyla_tree_loc + "tree-names.txt", "w") as f:
        f.write("LABELS" + "\n")
        f.write("SEPARATOR TAB" + "\n")
        f.write("DATA" + "\n")

        for acc in acc_signal_dict:
            taxonomy = acc_taxonomy_dict[acc]
            #node_name = acc_node_map[acc]
            node_name = acc
            f.write(node_name + "\t" + taxonomy + "\n")
    f.close()    

    return

# for each phyla, count up # of species underneath that we have data for

for root_phyla in phyla_name_node_dict:
    print("processing root phyla:", root_phyla, "phyla:", phyla_name_node_dict[root_phyla])
    phyla_of_interest = phyla_name_node_dict[root_phyla]

    clean_external_nodes = set()
    internal_labels = set()
        
    for phyla in phyla_of_interest:
        t = Tree(gtdb_full_tree_loc, format=1, quoted_node_names = True)
        clade_tree_search = t.search_nodes(name = phyla)[0] 
        downselected_node = clade_tree_search.detach()

        for node in downselected_node.search_nodes():
            if node.is_leaf() == False:
                if "__" in node.name:
                    # check if leaves in set
                    node_leaves = node.get_leaf_names()
                    for n in node_leaves:
                        if n in desired_nodes:
                            internal_labels.add(node.name)
            else:
                if node.name in desired_nodes:
                    clean_external_nodes.add(node.name)

    t = Tree(gtdb_full_tree_loc, format=1, quoted_node_names = True)
    clade_tree_search = t.search_nodes(name = "d__Archaea")[0]
    downselected_node = clade_tree_search.detach()

    downselected_node.prune(list(clean_external_nodes.union(internal_labels)), preserve_branch_length = True)

    # write out singleseq + rnaclust trees
    make_trees(singleseq_acc_signal_dict, root_phyla, singleseq_fldr, singleseq_acc_taxonomy_dict, True)
    make_trees(rnaclust_acc_signal_dict, root_phyla, rnaclust_fldr, rnaclust_acc_taxonomy_dict, False)
