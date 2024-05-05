import csv
from ete3 import Tree
import os


rRNA_type = "16S"
domain = "arc"

acc_count_dict = {}
acc_totcount_dict = {}
acc_species_dict = {}
with open(rRNA_type + "-species-uniq-cluster-analysis-data.csv", "r") as f:
    reader = csv.reader(f)
    next(reader, None)

    for row in reader:
        #acc, num_lts, num_uniq_clusters, num_predicted_yes, num_predicted_no = row
        acc, num_lts, num_predicted_yes, num_predicted_no, num_covariation_lts, num_uniq_clusters = row

        if int(num_covariation_lts) == 0:
            continue

        acc_count_dict[acc] = int(num_uniq_clusters)
        acc_totcount_dict[acc] = int(num_covariation_lts)
f.close()

with open("clusters-num-" + rRNA_type + ".csv", "r") as f:
    reader = csv.reader(f)
    for row in reader:
        acc, species = row[:2]
        acc_species_dict[acc] = species
f.close()

tree_fldr = domain + "_" + "trees_" + rRNA_type + "/"
os.makedirs(tree_fldr, exist_ok=True)

acc_fullnode_dict = {}
sp_clusters = "../../leader_trailer/resources/temp/sp_clusters.tsv"
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

if domain == "bac":
    gtdb_full_tree_loc = "../../codon_freqs/resources/bac120.tree"
elif domain == "arc":
    gtdb_full_tree_loc = "../../codon_freqs/resources/ar53_r207.tree"
t = Tree(gtdb_full_tree_loc, format=1, quoted_node_names = True)

if domain == "bac":
    clade_tree_search = t.search_nodes(name = "d__Bacteria")[0] 
elif domain == "arc":
    clade_tree_search = t.search_nodes(name = "d__Archaea")[0]
domain_node = clade_tree_search.detach()

# now prune the tree to only those nodes of interest
desired_nodes = set()
acc_node_map = {}
for acc in acc_count_dict:
    node_name = acc_fullnode_dict[acc]

    acc_node_map[acc] = node_name
    desired_nodes.add(node_name)

internal_labels = set()
clean_external_nodes = set()
for node in domain_node.search_nodes():
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

domain_node.prune(list(clean_external_nodes.union(internal_labels)), preserve_branch_length = True)

domain_node.write(format = 1, outfile = tree_fldr + "itol.tree")

# create a color gradient from min to max observed z-score
from colour import Color
light_red = Color("#e3cecb") 
dark_red = Color("#F91607")

with open(tree_fldr + "tree-labels.txt", "w") as f:
    f.write("TREE_COLORS" + "\n")
    f.write("SEPARATOR TAB" + "\n")
    f.write("DATA" + "\n")

    for acc in acc_count_dict:
        node_name = acc_node_map[acc]
        acc_signal = acc_count_dict[acc]
        fract = acc_signal / acc_totcount_dict[acc]

        # if only 1 total LT, make it grey
        if acc_totcount_dict[acc] == 1:
            color = "#808080" # grey

        elif fract <= 0.25 or acc_signal == 1:
            color = "#0000FF" # blue
        elif fract <= 0.5:
            color = "#FF2400" # light red
        else:
            color = "#800000" # dark red

        """
        if acc_totcount_dict[acc] == 1:
            color = "#808080" # grey

        elif acc_totcount_dict[acc] < 5:
            color = "#FF2400" # light red
        else:
            color = "#800000" # dark red
        if acc_totcount_dict[acc] == 1:
            color = "#808080" # grey
        elif acc_signal == 1:
            color = "#0000FF" # blue
        elif acc_signal < 5:
            color = "#FF2400" # light red
        else:
            color = "#800000" # dark red
        """

        """
        #color = colors[int(100 * (acc_signal - min_signal) / (max_signal - min_signal))]
        if acc_signal == "Unknown":
            color = "#808080" # Grey
        elif acc_signal == "TRUE":
            color = "#f72105" # red
        elif acc_signal == "FALSE":
            color = "#154c79" # dark blue
        """

        f.write(node_name + "\tlabel\t" + str(color) + "\n")
        f.write(node_name + "\tbranch\t" + str(color) + "\tnormal\t1\n")
f.close()

with open(tree_fldr + "tree-names.txt", "w") as f:
    f.write("LABELS" + "\n")
    f.write("SEPARATOR TAB" + "\n")
    f.write("DATA" + "\n")

    for acc in acc_count_dict:
        taxonomy = acc_species_dict[acc]
        node_name = acc_node_map[acc]
        f.write(node_name + "\t" + taxonomy + "\n")
f.close()    